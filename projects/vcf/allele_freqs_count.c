/* vcfallfreq.c - count the allele freq of VCF/BCF
 * SHI Quan (shiquan@genomics.cn)
 */
//#include <stdio.h>
#include "utils.h"
#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include "version.h"

enum filetype {
    filetype_promote_int = -1,
    file_is_vcf,
    file_is_bcf,
    file_is_compressed_bcf,
    file_is_bgzipped_vcf,
};

struct args {
    const char *fname_input;
    const char *fname_output;
    const char *ac_tag_string;
    const char *af_tag_string;
    const char *asam_tag_string;
    int output_filetype;
    htsFile *fp_input;
    htsFile *fp_output;
    bcf_hdr_t *hdr;
    int force_skip_uncover;
    int ntmp_arr;
    int32_t *tmp_arr;
    int n_cnt, m_cnt;
    int32_t n_samples;
    uint64_t *counts;
} args = {
    .fname_input = NULL,
    .fname_output = NULL,
    .ac_tag_string = NULL,
    .af_tag_string = NULL,
    .asam_tag_string = NULL,
    .output_filetype = file_is_vcf,
    .fp_input = NULL,
    .fp_output = NULL,
    .hdr = NULL,
    .force_skip_uncover = 0,
    .ntmp_arr = 0,
    .tmp_arr = NULL,
    .n_cnt = 0,
    .m_cnt = 0,
    .n_samples = 0,
    .counts = NULL,
};

static int smp_list = 0;

static const char *bcf_wmode(enum filetype file_type)
{
    if ( file_type == file_is_bcf )
        return "wbu";
    if ( file_type == file_is_compressed_bcf )
        return "wb";
    if ( file_type == file_is_bgzipped_vcf )
        return "wz";
    return "w";
}

int usage()
{
    fprintf(stderr, "A program to calculate the allele frequency in population for each position.\n"
            "Usage:\n"
            " vcf_allele_freq [options] input.vcf.gz\n"
            " Options:\n" 
            "     -ac   allele count tag name, default is AlleleCount\n"
            "     -af   allele frequency tag name, default is AlleleFreq\n"
            "     -hwe  tag name for Fisher exact test of Hardy-Weinberg Equilibrium\n"
            "     -sam  affected sample list tag name for each allele\n"
            "     -gen  sample counts for homozygous ref, heterozygous and homozygous alts, the format is [AA|Aa|aa]\n"
            "     -O    output format [b|z]\n"
            "     -o    output file\n"
            "     -force_skip_uncover  force skip the uncovered positions, default treat as reference\n"
            "\n"
            "Note:\n"
            "1. For parameter -gen, the tag format is [AA|Aa|aa], 'A' stand for wild allele, and 'a' stand for mutated\n"
            "   allele(s). All alternative alleles, if more than one alt allele like mnps, should be count into 'a'.\n"
            "Version: %s + htslib %s\n"
            "Homepage: https://github.com/shiquan/small_projects\n",
            PROJECTS_VERSION, hts_version()
        );
    return 1;
}

float hardy_weinberg_equilibrium_exact(int wild_hom, int het, int mutated_hom)
{
    
    return 0.0;
}
int parse_args(int argc, char **argv)
{
    int i;
    const char * output_filetype_string = 0;
    const char *version_string = 0;
    for ( i = 1; i < argc; ) {
        const char *a = argv[i++];
        if ( strcmp(a, "-h") == 0 || strcmp(a, "-help") == 0 )
            return usage();

        if ( strcmp(a, "-force_skip_uncover") == 0 ) {
            args.force_skip_uncover = 1;
            continue;
        }
        
        const char **var = 0;
        if ( strcmp(a, "-ac") == 0 )
            var = &args.ac_tag_string;
        else if ( strcmp(a, "-af") == 0 )
            var = &args.af_tag_string;
        else if ( strcmp(a, "-O") == 0 )
            var = &output_filetype_string;
        else if ( strcmp(a, "-o") == 0 )
            var = &args.fname_output;
        else if ( strcmp(a, "-sam") == 0 ) {
            var = &args.asam_tag_string;
            smp_list = 1;
        }
        else if ( strcmp(a, "-v") == 0 )
            var = &version_string;
        
        if ( var != 0 ) {
            if ( i == argc )
                error("Missing an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if ( args.fname_input == 0 ) {
            args.fname_input = a;
            continue;
        }

        error("Unknown argument : %s, use -h see help information.", a);
    }

    if ( version_string == 0 )
        error("Version of this database must be set by -v.");
    
    if ( args.fname_input == 0 && (!isatty(fileno(stdin))) )
        args.fname_input = "-";

    if ( args.fname_input == 0 )
        error("No input file.");

    args.fp_input = hts_open(args.fname_input, "r");
    if ( args.fp_input == 0 && (!isatty(fileno(stdin))) )
        args.fname_input = "-";
    
    if ( args.fp_input == NULL )
        error("Failed to open %s : %s.", args.fname_input, strerror(errno));
    
    htsFormat type = *hts_get_format(args.fp_input);
    if ( type.format != vcf && type.format != bcf )
        error("Unsupported input format.");

    args.output_filetype = file_is_vcf;
    if ( output_filetype_string != NULL ) {
        switch ( output_filetype_string[0] ) {
            case 'b':
                args.output_filetype = file_is_bcf;
                break;

            case 'u':
                args.output_filetype = file_is_vcf;
                break;

            case 'z':
                args.output_filetype = file_is_bgzipped_vcf;
                break;

            case 'v':
                args.output_filetype = file_is_vcf;
                break;

            default:
                error("The output type \"%c\" not recignised.", output_filetype_string[0]);
        }
    }

    args.fp_output = args.fname_output == 0 ? hts_open("-", bcf_wmode(args.output_filetype))
        : hts_open(args.fname_output, bcf_wmode(args.output_filetype));

    args.hdr = bcf_hdr_read(args.fp_input);

    args.n_samples = bcf_hdr_nsamples(args.hdr);

    if ( args.n_samples == 0 )
        error("No sample found.");
    
    if ( args.n_samples == 1 )
        error("Required more samples to calculate the allele frequency in population. Try to use bcftools merge more samples first.");
    
    if (args.ac_tag_string == NULL )
        args.ac_tag_string = "AlleleCount";

    if ( args.af_tag_string == NULL )
        args.af_tag_string = "AlleleFreq";

    // Update header IDs
    int id;
    id = bcf_hdr_id2int(args.hdr, BCF_DT_ID, args.ac_tag_string);
    if ( id == -1 ) {
        kstring_t str = { 0, 0, 0 };
        ksprintf(&str, "##INFO=<ID=%s,Number=A,Type=Integer,Description=\"Counts of each allele for all samples. Version=%s.\">", args.ac_tag_string, version_string);
        bcf_hdr_append(args.hdr, str.s);
        bcf_hdr_sync(args.hdr);
        id = bcf_hdr_id2int(args.hdr, BCF_DT_ID, args.ac_tag_string);
        assert(bcf_hdr_idinfo_exists(args.hdr, BCF_HL_INFO, id));
    }

    id = bcf_hdr_id2int(args.hdr, BCF_DT_ID, args.af_tag_string);
    if ( id == -1 ) {
        kstring_t str = { 0, 0, 0 };
        ksprintf(&str, "##INFO=<ID=%s,Number=A,Type=Float,Description=\"Frequency of each allele for all samples. Version=%s.\">", args.af_tag_string, version_string);
        bcf_hdr_append(args.hdr, str.s);
        bcf_hdr_sync(args.hdr);
        id = bcf_hdr_id2int(args.hdr, BCF_DT_ID, args.af_tag_string);
        assert(bcf_hdr_idinfo_exists(args.hdr, BCF_HL_INFO, id));
    }
    if ( smp_list == 1 ) {
        id = bcf_hdr_id2int(args.hdr, BCF_DT_ID, args.asam_tag_string);
        if ( id == -1 ) {
            kstring_t str = { 0, 0, 0 };
            ksprintf(&str, "##INFO=<ID=%s,Number=A,Type=Float,Description=\"Allele affected samples list.\">", args.asam_tag_string);
            bcf_hdr_append(args.hdr, str.s);
            bcf_hdr_sync(args.hdr);
            id = bcf_hdr_id2int(args.hdr, BCF_DT_ID, args.af_tag_string);
            assert(bcf_hdr_idinfo_exists(args.hdr, BCF_HL_INFO, id));
        }
    }
    bcf_hdr_write(args.fp_output, args.hdr);
    return 0;
}

int generate_freq(bcf1_t *line)
{
    uint64_t all_counts = 0;
    int i;
    
    args.n_cnt = line->n_allele;

    if ( args.m_cnt <= args.n_cnt ) {
        args.m_cnt = args.n_cnt + 1;
        args.counts = (uint64_t*)realloc(args.counts, args.m_cnt* sizeof(uint64_t));
    }
    for ( i = 0; i < args.n_cnt; ++i )
        args.counts[i] = 0;

    int ngt = bcf_get_genotypes(args.hdr, line, &args.tmp_arr, &args.ntmp_arr);

    // GT not present
    if ( ngt <= 0 )
        return 1;
    // if not diploid
    if ( ngt != args.n_samples*2)
        return 2;

    ngt /= args.n_samples;
    
    for ( i = 0; i < args.n_samples; ++i ) {
        int32_t *b = args.tmp_arr + i*ngt;
        int j;
        for ( j = 0; j < 2; ++j ) {
            if ( (bcf_gt_is_missing(b[j]) || b[j] == bcf_int32_vector_end) ) {
                if ( args.force_skip_uncover == 0 ) 
                    all_counts++;                
                continue;
            }
            if ( bcf_gt_allele(b[j]) == 0 ) {
                all_counts++;
                continue;
            }
            args.counts[bcf_gt_allele(b[j])-1] ++;
            all_counts++;
        }        
    }

    kstring_t str = { 0, 0, 0 };
    for ( i = 0; i < args.n_cnt-1; ++i ) {
        if ( i )
            kputc(',', &str);
        kputw(args.counts[i], &str);
    }
    bcf_update_info_string(args.hdr, line, args.ac_tag_string, str.s);
    str.l = 0;
    for ( i = 0; i < args.n_cnt-1; ++i ) {
        if ( i )
            kputc(',', &str);
        ksprintf(&str, "%f", (float)args.counts[i]/all_counts);        
    }
    bcf_update_info_string(args.hdr, line, args.af_tag_string, str.s);
    if ( str.m )
        free(str.s);
    return 0;
}

int calculate_frequence()
{
    bcf1_t *line = bcf_init();
    while ( bcf_read(args.fp_input, args.hdr, line) == 0 ) {
        if ( line->rid == -1 )
            continue;
        // calculation...
        generate_freq(line);

        // write
        bcf_write1(args.fp_output, args.hdr, line);
    }
    bcf_destroy(line);    
    return 0;
}

int release_memory()
{
    bcf_hdr_destroy(args.hdr);
    hts_close(args.fp_input);
    bcf_close(args.fp_output);
    if ( args.n_cnt )
        free(args.counts);
    if ( args.ntmp_arr )
        free(args.tmp_arr);
    return 0;
}


int main (int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    calculate_frequence();

    release_memory();
    
    return 0;
}


