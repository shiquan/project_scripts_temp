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
#include <math.h>

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
    const char *hwe_tag_string;
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
    .hwe_tag_string = NULL,
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
            "     -sam  affected sample list for each allele, sample seperated by /, format of this tag is [ Aa | aa ]\n"
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

#define SMALL_EPSILON 0.00000000000005684341886080801486968994140625
#define EXACT_TEST_BIAS 0.00000000000000000000000010339757656912845935892608650874535669572651386260986328125

// the SNPHWE2 function copy from plink_stats.c
double SNPHWE2(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, uint32_t midp) {
  // This function implements an exact SNP test of Hardy-Weinberg
  // Equilibrium as described in Wigginton, JE, Cutler, DJ, and
  // Abecasis, GR (2005) A Note on Exact Tests of Hardy-Weinberg
  // Equilibrium. American Journal of Human Genetics. 76: 000 - 000.
  //
  // The original version was written by Jan Wigginton.
  //
  // This version was written by Christopher Chang.  It contains the following
  // improvements over the original SNPHWE():
  // - Proper handling of >64k genotypes.  Previously, there was a potential
  //   integer overflow.
  // - Detection and efficient handling of floating point overflow and
  //   underflow.  E.g. instead of summing a tail all the way down, the loop
  //   stops once the latest increment underflows the partial sum's 53-bit
  //   precision; this results in a large speedup when max heterozygote count
  //   >1k.
  // - No malloc() call.  It's only necessary to keep track of a few partial
  //   sums.
  // - Support for the mid-p variant of this test.  See Graffelman J, Moreno V
  //   (2013) The mid p-value in exact tests for Hardy-Weinberg equilibrium.
  //
  // Note that the SNPHWE_t() function below is a lot more efficient for
  // testing against a p-value inclusion threshold.  SNPHWE2() should only be
  // used if you need the actual p-value.
    long obs_homc;
    long obs_homr;
    if (obs_hom1 < obs_hom2) {
        obs_homc = obs_hom2;
        obs_homr = obs_hom1;
    } else {
        obs_homc = obs_hom1;
        obs_homr = obs_hom2;
    }
    int64_t rare_copies = 2LL * obs_homr + obs_hets;
    int64_t genotypes2 = (obs_hets + obs_homc + obs_homr) * 2LL;
    int32_t tie_ct = 1;
    double curr_hets_t2 = obs_hets;
    double curr_homr_t2 = obs_homr;
    double curr_homc_t2 = obs_homc;
    double tailp = (1 - SMALL_EPSILON) * EXACT_TEST_BIAS;
    double centerp = 0;
    double lastp2 = tailp;
    double lastp1 = tailp;
    double curr_hets_t1;
    double curr_homr_t1;
    double curr_homc_t1;
    double preaddp;
    if (!genotypes2) {
        if (midp) {
            return 0.5;
        } else {
            return 1;
        }
    }
    
    if (obs_hets * genotypes2 > rare_copies * (genotypes2 - rare_copies)) {
        // tail 1 = upper
        while (curr_hets_t2 > 1.5) {
            // het_probs[curr_hets] = 1
            // het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
            curr_homr_t2 += 1;
            curr_homc_t2 += 1;
            lastp2 *= (curr_hets_t2 * (curr_hets_t2 - 1)) / (4 * curr_homr_t2 * curr_homc_t2);
            curr_hets_t2 -= 2;
            if (lastp2 < EXACT_TEST_BIAS) {
                if (lastp2 > (1 - 2 * SMALL_EPSILON) * EXACT_TEST_BIAS) {
                    tie_ct++;
                }
                tailp += lastp2;
                break;
            }
            centerp += lastp2;
            if (centerp == INFINITY) {
                return 0;
            }
        }
        if ((centerp == 0) && (!midp)) {
            return 1;
        }
        while (curr_hets_t2 > 1.5) {
            curr_homr_t2 += 1;
            curr_homc_t2 += 1;
            lastp2 *= (curr_hets_t2 * (curr_hets_t2 - 1)) / (4 * curr_homr_t2 * curr_homc_t2);
            curr_hets_t2 -= 2;
            preaddp = tailp;
            tailp += lastp2;
            if (tailp <= preaddp) {
                break;
            }
        }
        curr_hets_t1 = obs_hets + 2;
        curr_homr_t1 = obs_homr;
        curr_homc_t1 = obs_homc;
        while (curr_homr_t1 > 0.5) {
            // het_probs[curr_hets + 2] = het_probs[curr_hets] * 4 * curr_homr * curr_homc / ((curr_hets + 2) * (curr_hets + 1))
            lastp1 *= (4 * curr_homr_t1 * curr_homc_t1) / (curr_hets_t1 * (curr_hets_t1 - 1));
            preaddp = tailp;
            tailp += lastp1;
            if (tailp <= preaddp) {
                break;
            }
            curr_hets_t1 += 2;
            curr_homr_t1 -= 1;
            curr_homc_t1 -= 1;
        }
    } else {
        // tail 1 = lower
        while (curr_homr_t2 > 0.5) {
            curr_hets_t2 += 2;
            lastp2 *= (4 * curr_homr_t2 * curr_homc_t2) / (curr_hets_t2 * (curr_hets_t2 - 1));
            curr_homr_t2 -= 1;
            curr_homc_t2 -= 1;
            if (lastp2 < EXACT_TEST_BIAS) {
                if (lastp2 > (1 - 2 * SMALL_EPSILON) * EXACT_TEST_BIAS) {
                    tie_ct++;
                }
                tailp += lastp2;
                break;
            }
            centerp += lastp2;
            if (centerp == INFINITY) {
                return 0;
            }
        }
        if ((centerp == 0) && (!midp)) {
            return 1;
        }
        while (curr_homr_t2 > 0.5) {
            curr_hets_t2 += 2;
            lastp2 *= (4 * curr_homr_t2 * curr_homc_t2) / (curr_hets_t2 * (curr_hets_t2 - 1));
            curr_homr_t2 -= 1;
            curr_homc_t2 -= 1;
            preaddp = tailp;
            tailp += lastp2;
            if (tailp <= preaddp) {
                break;
            }
        }
        curr_hets_t1 = obs_hets;
        curr_homr_t1 = obs_homr;
        curr_homc_t1 = obs_homc;
        while (curr_hets_t1 > 1.5) {
            curr_homr_t1 += 1;
            curr_homc_t1 += 1;
            lastp1 *= (curr_hets_t1 * (curr_hets_t1 - 1)) / (4 * curr_homr_t1 * curr_homc_t1);
            preaddp = tailp;
            tailp += lastp1;
            if (tailp <= preaddp) {
                break;
            }
            curr_hets_t1 -= 2;
        }
    }
    if (!midp) {
        return tailp / (tailp + centerp);
    } else {
        return (tailp - ((1 - SMALL_EPSILON) * EXACT_TEST_BIAS * 0.5) * tie_ct) / (tailp + centerp);
    }
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
        else if ( strcmp(a, "-sam") == 0 ) 
            var = &args.asam_tag_string;
        else if ( strcmp(a, "-hwe") == 0 )
            var = &args.hwe_tag_string;
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

    if ( args.asam_tag_string == NULL )
        args.asam_tag_string = "SAMPLE_LIST";

    if ( args.hwe_tag_string == NULL )
        args.hwe_tag_string = "HWE_p";
    
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
    id = bcf_hdr_id2int(args.hdr, BCF_DT_ID, args.asam_tag_string);
    if ( id == -1 ) {
        kstring_t str = { 0, 0, 0 };
        ksprintf(&str, "##INFO=<ID=%s,Number=A,Type=String,Description=\"Allele affected samples list.\">", args.asam_tag_string);
        bcf_hdr_append(args.hdr, str.s);
        bcf_hdr_sync(args.hdr);
        id = bcf_hdr_id2int(args.hdr, BCF_DT_ID, args.af_tag_string);
        assert(bcf_hdr_idinfo_exists(args.hdr, BCF_HL_INFO, id));
    }

    id = bcf_hdr_id2int(args.hdr, BCF_DT_ID, args.hwe_tag_string);
    if ( id == -1 ) {
        kstring_t str = { 0, 0, 0 };
        ksprintf(&str, "##INFO=<ID=%s,Number=1,Type=Float,Description=\"Hardy weinberg equation test p value.\">", args.hwe_tag_string);
        bcf_hdr_append(args.hdr, str.s);
        bcf_hdr_sync(args.hdr);
        id = bcf_hdr_id2int(args.hdr, BCF_DT_ID, args.hwe_tag_string);
        assert(bcf_hdr_idinfo_exists(args.hdr, BCF_HL_INFO, id));
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
    kstring_t *allele_samples = (kstring_t*)malloc(line->n_allele* sizeof(kstring_t));
    for ( i = 0; i < line->n_allele; ++i ) {
        allele_samples[i].l = allele_samples[i].m = 0;
        allele_samples[i].s = NULL;
    }
    uint64_t n_ref_hom = 0;
    uint64_t n_het  = 0;
    uint64_t n_alt_hom = 0;
    int type;
    for ( i = 0; i < args.n_samples; ++i ) {
        int32_t *b = args.tmp_arr + i*ngt;
        int j;
        type = 0;
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
            int a = bcf_gt_allele(b[j]);
            args.counts[a-1] ++;
            kputs(args.hdr->samples[i], &allele_samples[a]);
            kputc('|', &allele_samples[a]);
            all_counts++;
            type++;
        }
        if ( j == 1 )
            error("Assume all position should be diploid. %s,%s,%d,%d",
                  args.hdr->samples[i], bcf_hdr_id2name(args.hdr, line->rid), line->pos+1,bcf_gt_allele(b[j]-1));
        if ( type == 0 )
            n_ref_hom++;
        else if ( type == 1 )
            n_het++;
        else
            n_alt_hom++;
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

    str.l = 0;
    for ( i = 1; i < line->n_allele; ++i ) {
        if ( i > 1)
            kputc(',', &str);
        if ( allele_samples[i].l ) {
            kputs(allele_samples[i].s, &str);
            free(allele_samples[i].s);
        } else {
            kputc('.', &str);
        }
    }
    free(allele_samples);
    bcf_update_info_string(args.hdr, line, args.asam_tag_string, str.s);

    double p_val = SNPHWE2(n_het, n_ref_hom, n_alt_hom, 0);
    bcf_update_info_float(args.hdr, line, args.hwe_tag_string, &p_val, 1);
    
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


