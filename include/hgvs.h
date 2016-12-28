#ifndef HGVSLIB_HEADER
#define HGVSLIB_HEADER

enum func_region_type {
    _func_region_promote_to_int = -1,
    func_region_unknown,
    //func_region_split_sites,
    func_region_cds,
    func_region_noncoding,
    func_region_intron,
    func_region_utr5,
    func_region_utr3,
    func_region_intergenic,
};

struct var_func_type {
    enum func_region_type func;
    enum var_type vartype;
    int start_flag; // region start flag for cds or intron count,
    // flag & 1 == 1 for intron else for cds/utr/exon
};

//  hgvs_core keeps transcript/protein name, the format of hgvs name is construct by
//                                    l_name   l_name2  l_type
//                                    |        |        |
//  [transcript|protein|gene|ensemble](gene_id):"prefix".postion...
//
//  l_name          l_type
//  |               |
//  [chrom]:"prefix".postion...
struct hgvs_core {    
    uint16_t l_name; // name offset in the data cache, for chrom l_name == 0;
    uint16_t l_name2; // name2 offset in data cache, if no name2, l_name2 should be 0
    uint16_t l_type; // the byte after type offset, usually offset of '.'
    kstring_t str;
    struct var_func_type type;   
};

// Dynamic allocate and init technology.
// For struct { int l, m, i; void *a }, a is the cached array of predefined type. l for used length, m for max length,
// i for inited length. And i always >= l, and m always >= i. If l == m, the cache array should be reallocated by a new
// memory size. This structure used to reuse the cache complex structures other than points.
struct hgvs {
    int l, m, i; 
    struct hgvs_core *a;
};
struct hgvs_cache {
    int l, m, i; // l == n_allele -1
    struct hgvs *a;
};
#define HGVS_CACHE_INIT { 0, 0, 0, 0 }
// HGVS nomenclature : 
// DNA recommandations 
// - substitution variant, 
// Format: “prefix”“position_substituted”“reference_nucleoride””>”new_nucleoride”, e.g. g.123A>G
// - deletion variant, 
// Format: “prefix”“position(s)_deleted”“del”, e.g. g.123_127del
// - duplication variant,
// Format: “prefix”“position(s)_duplicated”“dup”, e.g. g.123_345dup
// - insertion variant,
// Format: “prefix”“positions_flanking”“ins”“inserted_sequence”, e.g. g.123_124insAGC
// - inversion variant,
// Format: “prefix”“positions_inverted”“inv”, e.g. g.123_345inv
// - conversion variant,
// Format: “prefix”“positions_converted”“con”“positions_replacing_sequence”, e.g. g.123_345con888_1110
// - deletion-insertion variant,
// Format: “prefix”“position(s)_deleted”“delins”“inserted_sequence”, e.g. g.123_127delinsAG
// - repeated sequences variant,
// Format (repeat position): “prefix”“position_repeat_unit””["”copy_number””]”, e.g. g.123_125[36]
// 
// [ reference : http:// www.HGVS.org/varnomen ]
// hgvs_variant_type is different from var_type, for var_type is more concern about function, but
// hgvs_variant_type is just describe the change on DNA level
enum hgvs_variant_type {
    var_type_nonref = -1, // for gatk <NONREF> allele
    var_type_ref = 0,
    var_type_snp,
    var_type_dels,
    var_type_ins,
    var_type_delins,
    var_type_copy, // dup for ins
    var_type_complex, 
    var_type_unknow,
};
// descriptions struct of hgvs name
struct hgvs_des {
    // if type ==  var_type_ref, hgvs name generater will skip to construct a name, only type will be inited, 
    // DONOT use any other value in this struct then 
    enum hgvs_variant_type type;
    uint8_t strand; // for plus strand start <= end, for minus start >= end
    int32_t start; // position on ref
    // for var_type_snp end == start, for var_type_dels end > start, for var_type_ins end = start +1, 
    // for delvar_type_ins end > start 
    int32_t end;
    // for var_type_snp ref_length == 1, for var_type_dels ref_length == length(dels),
    // for var_type_ins ref_length == 0, for var_type_delins ref_length == length(var_type_dels) 
    int ref_length;
    char *ref;
    // for var_type_snp alt_length == 1, for var_type_dels alt_length == 0, for var_type_ins alt_length ==
    // length(var_type_ins), for delvar_type_ins alt_length == length(var_type_ins) 
    int alt_length;
    char *alt;
    // the copy number only valid if variants type is var_type_copy, for most case, my algrithm only check
    // mark var_type_dels or var_type_ins in the first step, and check the near sequences from refseq 
    // databases and if there is copy number changed (more or less), the variants will update to var_type_copy 
    int tandam_repeat_number;
    int tandam_repeat_number_ori;
};

#endif
