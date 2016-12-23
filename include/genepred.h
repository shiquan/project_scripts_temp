#ifndef GENEPRED_HEADER
#define GENEPRED_HEADER
#include "htslib/tbx.h"
#include "sort_list.h"

struct genepred_line {
    struct genepred_line *next;
    // Chromosome id.
    int tid;
    int txstart;
    int txend;
    char strand;
    // Usually the transcript name.
    char *name1;
    // Gene name.
    char *name2;

    int cdsstart;
    int cdsend;

    // Forward length and backward length are the UTRs in the transcript, for different
    // strand the forward and backward could be UTR5' or UTR3'; and most importantly,
    // for noncoding transcript, no UTRs, so the forward and backward length should
    // always be 0.
    int forward_length;
    int backward_length;

    // Length of this transcript.
    int reference_length;

    int exon_count;
    // [start, end]
    int *exons[2];
    // Offsets location on gene.
    int *dna_ref_offsets[2];
    // Transcript locations of each block, coding transcripts consist of UTRs and CDS,
    // so loc[0,1] is not the edge of coding sequences.
    int *loc[2];    
};

struct genepred_format {
    int chrom;
    // Usually name1 is transcript, name2 is gene.
    int name1;
    int name2;
    int strand;
    int txstart;
    int txend;
    int cdsstart;
    int cdsend;
    int exoncount;
    int exonstarts;
    int exonends;
};

struct genepred_spec {
    const char *data_fname;
    tbx_t *idx;
    void *chr_hash;
    char **names;
    void *tran_hash;
    void *gene_hash;    
};

struct genepred_spec *genepred_spec_init();
void genepred_spec_destroy(struct genepred_spec *spec);

extern int genepred_load_data(struct genepred_spec *spec, const char *fname);
extern struct genepred_line *genepred_retrieve_gene(struct genepred_spec *spec, const char *name);
extern struct genepred_line *genepred_retrieve_trans(struct genepred_spec *spec, const char *name);


#endif
