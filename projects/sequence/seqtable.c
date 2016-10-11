#include <stdio.h>
#include <stdlib.h>

char *codon_array[][2] = {
    { "UUU", "Phe" },
    { "UUC", "Phe" },
    { "UUA", "Leu" },
    { "UUG", "Leu" },
    { "UCU", "Ser" },
    { "UCC", "Ser" },
    { "UCA", "Ser" },
    { "UCG", "Ser" },
    { "UAU", "Tyr" },
    { "UAC", "Tyr" },
    { "UAA", "Stop" },
    { "UAG", "Stop" },
    { "UGU", "Cys" },
    { "UGC", "Cys" },
    { "UGA", "Stop" },
    { "UGG", "Trp" },
    { "CUU", "Leu" },
    { "CUC", "Leu" },
    { "CUA", "Leu" },
    { "CUG", "Leu" },
    { "CCU", "Pro" },
    { "CCC", "Pro" },
    { "CCA", "Pro" },
    { "CCG", "Pro" },
    { "CAU", "His" },
    { "CAC", "His" },
    { "CAA", "Gln" },
    { "CAG", "Gln" },
    { "CGU", "Arg" },
    { "CGC", "Arg" },
    { "CGA", "Arg" },
    { "CGG", "Arg" },
    { "AUU", "Ile" },
    { "AUC", "Ile" },
    { "AUA", "Ile" },
    { "AUG", "Met" },
    { "ACU", "Thr" },
    { "ACC", "Thr" },
    { "ACA", "Thr" },
    { "ACG", "Thr" },
    { "AAU", "Asn" },
    { "AAC", "Asn" },
    { "AAA", "Lys" },
    { "AAG", "Lys" },
    { "AGU", "Ser" },
    { "AGC", "Ser" },
    { "AGA", "Arg" },
    { "AGG", "Arg" },
    { "GUU", "Val" },
    { "GUC", "Val" },
    { "GUA", "Val" },
    { "GUG", "Val" },
    { "GCU", "Ala" },
    { "GCC", "Ala" },
    { "GCA", "Ala" },
    { "GCG", "Ala" },
    { "GAU", "Asp" },
    { "GAC", "Asp" },
    { "GAA", "Glu" },
    { "GAG", "Glu" },
    { "GGU", "Gly" },
    { "GGC", "Gly" },
    { "GGA", "Gly" },
    { "GGG", "Gly" },
};

#define C4_Stop 0
#define C4_Phe  1
#define C4_Leu  2
#define C4_Ser  3
#define C4_Tyr  4
#define C4_Cys  5
#define C4_Trp  6
#define C4_Pro  7
#define C4_His  8
#define C4_Gln  9
#define C4_Arg 10
#define C4_Ile 11
#define C4_Met 12
#define C4_Thr 13
#define C4_Asn 14
#define C4_Lys 15
#define C4_Val 16
#define C4_Ala 17
#define C4_Asp 18
#define C4_Glu 19
#define C4_Gly 20

const char *codon_names[] = {
    "Stop", "Phe", "Leu", "Ser", "Tyr", "Cys", "Trp", "Pro", "His", "Gln",
    "Arg", "Ile", "Met", "Thr", "Asn", "Lys", "Val", "Ala", "Asp", "Glu", "Gly",
};

const int codon_matrix[4][4][4] =  {
    {
        { C4_Lys, C4_Asn, C4_Lys, C4_Asn, },
        { C4_Thr, C4_Thr, C4_Thr, C4_Thr, },
        { C4_Arg, C4_Ser, C4_Arg, C4_Ser, },
        { C4_Ile, C4_Ile, C4_Met, C4_Ile, },
    },
    {
        { C4_Gln, C4_His, C4_Gln, C4_His, },
        { C4_Pro, C4_Pro, C4_Pro, C4_Pro, },
        { C4_Arg, C4_Arg, C4_Arg, C4_Arg, },
        { C4_Leu, C4_Leu, C4_Leu, C4_Leu, },
    },
    {
        { C4_Glu, C4_Asp, C4_Glu, C4_Asp, },
        { C4_Ala, C4_Ala, C4_Ala, C4_Ala, },
        { C4_Gly, C4_Gly, C4_Gly, C4_Gly, },
        { C4_Val, C4_Val, C4_Val, C4_Val, },
    },
    {
        { C4_Stop, C4_Tyr, C4_Stop, C4_Tyr, },
        { C4_Ser, C4_Ser, C4_Ser, C4_Ser, },
        { C4_Stop, C4_Cys, C4_Trp, C4_Cys, },
        { C4_Leu, C4_Phe, C4_Leu, C4_Phe, },
    },
};

static inline uint8_t seq2code4(uint8_t seq)
{
    static const uint8_t seq2num_table[256] = {
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    };
    return seq2num_table[seq];
}

int main()
{
    int i;
    for ( i = 0; i < 256; ++i ) {
        switch(i) {
            case 'a':
            case 'A':
                fputs("0, ", stdout);
                break;
                
            case 'c':
            case 'C':
                fputs("1, ", stdout);
                break;

            case 'g':
            case 'G':
                fputs("2, ", stdout);
                break;
            case 'u':
            case 'U':
            case 't':
            case 'T':
                fputs("3, ", stdout);
                break;

            default:
                fputs("4, ", stdout);
                break;
        }
        if ( i >= 15 && (i+1) %16 == 0 )
            fputc('\n', stdout);            
    }
    
    fputc('\n', stdout);
    fputc('\n', stdout);


    char *codon[4][4][4];

    int length = sizeof(codon_array)/ sizeof(codon_array[0]);

    int j;
    for (i = 0; i < length; ++i) {
        char *name1 = codon_array[i][0];
        char *name2 = codon_array[i][1];
        // codon[seq2code4(name1[0])][seq2code4(name1[1])][seq2code4(name1[2])] = name2;
        printf("%s\t",codon_names[codon_matrix[seq2code4(name1[0])][seq2code4(name1[1])][seq2code4(name1[2])]]);
        fputs(name2, stdout);
        fputs("\n", stdout);
    }

    /* int k; */
    /* for ( i = 0; i < 4; i++) { */
    /*     fputs("{\n", stdout); */
    /*     for (j = 0; j < 4; j ++ ) { */
    /*         fputs("{ ", stdout); */
    /*         for (k = 0; k < 4; k++ ) { */
    /*             fputs("C4_",stdout); */
    /*             fputs(codon[i][j][k], stdout); */
    /*             fputs(", ",stdout); */
    /*         } */
    /*         fputs("},\n", stdout); */
    /*     } */
    /*     fputs("},", stdout); */
    /*     fputc('\n', stdout); */
    /* } */
    return 0;
}
