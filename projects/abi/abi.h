// parse ABIF format
// ref :
// 1) Applied Biosystems Genetic Analysis Data File Format
// 2) io_lib

#ifndef ABI_HEADER
#define ABI_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define word_pack_4bits(x) ( ((x[0] & MBIT)<<24) | ((x[1] & MBIT)<<16) | ((x[2] & MBIT)<<8) | (x[3] & MBIT))

#define ABI_MAGICKEY word_pack_4bits("ABIF")

#define TYPE_N 0x0
#define TYPE_A 0x1
#define TYPE_C 0x2
#define TYPE_G 0x4
#define TYPE_T 0x8

#define MBIT 0xF

struct prob_base {
    uint8_t type;
    uint8_t prob;
};

struct abi_read {
    // bases Number    
    int n_bases;
    // peaks number
    int n_peaks; 
    // bases
    uint8_t *bases;
    // base = MAX[trace[base_pos]]
    uint32_t *bases_pos;
    
    // peaks
    uint8_t *traceA;
    uint8_t *traceC;
    uint8_t *traceG;
    uint8_t *traceT;
    // maximun trace value of any trace
    unit8_t max_trace_value;
    // baseline offset, usually be zero
    int32_t baseline;

    // miscellaneous information
    uint8_t *comments;

    // probability
    struct prob_base *probs;
};


struct abi_read *abi_read(const char *fn);

void abi_destroy(struct abi_read *read);

#endif
