// vcfspecs - collections of plugins to handle vcf files
// The idea of vcfspecs is construct a stable program body and load dynamic libraries
// to handle different functions. The environment path VCFSPECS is required to run
// vcfspecs. The toolkits specified by users should be defined by dynamics libraries and
// keep in VCFSPECS path. Here is a demo.
//
//  The cleantags.so should be kept in VCFSPECS before run this command,
//
//  vcfspecs cleantags in.vcf.gz -O z -o out.vcf.gz
//
//  The vcfspecs will check cleantags.so first, and load it, and delieve arguemnts to
//  parse_args();

#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>
#include <>
