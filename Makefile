PROG= allele_freqs \
	seqtrim \
	split_barcode \
	umi_parser \
	dyncut_adaptor \
	sam_parse_uid \
	retrievebed \
	vcfeva	\
	CNV_frequency_from_samples \
	CNV_regions_format_per_sample \
	comp_ref_trans \
	rs_finder \
	bamdst_depth_retrieve

all: $(PROG)

HTSDIR = htslib-1.5
include $(HTSDIR)/htslib.mk
include $(HTSDIR)/htslib_static.mk
HTSLIB = $(HTSDIR)/libhts.a
BGZIP  = $(HTSDIR)/bgzip
TABIX  = $(HTSDIR)/tabix
HTSLIB_LDFLAGS = $(HTSLIB_static_LDFLAGS)
HTSLIB_LIBS = $(HTSLIB_static_LIBS)

ifeq "$(shell uname -s)" "Darwin"
DYNAMIC_FLAGS = -Wl,-export_dynamic
else
DYNAMIC_FLAGS = -rdynamic
endif

CC       = gcc
CFLAGS   = -Wall -Wc++-compat -O2
DEBUG_CFLAGS   = -g -Wall -O0
DFLAGS   = -lz -lm -lbz2 -llzma -lcrypto -lcurl -pthread $(DYNAMIC_FLAGS)
INCLUDES = -Iinclude/ -I. -I$(HTSDIR)/

all:$(PROG)

ifneq "$(wildcard .git)" ""
PACKAGE_VERSION := $(shell git describe --tags)
DOC_VERSION :=  $(shell git describe --tags)+
DOC_DATE := $(shell date +'%Y-%m-%d %R %Z')
pkg_version.h: $(if $(wildcard pkg_version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat pkg_version.h)),,force))
endif
pkg_version.h:
	echo '#define PROJECTS_VERSION "$(PACKAGE_VERSION)"' > $@


.SUFFIXES:.c .o
.PHONY:all clean clean-all clean-plugins distclean install lib tags test testclean force plugins docs

force:

.c.o:
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

mk: pkg_version.h $(HTSLIB)
	-mkdir -p bin

CNV_frequency_from_samples: mk
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o bin/$@  projects/CNV_database/cnv_frequency_from_samples.c lib/number.c lib/sort_list.c lib/cnv_bed.c  $(HTSLIB)

CNV_regions_format_per_sample: mk
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o bin/$@  projects/CNV_database/CNV_regions_format_per_sample.c lib/number.c lib/sort_list.c lib/cnv_bed.c  $(HTSLIB)

rs_finder: mk
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o bin/$@  projects/rs/rs_finder.c lib/number.c $(HTSLIB)

allele_freqs: mk
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o bin/allele_freqs_count projects/vcf/allele_freqs_count.c  $(HTSLIB)

seqtrim: mk
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o bin/seqtrim projects/sequence/seqtrim/seqtrim.c lib/sequence.c $(HTSLIB)

split_barcode: mk
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o bin/split_barcode projects/sequence/split_barcode/split_barcode.c lib/number.c lib/fastq.c $(HTSLIB)

umi_parser: mk
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o bin/umi_parser projects/sequence/umi_parser/umi_parser.c lib/number.c $(HTSLIB)	

dyncut_adaptor: mk
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o bin/dyncut_adaptor projects/sequence/dyncut_adaptor/dyncut_adaptor_trim_uid.c lib/number.c lib/fastq.c $(HTSLIB)

sam_parse_uid: mk
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o bin/$@ projects/bam/parse_UID_tag.c lib/number.c lib/sequence.c $(HTSLIB)

retrievebed: mk
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -DGENEPRED_TEST_MAIN -o bin/$@ lib/genepred.c lib/number.c lib/sort_list.c  $(HTSLIB)

vcfeva: mk
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o bin/$@ projects/vcf/vcfeva.c $(HTSLIB)

comp_ref_trans: mk
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o bin/$@ projects/gene_regions/check_genepred_transcripts.c lib/ksw.c lib/genepred.c lib/sequence.c lib/number.c lib/kthread.c lib/faidx_def.c $(HTSLIB)

bamdst_depth_retrieve: mk
	$(CC) $(DEBUG_CFLAGS) $(DFLAGS) $(INCLUDES) -o bin/$@ projects/depths/bamdst_depth_retrieve.c lib/number.c $(HTSLIB)

clean: testclean
	-rm -f gmon.out *.o *~ $(PROG) pkg_version.h  version.h
	-rm -rf bin/*.dSYM test/*.dSYM
	-rm -rf bin/

testclean:
	-rm -f test/*.o test/*~ $(TEST_PROG)

distclean: clean
	-rm -f TAGS

clean-all: clean clean-htslib

tags:
	ctags -f TAGS lib/*.[ch]
