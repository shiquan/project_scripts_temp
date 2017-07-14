PROG=  mk allele_freqs seqtrim split_barcode umi_parser

all: $(PROG)

HTSDIR = htslib-1.4.1
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
PACKAGE_VERSION := $(shell git describe --always --dirty)
DOC_VERSION :=  $(shell git describe --always)+
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

allele_freqs: allele_freqs.o
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o bin/allele_freqs_count projects/vcf/allele_freqs_count.c  $(HTSLIB)

seqtrim:
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o bin/seqtrim projects/sequence/seqtrim.c lib/sequence.c $(HTSLIB)

split_barcode:
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o bin/split_barcode projects/sequence/split_barcode.c lib/number.c $(HTSLIB)

umi_parser:
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o bin/umi_parser projects/sequence/umi_parser.c lib/number.c $(HTSLIB)	

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
