// convert ab1 file to base and trace value table
// reference : http://bioinformatics.oxfordjournals.org/content/14/1/92.full.pdf

#include <stdio.h>
#include <stdlib.h>
// require io_lib
#include <io_lib/Read.h>
#include <errno.h>
#include "utils.h"

int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) == 1)
	return 1;

    

    int i;
    Read *read = read_reading( argv[1], TT_ABI);
    if ( read == NULL )
	error("Failed to open %s : %s", argv[1], strerror(errno));

    char iupac(char A1, char A2)
{
    char iu = 'N';
    if (A1 == 'A' || A1 == 'a') {
	switch (A2) {
	    case 'C':
	    case 'c':
		iu = 'M';
		break;
	    case 'G':
	    case 'g':
		iu = 'R';
		break;
	    case 'T':
	    case 't':
		iu = 'W';
		break;
	    default:
		iu = 'A';
		break;
	}
    } else if (A1 == 'C' || A1 == 'c') {
	switch (A2) {
	    case 'A':
	    case 'a':
		iu = 'M';
		break;
	    case 'G':
	    case 'g':
		iu = 'S';
		break;
	    case 'T':
	    case 't':
		iu = 'Y';
		break;
	    default:
		iu = 'C';
		break;
	}
    } else if (A1 == 'G' || A1 == 'g') {
	switch (A2) {
	    case 'A':
	    case 'a':
		iu = 'R';
		break;
	    case 'C':
	    case 'c':
		iu = 'S';
		break;
	    case 'T':
	    case 't':
		iu = 'K';
		break;
	    default:
		iu = 'G';
		break;
	}
    } else if (A1 == 'T' || A1 == 't') {
	switch (A2) {
	    case 'A':
	    case 'a':
		iu = 'W';
		break;
	    case 'C':
	    case 'c':
		iu = 'Y';
		break;
	    case 'G':
	    case 'g':
		iu = 'K';
		break;
	    default:
		iu = 'T';
		break;
	}
    } else {
	switch (A2) {
	    case 'A':
	    case 'a':
		iu = 'A';
		break;
	    case 'C':
	    case 'c':
		iu = 'C';
		break;
	    case 'G':
	    case 'g':
		iu = 'G';
		break;
	    case 'T':
	    case 't':
		iu = 'T';
		break;
	    default:
		iu = 'N';
		break;
	}
    }
    return iu;
}
int main(int argc, char **argv)
{
    Read* read;
    int i;

    if (argc != 2) {
    fprintf(stderr, "Usage: trace_dump <trace file>\n");
    return 1;
    }


    read = read_reading( argv[1], TT_ANY );


    if (read == NULL) {
    fprintf(stderr, "Tracedump was unable to open file %s\n", argv[1] );
    return 1;
    }

    struct peak {
	int value;
	int type;
    };
    struct peak last = { 0, 'N'};
    for (i = 0; i < read->NBases; i++) {
	int peak_A = read->traceA[read->basePos[i]] - read->baseline;
	int peak_C = read->traceC[read->basePos[i]] - read->baseline;
	int peak_G = read->traceG[read->basePos[i]] - read->baseline;
	int peak_T = read->traceT[read->basePos[i]] - read->baseline;

	char base;
	struct peak max_peak = { 0, 'N' };
	struct peak second_peak = {0, 'N'};
	if ( max_peak.value < peak_A ) {
	    second_peak.value = max_peak.value;
	    second_peak.type = max_peak.type;
	    max_peak.value = peak_A;
	    max_peak.type = 'A';
	} else if ( second_peak.value < peak_A ) {
	    second_peak.value = peak_A;
	    second_peak.type = 'A';
	}

	if ( max_peak.value < peak_C ) {
	    second_peak.value = max_peak.value;
	    second_peak.type = max_peak.type;
	    max_peak.value = peak_C;
	    max_peak.type = 'C';
	} else if ( second_peak.value < peak_C ) {
	    second_peak.value = peak_C;
	    second_peak.type = 'C';
	}
	
	if ( max_peak.value < peak_G ) {
	    second_peak.value = max_peak.value;
	    second_peak.type = max_peak.type;
	    max_peak.value = peak_G;
	    max_peak.type = 'G';
	} else if ( second_peak.value < peak_G ) {
	    second_peak.value = peak_G;
	    second_peak.type = 'G';
	}

	if ( max_peak.value < peak_T ) {
	    second_peak.value = max_peak.value;
	    second_peak.type = max_peak.type;
	    max_peak.value = peak_T;
	    max_peak.type = 'T';
	} else if ( second_peak.value < peak_T ) {
	    second_peak.value = peak_T;
	    second_peak.type = 'T';
	}

	if ( last.type == second_peak.type && second_peak.value < last.value) {
	    second_peak.type = 'N';
	    second_peak.value = 0;
	}
	if ( max_peak.value == 0)
	    base = 'N';
	else if ( second_peak.value < 100)
	    base = max_peak.type;
	else
	    base = iupac(max_peak.type,second_peak.type);	
	printf("%4d %c %c %10d %10d %10d %10d %10d\n", i, base, last.type, last.value, peak_A, peak_C, peak_G, peak_T);
	last.value = max_peak.value;
	last.type = max_peak.type;
    }

    read_deallocate(read);

    return 0;

}
