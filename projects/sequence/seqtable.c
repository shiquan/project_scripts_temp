#include <stdio.h>
#include <stdlib.h>

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
}
