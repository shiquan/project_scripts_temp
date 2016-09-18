#include <stdio.h>
#include <stdlib.h>

unsigned is_litter_endian(void)
{
    const union {
	uint32_t i;
	uint8_t c[4];
    } one = { 1 };
    return one.c[0];
}

unsigned is_64bits(void)
{
    return sizeof(void*) == 8;
}
