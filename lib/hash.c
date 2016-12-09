#include <stdio.h>
#include <time.h>

/* Can be obtained from http://www.pobox.com/~qed/pstdint.h , but can be
   replaced with <stdint.h> if appropriate */
#include <stdint.h>

#ifndef CLOCKS_PER_SEC
# ifdef CLK_TCK
#  define CLOCKS_PER_SEC (CLK_TCK)
# endif
#endif

/* Hash performance test for various hash functions */

/* http://burtleburtle.net/bob/hash/doobs.html */

typedef uint32_t ub4;   /* unsigned 4-byte quantities */
typedef uint8_t  ub1;   /* unsigned 1-byte quantities */

#undef get32bits
#if (defined(__GNUC__) && defined(__i386__)) || defined(__WATCOMC__) \
  || defined(_MSC_VER) || defined (__BORLANDC__) || defined (__TURBOC__)
#define get32bits(d) (*((const uint32_t *) (d)))
#endif

#if !defined (get32bits)
#define get32bits(d) (  (uint32_t)(d)[0]                \
                     + ((uint32_t)(d)[1]<<UINT32_C(8))  \
                     + ((uint32_t)(d)[2]<<UINT32_C(16)) \
                     + ((uint32_t)(d)[3]<<UINT32_C(24)) )
#endif

#define hashsize(n) ((ub4)1<<(n))
#define hashmask(n) (hashsize(n)-1)

/*
--------------------------------------------------------------------
mix -- mix 3 32-bit values reversibly.
For every delta with one or two bits set, and the deltas of all three
  high bits or all three low bits, whether the original value of a,b,c
  is almost all zero or is uniformly distributed,
* If mix() is run forward or backward, at least 32 bits in a,b,c
  have at least 1/4 probability of changing.
* If mix() is run forward, every bit of c will change between 1/3 and
  2/3 of the time.  (Well, 22/100 and 78/100 for some 2-bit deltas.)
mix() was built out of 36 single-cycle latency instructions in a 
  structure that could supported 2x parallelism, like so:
      a -= b; 
      a -= c; x = (c>>13);
      b -= c; a ^= x;
      b -= a; x = (a<<8);
      c -= a; b ^= x;
      c -= b; x = (b>>13);
      ...
  Unfortunately, superscalar Pentiums and Sparcs can't take advantage 
  of that parallelism.  They've also turned some of those single-cycle
  latency instructions into multi-cycle latency instructions.  Still,
  this is the fastest good hash I could find.  There were about 2^^68
  to choose from.  I only looked at a billion or so.
--------------------------------------------------------------------
*/
#define mix(a,b,c) \
{ \
  a -= b; a -= c; a ^= (c>>13); \
  b -= c; b -= a; b ^= (a<<8); \
  c -= a; c -= b; c ^= (b>>13); \
  a -= b; a -= c; a ^= (c>>12);  \
  b -= c; b -= a; b ^= (a<<16); \
  c -= a; c -= b; c ^= (b>>5); \
  a -= b; a -= c; a ^= (c>>3);  \
  b -= c; b -= a; b ^= (a<<10); \
  c -= a; c -= b; c ^= (b>>15); \
}

/*
--------------------------------------------------------------------
hash() -- hash a variable-length key into a 32-bit value
  k       : the key (the unaligned variable-length array of bytes)
  len     : the length of the key, counting by bytes
  initval : can be any 4-byte value
Returns a 32-bit value.  Every bit of the key affects every bit of
the return value.  Every 1-bit and 2-bit delta achieves avalanche.
About 6*len+35 instructions.

The best hash table sizes are powers of 2.  There is no need to do
mod a prime (mod is sooo slow!).  If you need less than 32 bits,
use a bitmask.  For example, if you need only 10 bits, do
  h = (h & hashmask(10));
In which case, the hash table should have hashsize(10) elements.

If you are hashing n strings (ub1 **)k, do it like this:
  for (i=0, h=0; i<n; ++i) h = hash( k[i], len[i], h);

By Bob Jenkins, 1996.  bob_jenkins@burtleburtle.net.  You may use this
code any way you wish, private, educational, or commercial.  It's free.

See http://burtleburtle.net/bob/hash/evahash.html
Use for hash table lookup, or anything where one collision in 2^^32 is
acceptable.  Do NOT use for cryptographic purposes.
--------------------------------------------------------------------
*/

static uint32_t hashBobJenkinsUpdate(const char * k, int length, unsigned int initval) {
   register ub4 a,b,c,len;

   /* Set up the internal state */
   len = length;
   a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
   c = initval;         /* the previous hash value */

   /*---------------------------------------- handle most of the key */
   while (len >= 12) {
      a += get32bits (k);
      b += get32bits (k+4);
      c += get32bits (k+8);
      mix(a,b,c);
      k += 12; len -= 12;
   }

   /*------------------------------------- handle the last 11 bytes */
   c += length;
   switch(len) {            /* all the case statements fall through */
      case 11: c+=((ub4)k[10]<<24);
      case 10: c+=((ub4)k[9]<<16);
      case 9 : c+=((ub4)k[8]<<8);
        /* the first byte of c is reserved for the length */
      case 8 : b+=((ub4)k[7]<<24);
      case 7 : b+=((ub4)k[6]<<16);
      case 6 : b+=((ub4)k[5]<<8);
      case 5 : b+=k[4];
      case 4 : a+=((ub4)k[3]<<24);
      case 3 : a+=((ub4)k[2]<<16);
      case 2 : a+=((ub4)k[1]<<8);
      case 1 : a+=k[0];
      /* case 0: nothing left to add */
   }
   mix(a,b,c);
   /*-------------------------------------------- report the result */
   return c;
}

uint32_t hashBobJenkins (const char * k, int length) {
	return hashBobJenkinsUpdate (k, length, 0);
}

/* ======================================================================== */

static uint32_t crc_16_table[16] = {
  0x0000, 0xCC01, 0xD801, 0x1400, 0xF001, 0x3C00, 0x2800, 0xE401,
  0xA001, 0x6C00, 0x7800, 0xB401, 0x5000, 0x9C01, 0x8801, 0x4400 
};

/* 
 *	This code was found at: http://wannabe.guru.org/alg/node191.html
 *  and still exists here: http://www.fearme.com/misc/alg/node191.html
 * 
 *  this source code is based on Rex and Binstock which, in turn,
 *  acknowledges William James Hunt.
 *
 *	According to the site this CRC uses the polynomial x^16+x^5+x^2+1.
 *	Unfortunately, DOCSIS uses x^16+x^12+x^5+1.  D'oh!
 */

static uint32_t GetCRC16Update (uint32_t start_crc, const char * data_stream, int length) {
uint32_t crc = start_crc;
uint32_t r;

  /* while there is more data to process */
  while (length-- > 0) {

    /* compute checksum of lower four bits of *data_stream */
    r = crc_16_table[crc & 0xF];
    crc = (crc >> 4) & 0x0FFF;
    crc ^= r ^ crc_16_table[*data_stream & 0xF];

    /* now compute checksum of upper four bits of *data_stream */
    r = crc_16_table[crc & 0xF];
    crc = (crc >> 4) & 0x0FFF;
    crc ^= r ^ crc_16_table[(*data_stream >> 4) & 0xF];

    /* next... */
    data_stream++;
  }
 
  return crc;
}

uint32_t GetCRC16 (const char * data_stream, int length) {
	return GetCRC16Update (0, data_stream, length);
}

/* ======================================================================== */

static uint32_t crc_table[256];

/*  This code was found at: 
 *  http://cell.onecall.net/cell-relay/publications/software/CRC/32bitCRC.c.html   
 */

/*                                                                      */
/* crc32h.c -- package to compute 32-bit CRC one byte at a time using   */
/*             the high-bit first (Big-Endian) bit ordering convention  */
/*                                                                      */
/* Synopsis:                                                            */
/*  gen_crc_table() -- generates a 256-word table containing all CRC    */
/*                     remainders for every possible 8-bit byte.  It    */
/*                     must be executed (once) before any CRC updates.  */
/*                                                                      */
/*  unsigned update_crc(crc_accum, data_blk_ptr, data_blk_size)         */
/*           unsigned crc_accum; char *data_blk_ptr; int data_blk_size; */
/*           Returns the updated value of the CRC accumulator after     */
/*           processing each byte in the addressed block of data.       */
/*                                                                      */
/*  It is assumed that an unsigned long is at least 32 bits wide and    */
/*  that the predefined type char occupies one 8-bit byte of storage.   */
/*                                                                      */
/*  The generator polynomial used for this version of the package is    */
/*  x^32+x^26+x^23+x^22+x^16+x^12+x^11+x^10+x^8+x^7+x^5+x^4+x^2+x^1+x^0 */
/*  as specified in the Autodin/Ethernet/ADCCP protocol standards.      */
/*  Other degree 32 polynomials may be substituted by re-defining the   */
/*  symbol POLYNOMIAL below.  Lower degree polynomials must first be    */
/*  multiplied by an appropriate power of x.  The representation used   */
/*  is that the coefficient of x^0 is stored in the LSB of the 32-bit   */
/*  word and the coefficient of x^31 is stored in the most significant  */
/*  bit.  The CRC is to be appended to the data most significant byte   */
/*  first.  For those protocols in which bytes are transmitted MSB      */
/*  first and in the same order as they are encountered in the block    */
/*  this convention results in the CRC remainder being transmitted with */
/*  the coefficient of x^31 first and with that of x^0 last (just as    */
/*  would be done by a hardware shift register mechanization).          */
/*                                                                      */
/*  The table lookup technique was adapted from the algorithm described */
/*  by Avram Perez, Byte-wise CRC Calculations, IEEE Micro 3, 40 (1983).*/

/* generate the table of CRC remainders for all possible bytes */

#define CRC32POLYNOMIAL 0x04c11db7L

static void GenerateCRC32Table (void) { 
register int i, j;  
register uint32_t crc_accum;

	for ( i = 0;  i < 256;  i++ ) { 
		crc_accum = ( (unsigned long) i << 24 );
        for ( j = 0;  j < 8;  j++ ) { 
			if ( crc_accum & 0x80000000L ) {
				crc_accum = ( crc_accum << 1 ) ^ CRC32POLYNOMIAL;
			} else {
				crc_accum = ( crc_accum << 1 );
			}
		}
		crc_table[i] = crc_accum; 
	}
	return;
}

/* update the CRC on the data block one byte at a time */

static uint32_t UpdateCRC32 (uint32_t crc_accum, const char *data_blk_ptr, int data_blk_size) {
register int j;
register uint8_t i;

	for (j = 0; j < data_blk_size; j++) {
		i = (crc_accum >> 24) ^ *data_blk_ptr++;
		crc_accum = (crc_accum << 8) ^ crc_table[i]; 
	}
	return crc_accum;
}

uint32_t GetCRC32 (const char * data_stream, int length) {
	return UpdateCRC32 (0, data_stream, length);
}

/* ======================================================================== */

/* Performs two parallel CRC-32 on even and odd bytes of the input, then
   combines the two in a further CRC-32 calculation */

uint32_t GetCRC32PH (const char *data_blk_ptr, int data_blk_size) { 
int j;
uint8_t i0, i1;
uint32_t crc_accum0 = 0, crc_accum1 = UINT32_C (0x23456789);

	if (data_blk_size & 1) crc_accum0 ^= *data_blk_ptr++;
	for (j = 1; j < data_blk_size; j+=2) {
		i0 = ((crc_accum0 >> 24) ^ *data_blk_ptr++);
		i1 = ((crc_accum1 >> 24) ^ *data_blk_ptr++);
		crc_accum0 = (crc_accum0 << 8) ^ crc_table[i0]; 
		crc_accum1 = (crc_accum1 << 8) ^ crc_table[i1]; 
	}
	return crc_accum0 + crc_accum1;
}

/* ======================================================================== */

/* Fowler / Noll / Vo (FNV) Hash 

   http://www.isthe.com/chongo/tech/comp/fnv/ */

uint32_t FNVHash (const char * data, int len) {
int i;
uint32_t hash;

	hash = UINT32_C(2166136261);
	for (i=0; i < len; i++) {
		hash = (UINT32_C(16777619) * hash) ^ data[i];
	}
	return hash;
}

/* ======================================================================== */

/*
 * http://burtleburtle.net/bob/hash/doobs.html
 */

uint32_t oneAtATimeHash (const char * s, int len) {
int32_t hash;
int i;

	for (hash = 0, i = 0; i < len; i++) {
		hash += s[i];
		hash += (hash << 10);
		hash ^= (hash >>  6);	/* Non-portable due to ANSI C */
	}
	hash += (hash <<  3);
	hash ^= (hash >> 11);		/* Non-portable due to ANSI C */
	hash += (hash << 15);
	return (uint32_t) hash;
} 

/* ======================================================================== */

uint32_t oneAtATimeHashPH (const char * s, int len) {
int32_t hash0 = 0, hash1 = INT32_C(0x23456789);
int i;

	if (len & 1) hash1 ^= *s++;

	for (i = 1; i < len; i+=2) {
		hash0 += *s++;
		hash1 += *s++;
		hash0 += (hash0 << 10);
		hash1 += (hash1 << 10);
		hash0 ^= (hash0 >>  6);	/* Non-portable due to ANSI C */
		hash1 ^= (hash1 >>  6);	/* Non-portable due to ANSI C */
	}

	hash0 += hash1;

	hash0 += (hash0 <<  3);
	hash0 ^= (hash0 >> 11);		/* Non-portable due to ANSI C */
	hash0 += (hash0 << 15);
	return (uint32_t) hash0;
} 

/* ======================================================================== */

/* By Paul Hsieh (C) 2004, 2005.  Covered under the Paul Hsieh derivative 
   license. See: 
   http://www.azillionmonkeys.com/qed/weblicense.html for license details.

   http://www.azillionmonkeys.com/qed/hash.html */

#undef get16bits
#if (defined(__GNUC__) && defined(__i386__)) || defined(__WATCOMC__) \
  || defined(_MSC_VER) || defined (__BORLANDC__) || defined (__TURBOC__)
#define get16bits(d) (*((const uint16_t *) (d)))
#endif

#if !defined (get16bits)
#define get16bits(d) ((((uint32_t)(((const uint8_t *)(d))[1])) << 8)\
                       +(uint32_t)(((const uint8_t *)(d))[0]) )
#endif

uint32_t SuperFastHash (const char * data, int len) {
uint32_t hash = 0, tmp;
int rem;

	if (len <= 0 || data == NULL) return 0;

	rem = len & 3;
	len >>= 2;

	/* Main loop */
	for (;len > 0; len--) {
		hash  += get16bits (data);
		tmp    = (get16bits (data+2) << 11) ^ hash;
		hash   = (hash << 16) ^ tmp;
		data  += 2*sizeof (uint16_t);
		hash  += hash >> 11;
	}

	/* Handle end cases */
	switch (rem) {
		case 3:	hash += get16bits (data);
				hash ^= hash << 16;
				hash ^= data[sizeof (uint16_t)] << 18;
				hash += hash >> 11;
				break;
		case 2:	hash += get16bits (data);
				hash ^= hash << 11;
				hash += hash >> 17;
				break;
		case 1: hash += *data;
				hash ^= hash << 10;
				hash += hash >> 1;
	}

	/* Force "avalanching" of final 127 bits */
	hash ^= hash << 3;
	hash += hash >> 5;
	hash ^= hash << 4;
	hash += hash >> 17;
	hash ^= hash << 25;
	hash += hash >> 6;

	return hash;
}

uint32_t SuperFastHashAsm (const char * data, int len) {
uint32_t hash = 0;

	if (len <= 0 || data == NULL) return 0;

#if defined(__WATCOMC__) || defined(_MSC_VER)
	__asm {
		xor   ebx, ebx
		mov   esi, data
		mov   edi, len
		mov   eax, edi
		mov   ecx, edi
		and   edi, 3
		shr   ecx, 2
		jz    L2

		L1:
		movzx ebx, word ptr [esi]   ; 0
		add   eax, ebx              ; 1
		movzx ebx, word ptr [esi+2] ; 0

		shl   ebx, 11               ; 1
		xor   ebx, eax              ; 2
		shl   eax, 16               ; 2

		add   esi, 4                ; 0
		xor   eax, ebx              ; 3
		mov   edx, eax              ; 4

		shr   eax, 11               ; 4
		add   eax, edx              ; 5
		dec   ecx                   ; 2
		jnz   L1

		L2:

		mov   ecx, edi
		cmp   ecx, 1
		jz    lcase1
		cmp   ecx, 2
		jz    lcase2
		cmp   ecx, 3
		jnz   L3

		mov    bx, word ptr [esi]
		add   eax, ebx
		xor   ebx, ebx
		mov   edx, eax
		shl   eax, 16
		xor   eax, edx
		movsx ebx, byte ptr [esi+2]
		shl   ebx, 18
		xor   eax, ebx
		mov   edx, eax
		shr   eax, 11
		add   eax, edx
		jmp   L3

		lcase2:
		mov    bx, word ptr [esi]
		add   eax, ebx
		mov   edx, eax
		shl   eax, 11
		xor   eax, edx
		mov   edx, eax
		shr   eax, 17
		add   eax, edx
		jmp   L3

		lcase1:
		movsx ebx, byte ptr [esi]
		add   eax, ebx
		mov   edx, eax
		shl   eax, 10
		xor   eax, edx
		mov   edx, eax
		shr   eax, 1
		add   eax, edx

		L3:

		mov   ebx, eax
		shl   eax, 3
		xor   eax, ebx
		mov   ebx, eax
		shr   eax, 5
		add   eax, ebx
		mov   ebx, eax
		shl   eax, 4
		xor   eax, ebx
		mov   ebx, eax
		shr   eax, 17
		add   eax, ebx
		mov   ebx, eax
		shl   eax, 25
		xor   eax, ebx
		mov   ebx, eax
		shr   eax, 6
		add   eax, ebx

		mov   hash, eax
	}
#endif

	return hash;
}

/* ======================================================================== */

/* A hashing function that I believe was created by either Chris Torek or 
   Dan Bernstein */

uint32_t alphaNumHash (const char * s, int len) {
uint32_t h; 
int i;

    for (h = 0, i = 0; i < len; i++) {
        h = (h << 5) + (h * 5) + (unsigned char) s[i];
    }
    return h;
} 

/* ======================================================================== */

uint32_t foldBits (uint32_t k, uint32_t n) {
uint32_t p = 0;

	if ((n + ~UINT32_C(0)) >= 31) return k;
	while (k) {
		p ^= k;
		k = k >> n;
	}
	return p & ((1 << n) - 1);
}

typedef uint32_t (* hashFn) (const char * s, int len);

#define BUFF_SZ (128*2)
#define NTESTS INT32_C(5000000)

double test (hashFn hash) {
static char buff[BUFF_SZ];
clock_t c0, c1;
int32_t i;

	for (buff[0]=0, i=1; i < BUFF_SZ; i++) buff[i] = (char) (i + buff[i-1]);

	c0 = clock ();
	for (i=0; i < NTESTS; i++) hash (buff, BUFF_SZ);
	c1 = clock ();
	return (c1 - c0)*(1.0 / (double)CLOCKS_PER_SEC);
}

struct tagtest {
	double res;
	char * name;
	hashFn hash;
} tests[] = {
#if 0
	{ 0.0, "CRC16\t\t", GetCRC16                  },
	{ 0.0, "oneAtATimeHashPH", oneAtATimeHashPH   },
	{ 0.0, "CRC32PH\t\t", GetCRC32PH              },
	{ 0.0, "SFHAsm\t\t", SuperFastHashAsm         },
#endif
	{ 0.0, "CRC32\t\t", GetCRC32                  },
	{ 0.0, "oneAtATimeHash\t", oneAtATimeHash     },
	{ 0.0, "alphaNumHash\t",  alphaNumHash        },
	{ 0.0, "FNVHash\t\t", FNVHash                 },
	{ 0.0, "BobJenkins\t", hashBobJenkins         },
	{ 0.0, "SuperFastHash\t", SuperFastHash       },
	{ 0.0, NULL, NULL                             }
};

int main () {
int i, j;
	GenerateCRC32Table ();

	for (j=0; tests[j].name != NULL; j++) {
		for (i=0; i < 3; i++) {
			double res = test (tests[j].hash);
			if (tests[j].res == 0.0 || tests[j].res > res) tests[j].res = res;
		}
		printf ("%s:%8.4fs\n", tests[j].name, tests[j].res);
	}

	return 0;
}
