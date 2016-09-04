#include "abi.h"
#include "utils.h"


struct abi_read *abi_read_init()
{
}

void abi_read_destroy(struct abi_read *read)
{
}
static int check_magic_key(FILE *fp)
{
    rewind(fp);
    uint32_t magic;
    
}
struct abi_read *abi_read_core(FILE *fp)
{
    abi_read *abi_read = abi_read_init();
    
    if ( check_magic_key(fp) == 1) {
	error_print("This is not a AB1 file.");
	goto bail_out;
    }


    return abi_read;
    
  bail_out:
    abi_read_destroy(abi_read);
    return NULL;
}

struct abi_read *abi_read(const char *fn)
{
    FILE *fp;
    fp = fopen(fn, "rb");
    if (fp == NULL)
	error("Failed to open %s : %s.", fn, strerror(errno));
    
    struct abi_read *abi_read = abi_read_core(fp);
    fclose(fp);
    
    if ( abi_read == NULL )
	error("Failed to read %s.", fn);

    return abi_read;	
}

