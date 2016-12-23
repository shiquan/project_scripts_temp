#include "utils.h"
#include "genepred.h"



#ifdef GENEPRED_TEST_MAIN

struct args {
    const char *input_fname;
} args = {
    .input_fname = 0,
};

int parse_args(int ac, char **av)
{
    return 0;
}
int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    test_genepred();
    memory_release();
    return 0;
}

#endif
