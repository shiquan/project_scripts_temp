#ifndef OPTIONS_HEADER
#define OPTIONS_HEADER

#include <stdio.h>
#include <stdlib.h>


enum argument_type {
    _argument_type_promote_to_int = -1,
    argument_bool, // no argument required
    argument_string,
    argument_int,
    argument_float,
    argument_path,
    argument_string_or_path,
    _unknown,
};

enum argument_required {
    argument_must_set,
    argument_optional,
    _unknown,
};

struct option_spec {
    char *name;
    enum argument_type argument_type;
    enum argument_required argument_required;
    char *comment;
};

// print verbose information for arguments
void set_verbose_mode();

// init arguments
int init_args(struct option_spec *opts, int argc, char **argv);

// argument exists, 1 for yes, 0 for empty
int argument_exist(char *name);

// get the value of arguments
// get_argument_string() for string
// get_argument_path() for file path or directory
// get_argument_int() for int value
// get_argument_float() for float value
char *get_argument_string(char *name);
int get_argument_int(char *name);
float get_argument_float(char *name);

// release the memory
void arguments_destroy();


#endif

