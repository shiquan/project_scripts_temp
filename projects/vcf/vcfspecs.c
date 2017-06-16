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
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <dlfcn.h>
#include <unistd.h> // access() function
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/kstring.h>

#ifndef KSTRING_INIT
#define KSTRING_INIT { 0, 0, 0, }
#endif

#define VCFSPECENV "VCFSPEC"

static int error_flag = 0;

struct lpaths {
    int n, m;
    kstring_t *paths;
};

typedef void (*dl_version_f)(const char **, const char **);
typedef int (*dl_run_f)(int, char **);
typedef int (*dl_init_f)();
struct plugin_spec {
    char *name;
    int argc;
    char **argv;
    // output the version number
    dl_version_f version;
    // init the output header, and argvs
    dl_init_f init;
    // handle each line
    dl_run_f run;
    // print usage information
    dl_usage_f usage;
    // print summary information, this is usually a short introduction
    dl_summary_f summary;
    // release memeory, and export final report
    dl_final_f final;
    
    void *handle;
};
void *dlopen_plugin(const char *path)
{
    return dlopen(path,  RTLD_NOW);
}

// return 1 on failure, return 0 on success
int load_plugin(struct plugin_spec *plugin, const char *path)
{
    plugin->name = strdup(path);
    plugin->handle = dlopen_plugin(path);
    if ( plugin->handle == NULL)
        return 1;

    dlerror();
    plugin->init = (dl_init_f)dlsym(plugin->handle, "init");
    char *ret = dlerror();
    if ( ret )
        error("%s : %s", path, ret);
    plugin->run = (dl_init_f)dlsym(plugin->handle, "run");
    ret = dlerror();
    if ( ret )
        error("%s : %s", path, ret);

    plugin->version = (dl_version_f)dlsym(plugin->handle, "version");
    ret = dlerror();
    if ( ret )
        error("%s : %s", path, ret);

    plugin->usage = (dl_usage_f)dlsym(plugin->handle, "usage");
    ret = dlerror();
    if ( ret )
        error("%s : %s", path, ret);
    
    plugin->summary = (dl_version_f)dlsym(plugin->handle, "summary");
    ret = dlerror();
    if ( ret )
        error("%s : %s", path, ret);

    plugin->run = (dl_run_f)dlsym(plugin->handle, "run");
    ret = dlerror();
    if ( ret )
        error("%s : %s", path, ret);

    plugin->final = (dl_final_f)dlsym(plugin->handle, "final");
    ret = dlerror();    
    if ( ret )
        error("%s : %s", path, ret);    
    
    return 0;
}
// return number of path directroies
struct lpaths *library_path_init(const char *env)
{
    char *path = getenv(env);
    if (path == NULL) {
        error_flag = 1;
        return NULL;
    }

    struct lpaths *lpaths = (struct lpaths *)malloc(sizeof(struct lpaths));
    memset(lpaths, 0, sizeof(struct lpaths));
    
    char *ss = path;
    char *se = NULL;
    do {
        se = strchr(ss, ':');
        if ( lpaths->n == lpaths->m ) {
            lpaths->m = lpaths->n == m ? 1 : lpaths->m += 2;
            lpaths->paths = (kstring_t*)realloc(lpaths->paths, lpaths->m * sizeof(kstring_t));
        }
        kstring_t *str = &lpaths->paths[lpaths->n-1];
        memset(str, 0, sizeof(kstring_t));
        kputsn(ss, se -ss, str);
        ss = se + 1;
        lpaths->n ++;
    } while (se != NULL);
    return lpaths;
}

int show_plugins()
{
    kstring_t string = KSTRING_INIT;
    struct lpaths * lpaths = library_path_init(VCFSPECENV);
    int i;
    for (i = 0; i < lpaths->n; ++i ) {
        DIR *d = opendir(lpaths->paths[i].s);
        if (d == NULL)
            continue;
        struct dirent *e;
        while ( (e = readdir(d)) ) {
            int length = strlen(e->d_name);
            if ( strcasecmp(".so", e->d_name+length-3) )
                continue;
            string.l = 0;
            ksprintf(&string, "%s/%s", lpaths->paths[i].s, e->d_name);
            fprintf(stderr, "%s\n", string.s);
        }
        closedir(d);
    }
    free(string.s);
    return 0;
}

int usage()
{
    fprintf(stderr,
            "Usage:\n"
            "  vcfspecs <plugins>  [ options ]\n"
            "plugins:\n");

    show_plugins();

    fprintf(stderr,
            "Homepage:\n"
            "https://github.com/shiquan/small_projects_collections\n");
    return 1;
}


char *find_main_func(char *name)
{
    
}

int parse_args(int argc, char **argv)
{
    struct lpaths *paths = library_path_init(VCFSPECENV);
    if ( paths == NULL) 
        error("No $VCFSPEC environment set.");
    
    if ( argc == 1)
        return usage();

    char *dyn = find_main_func(argv[1]);
    if ( dyn )
        error("Failed to load %s", dyn);

    void *handle;
    handle = dl_open(dyn, RTLD_NOW);
    free(dyn);
    
    return 0;
}

void memory_release()
{
}

int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    memory_release();
    
    return 0;
}
