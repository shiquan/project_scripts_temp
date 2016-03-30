/* Construct common sequences from multi alignment results.
 * All sequences in fasta file should be equally length.
 * 
 * 2016 Copyright by Shi Quan (shiquan@genomics.cn)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/wait.h>
#include <unistd.h>
#include <signal.h>
#include <assert.h>
#include <fcntl.h>
#include <pthread.h>
#include <htslib/faidx.h>

#define ERR_PROC_FINISHED  0
#define ERR_OS_PLATFORM   -1
#define ERR_MISSING_INFO  -2

/* memory stat 
 * This part is adapt from Li Heng's runit program.
 */

// static system information
struct run_sys_static {
    size_t mem_total;
    size_t page_size;
    size_t swap_total;
};
// dynamic system information
struct run_sys_dyn {
    double wall_clock;
    size_t mem_free;
    size_t mem_available;
};
// dynamic process information
struct run_proc_dyn {
    size_t rss;
    size_t vsize;
    double utime;
    double stime;
};

int get_static_sys_info(struct run_sys_static *rss);
int get_dyn_sys_info(struct run_sys_dyn *rsd);
int get_dyn_proc_info(pid_t pid, struct run_proc_dyn *rpd);

#ifdef __linux__
int get_static_sys_info(struct run_sys_static *rss)
{
    FILE *fp = NULL;
    char buffer[64];
    unsigned flag = 0;
    rss->page_size = sysconf(_SC_PAGESIZE);
    fp = fopen("/proc/meminfo", "r");
    if (fp == NULL) return ERR_OS_PLATFORM;

    while (fscanf(fp, "%s", buffer) > 0) {
	if (strstr(buffer, "MemTotal") == buffer) {
	    fscanf(fp, "%lu", &(rss->mem_total));
	    rss->mem_total *= 1024;
	    flag |= 0x1;	    
	} else if (strstr(buffer, "SwapTotal") == buffer) {
	    fscanf(fp, "%lu", &(rss->swap_total));
	    rss->swap_total *= 1024;
	    flag |= 0x2;
	}
    }
    fclose(fp);
    return (flag == 0x3) ? 0 : ERR_MISSING_INFO;
}

int get_dyn_sys_info(struct run_sys_dyn *rsd)
{
    FILE *fp = NULL;
    size_t mem_buffer, mem_cache;
    char buffer[64];
    unsigned flag = 0;
    struct timeval tp;
    struct timezone tzp;

    gettimeofday(&tp, &tzp);
    rsd->wall_clock = tp.tv_sec + tp.tv_usec*1e-6;
    fp = fopen("/proc/meminfo", "r");
    if (fp == NULL) return ERR_OS_PLATFORM;
    while (fscanf(fp, "%s", buffer)>0) {
	if (strstr(buffer, "MemFree") == buffer) {
	    fscanf(fp, "%lu", &(rsd->mem_free));
	    flag |= 0x1;
	} else if (strstr(buffer, "Buffers") == buffer){
	    fscanf(fp, "%lu", &mem_buffer);
	    flag |= 0x2;
	} else if (strstr(buffer, "Cached") == buffer) {
	    fscanf(fp, "%lu", &mem_cache);
	    flag |= 0x4;
	}
    }
    rsd->mem_available = rsd->mem_free + mem_buffer + mem_cache;
    rsd->mem_free *= 1024;
    rsd->mem_available *= 1024;
    fclose(fp);
    return (flag == 0x7) ? 0 : ERR_MISSING_INFO;
}

int get_dyn_proc_info(pid_t pid, struct run_proc_dyn *rpd)
{
    int c, n_spc;
    char str[64];
    FILE *fp = NULL;
    unsigned long tmp, tmp2;
    size_t page_size;

    page_size = sysconf(_SC_PAGESIZE_);
    sprintf(str, "/proc/%u/stat",pid);
    fp = fopen(str, "r");
    if (fp == NULL) return ERR_PROC_FINISHED;
    n_spc = 0;
    while((c = fgetc(fp)) != EOF) {
	if (c == ' ') ++n_spc;
	if (n_spc == 22) break;
    }
    fscanf(fp, "%lu%lu", &tmp, &tmp2);
    fclose(fp);
    rpd->vsize = tmp/1024;
    rpd->rss = tmp2 *(page_size/1024);
    rpd->rss *= 1024;
    rpd->vsize *= 1024;
    return 0;
}

#endif  /* __linux */

#ifdef __APPLE__

#include <sys/types.h>
#include <sys/sysctl.h>
#include <sys/vmmeter.h>
#include <sys.>

#endif /* __apple */

int main(int argc, char **argv)
{
    struct run_sys_static rss;
    struct run_sys_dyn rsd;
    get_static_sys_info(&rss);
    get_dyn_sys_info(&rsd);
    fprintf(stderr, "totalmem : %16.3f KB\n", rss.mem_total/ 1024.0);
    fprintf(stderr, "available: %16.3f KB\n", rsd.mem_available/1024.0);
    fprintf(stderr, "free     : %16.3f KB\n", rsd.mem_free /1024.0);
    return EXIT_SUCCESS;
}
