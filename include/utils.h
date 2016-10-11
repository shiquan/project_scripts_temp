#ifndef UTILS_COMMON_HEADER
#define UTILS_COMMON_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>
#include <time.h>

#if defined(_M_AMD_64) || defined(__x86_64__) || defined(__arm64__)
#define OS_BITS_64  1
#else
#define OS_BITS_64  0
#endif

#if defined(_MSC_VER) && !defined(__clang__)
# define FORCE_INLINE __forceinline
# define inline __inline
#else
# define FORCE_INLINE static inline __attribute__((__always_inline__))
#endif


#define check_mem(p) do {\
        void **_pp = (void**)&p;\
        if (_pp == NULL || *_pp == NULL) {\
            fprintf(stderr, "[memory out] func: %s, line: %d\n", __FUNCTION__, __LINE__);\
            exit(EXIT_FAILURE);\
        }\
    } while(0)

#define str_errno() (errno == 0 ? "None" : strerror(errno))

#define clear_errno() do \
    {\
        fprintf(stderr, "%s\n", str_errno());\
        errno = 0;\
    } while(0)

#define error(line, ...) do \
    {\
        fprintf(stderr, "[error] [func: %s, line: %d] " line "\n", __FUNCTION__, __LINE__, ##__VA_ARGS__);\
        errno = 0;\
        exit(EXIT_FAILURE); \
    } while(0)

#define error_print(line, ...) do \
    { \
        fprintf(stderr, "[error] [func: %s, line: %d] " line "\n", __FUNCTION__, __LINE__, ##__VA_ARGS__); \
    } while(0)

#define warnings(line, ...) do \
    { \
        fprintf(stderr, "[warnings] " line "\n", ##__VA_ARGS__);\
    } while(0)

#define debug_print(line, ...) do {\
        fprintf(stderr, "[ ** DEBUG ** func: %s, line: %d ] " line "\n", __FUNCTION__, __LINE__, ##__VA_ARGS__); \
    } while(0)

#define LOG_print(line, ...) do {\
        time_t second;\
        time(&second);\
        char _time_buff[100];\
        strftime (_time_buff, 100, "%Y-%m-%d %H:%M:%S", localtime (&second)); \
        fprintf(stderr, "[%s] " line "\n", _time_buff, ##__VA_ARGS__);\
    } while(0)

#define BE_SMART_STRING "Please DONOT post this error message on the forum or copy it into the emails."\
    "Try to figure out this issue by youself by reading the log information carefully and checking you input arguments."

#define CHECK_STDIN (!isatty(fileno(stdin)))

#endif
