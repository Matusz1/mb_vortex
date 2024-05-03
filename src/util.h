#ifndef _UTIL_H_
#define _UTIL_H_

#include <stdio.h>
#include <stdlib.h>

#include "elementary.h"

#define util_error(cond, ...) {                                                     \
        if (cond) {                                                                 \
                fprintf(stderr, "# ERROR: file %s, line %u\n", __FILE__, __LINE__); \
                fprintf(stderr, __VA_ARGS__);                                       \
                exit(EXIT_FAILURE);                                                 \
        }                                                                           \
}

#define util_malloc(size)  __util_malloc_check_exit(size, __FILE__, __LINE__)
#define util_realloc(ptr, size) __util_realloc_check_exit(ptr, size, __FILE__, __LINE__)
#define util_free(file) free(file)

#define util_fopen(fname, mode) __util_fopen_check(fname, mode, __FILE__, __LINE__)
#define util_fclose(file) __util_fclose(file)
#define util_fwrite(ptr, size, num, file) __util_fwrite_check(ptr, size, num, file, __FILE__, __LINE__)
#define util_fread(ptr, size, num, file) __util_fread_check(ptr, size, num, file, __FILE__, __LINE__)


void* __util_malloc_check_exit(size_t size, const char* file, uint line);
void* __util_realloc_check_exit(void* ptr, size_t size, const char* file, uint line);

FILE* __util_fopen_check(const char* filename, const char* mode, const char* f, uint line);
void __util_fclose(FILE* file);
void __util_fwrite_check(const void* ptr, size_t size, size_t num, FILE* file, const char* f, uint line);
void __util_fread_check(void* ptr, size_t size, size_t num, FILE* file, const char* f, uint line);

#endif // !_UTIL_H_
