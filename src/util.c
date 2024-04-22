#include "util.h"

void* __util_malloc_check_exit(size_t size, const char* file, uint line) {
        void* ptr = malloc(size);
        if (ptr == NULL) {
                fprintf(stderr, "# ERROR: Failed to malloc, file %s, line %u.\n", file, line);
                exit(EXIT_FAILURE);
        }
        return ptr;
}

void* __util_realloc_check_exit(void* ptr, size_t size, const char* file, uint line) {
        ptr = realloc(ptr, size);
        if (ptr == NULL) {
                fprintf(stderr, "# ERROR: Failed to realloc, file %s, line %u.\n", file, line);
                exit(EXIT_FAILURE);
        }
        return ptr;
}


FILE* __util_fopen_check(const char* filename, const char* mode, const char* f, uint line) {
	FILE* ret = fopen(filename, mode);
	if (ret == NULL) {
		fprintf(stderr, "# INFO: file %s, line %u.\n", f, line);
		fprintf(stderr, "  Failed trying to open file %s.\n", filename);
	}
	return ret;
}

void __util_fclose(FILE* file) {
	fclose(file);
}

void __util_fwrite_check(const void* ptr, size_t size, size_t num, FILE* file, const char* f, uint line) {
	const size_t w = fwrite(ptr, size, num, file);
	if (w != num) {
		fprintf(stderr, "# INFO: file %s, line %u.", f, line);
		fprintf(stderr, "  Failed to read %lu bytes, read %lu bytes.\n", size*num, w*size);
	}
}

void __util_fread_check(void* ptr, size_t size, size_t num, FILE* file, const char* f, uint line) {
	const size_t read = fread(ptr, size, num, file);
	if (read != num) {
		fprintf(stderr, "# INFO: file %s, line %u.", f, line);
		fprintf(stderr, "  Failed to read %lu bytes, read %lu bytes.\n", size*num, read*size);
	}
}
