#ifndef _ELEMENTARY_H_
#define _ELEMENTARY_H_

#include <stdint.h>
#include <stdbool.h>

typedef unsigned int uint;

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef int8_t s8;
typedef int16_t s16;
typedef int32_t s32;
typedef int64_t s64;


#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define BOUND(x, a, b) ((x) < (a) ? (a) : ((b) < (x) ? (b) : (x)))

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#endif // !_ELEMENTARY_H_
