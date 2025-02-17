#ifndef C_UTILS_DEFINES_H
#define C_UTILS_DEFINES_H

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <stdalign.h>

#include "build.h"
#include "builtin.h"

C_GUARD_BEGINN()

#define MAX(x, y) \
    (((x) > (y)) ? (x) : (y))

#define MIN(x, y) \
    (((x) < (y)) ? (x) : (y))

#define TO_RADIANS(deg) \
    ((deg) * M_PI / 180.0F)

#define TO_DEGREES(rad) \
    ((rad) * 180.0F / M_PI)

#define FAST_MOD_POW2(num, mod) \
    ((num) & (mod - 1))

#define RND_POW_2(_n, _p) \
    (((_n) + (_p) - 1) & ~((_p) - 1))

#define BITWISE_U32(_type, _t) \
    ((union { u32 u; _type t; }) { .t = _t }.u)

#define BITWISE_TYPE(_type, _u) \
    ((union { u32 u; _type t; }) { .u = _u }.t)

C_GUARD_END()

#endif //C_UTILS
