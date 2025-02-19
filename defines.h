#ifndef C_UTILS_DEFINES_H
#define C_UTILS_DEFINES_H

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <stdalign.h>
#include <signal.h>

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

#define TRAP() \
    raise(SIGTRAP)

#define __NARG__(...) \
     __NARG_I_(__VA_ARGS__,__RSEQ_N())

#define __NARG_I_(...) \
    __ARG_N(__VA_ARGS__)

#define __ARG_N(                                                                                    \
      _1, _2, _3, _4, _5, _6, _7, _8, _9,_10,                                                       \
     _11,_12,_13,_14,_15,_16,_17,_18,_19,_20,                                                       \
     _21,_22,_23,_24,_25,_26,_27,_28,_29,_30,                                                       \
     _31,_32,_33,_34,_35,_36,_37,_38,_39,_40,                                                       \
     _41,_42,_43,_44,_45,_46,_47,_48,_49,_50,                                                       \
     _51,_52,_53,_54,_55,_56,_57,_58,_59,_60,                                                       \
     _61,_62,_63,N,...) N

#define __RSEQ_N()                                                                                  \
     63, 62, 61, 60,                                                                                \
     59, 58, 57, 56, 55, 54, 53, 52, 51, 50,                                                        \
     49, 48, 47, 46, 45, 44, 43, 42, 41, 40,                                                        \
     39, 38, 37, 36, 35, 34, 33, 32, 31, 30,                                                        \
     29, 28, 27, 26, 25, 24, 23, 22, 21, 20,                                                        \
     19, 18, 17, 16, 15, 14, 13, 12, 11, 10,                                                        \
     9,  8,  7,  6,  5,  4,  3,  2,  1,  0

#define __VFUNC_(name, n) name_##n

#define __VFUNC(name, n) \
    __VFUNC_(name, n)

#define VFUNC(func, ...) \ 
    __VFUNC(func, __NARG__(__VA_ARGS__)) (__VA_ARGS__)

C_GUARD_END()

#endif //C_UTILS
