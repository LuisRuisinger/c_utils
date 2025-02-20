//
// Created by Luis Ruisinger on 27.12.24.
//

#ifndef C_UTILS_BUILTIN_H
#define C_UTILS_BUILTIN_H

#include <stdbool.h>
#include <stdio.h>

#include "build.h"

C_GUARD_BEGINN()

typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef int8_t   i8;
typedef int16_t  i16;
typedef int32_t  i32;
typedef int64_t  i64;

typedef float_t  f32;
typedef double_t f64;

typedef size_t   usize;
typedef ssize_t  ssize;

#if defined(C_UTILS_X86)
#if defined(__AVX2__) || defined(__AVX__)
#if defined(__AVX512F__)
    #define C_UTILS_AVX512F
    #define ALIGNMENT_512 sizeof(__m512)
#endif
    #include <immintrin.h>
    #define C_UTILS_AVX
    #define ALIGNMENT_256 sizeof(__m256)
#endif

#if defined(__SSE4_2__)
    #include <nmmintrin.h>
    #define C_UTILS_SSE4_2
#endif

#if defined(__SSE4_1__)
    #include <smmintrin.h>
    #define C_UTILS_SSE4_2
#endif

#if defined(__SSSE3__)
    #include <tmmintrin.h>
    #define C_UTILS_SSSE3

    // src https://stackoverflow.com/a/35270026
    C_UTILS_INLINE static inline f32 _mm_hsum_ps(__m128 _a) {
        __m128 _tmp_0 = _mm_movehdup_ps(_a);
        __m128 _tmp_1 = _mm_add_ps(_a, _tmp_0);
        __m128 _tmp_2 = _mm_movehl_ps(_tmp_0, _tmp_1);
        _tmp_2 = _mm_add_ss(_tmp_1, _tmp_2);

        return _mm_cvtss_f32(_tmp_2);
    }
#endif

#if defined(__SSE3__)
    #include <pmmintrin.h>
    #define C_UTILS_SSE3
#endif

#if defined(__SSE2__)
    #include <emmintrin.h>
    #define C_UTILS_SSE2
    #define ALIGNMENT_128 sizeof(__m128)
#endif
#endif

#if defined(__has_builtin)
    #define C_UTILS_HAS_BUILTIN(intrin) __has_builtin(intrin)
#elif
    #define C_UTILS_HAS_BUILTIN(intrin) true
#endif

// unreachable
#if defined(C_UTILS_GNUC_CLANG) && C_UTILS_HAS_BUILTIN(__builtin_unreachable)
    #define C_UTILS_UNREACHABLE() __builtin_unreachable()
#elif defined(C_UTILS_MSC)
    #define C_UTILS_UNREACHABLE(stmt) __assume(false)
#else
    #define C_UTILS_UNREACHABLE()
#endif

// assume
#if defined(C_UTILS_GNUC_CLANG) && C_UTILS_HAS_BUILTIN(__builtin_assume)
    #define C_UTILS_ASSUME(stmt) __builtin_assume(stmt)
#elif defined(C_UTILS_MSC)
    #define C_UTILS_ASSUME(stmt) __assume(stmt)
#else
    #define C_UTILS_ASSUME(stmt)
#endif

#if defined(C_UTILS_GNUC_CLANG) && C_UTILS_HAS_BUILTIN(__builtin_bit_cast)
    #define C_UTILS_BITCAST(_t, _f) __builtin_bit_cast(_t, _f)
#endif

#if defined(C_UTILS_GNUC_CLANG) && C_UTILS_HAS_BUILTIN(__builtin_assume_aligned)
    #define C_UTILS_ASSUME_ALIGNED(_exp, _align, ...) \
        __builtin_assume_aligned(_exp, _align, ##__VA_ARGS__)
#else 
    #define C_UTILS_ASSUME_ALIGNED(_exp, _align, ...) _exp
#endif

#if defined(C_UTILS_GNUC_CLANG) && \
    C_UTILS_HAS_BUILTIN(__builtin_clz) && C_UTILS_HAS_BUILTIN(__builtin_clzll)

    C_UTILS_INLINE static inline i32 cr_clz_u8(u8 x) {
        return __builtin_clz(x) - (i32) ((sizeof(u32) - sizeof(u8)) * 8);
    }

    C_UTILS_INLINE static inline i32 cr_clz_u16(u16 x) {
        return __builtin_clz(x) - (i32) ((sizeof(u32) - sizeof(u16)) * 8);
    }

    C_UTILS_INLINE static inline i32 cr_clz_u32(u32 x) {
        return __builtin_clz(x);
    }

    C_UTILS_INLINE static inline i32 cr_clz_u64(u64 x) {
        return __builtin_clzll(x);
    }

#else
    C_UTILS_INLINE static inline i32 clz(u32 x) {
        union c {
            u32 u[2]; f64 f;
        } cast = { .f = (f32) x };

        i32 tmp = ((cast.u[!((union c){.f = 1.0F}.u[0])] >> 20) & 0x7FF) - 1023;

        return ((i32) sizeof(u32) * __CHAR_BIT__) - tmp + 1;
    }

    C_UTILS_INLINE static inline i32 cr_clz_u8(u8 x) {
        return clz(x) - (i32) ((sizeof(u32) - sizeof(u8)) * 8);
    }

    C_UTILS_INLINE static inline i32 cr_clz_u16(u16 x) {
        return clz(x) - (i32) ((sizeof(u32) - sizeof(u16)) * 8);
    }

    C_UTILS_INLINE static inline i32 cr_clz_u32(u32 x) {
        return clz(x);
    }

    C_UTILS_INLINE static inline i32 cr_clz_u64(u64 x) {
        i32 upper = clz(x >> 32);
        return upper ? upper : clz((u32) x);
    }
#endif

// TODO: add support for other compilers ?

#if defined(C_UTILS_C11_LEAST)
    #define C_UTILS_CLZ(x) \
        _Generic((x), u8: cr_clz_u8, u16: cr_clz_u16, u32: cr_clz_u32, u64: cr_clz_u64)(x)
#endif

#if defined(C_UTILS_GNUC_CLANG) && C_UTILS_HAS_BUILTIN(__builtin_expect)
    #define LIKELY(stmt)   __builtin_expect(stmt, 1)
    #define UNLIKELY(stmt) __builtin_expect(stmt, 0)
#else
    #define LIKELY(stmt)   (stmt)
    #define UNLIKELY(stmt) (!(stmt))
#endif



C_GUARD_END()

#endif //C_UTILS_BUILTIN_H
