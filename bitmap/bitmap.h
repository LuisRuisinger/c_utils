#ifndef C_UTILS_BITMAP_H
#define C_UTILS_BITMAP_H

#include <stdlib.h>
#include <string.h>

#include "../defines.h"

C_GUARD_BEGINN()

typedef u32 BITMAP_DWORD;

#define BITMAP_DWORD_SIZE (sizeof(BITMAP_DWORD) * 8)

// amount of dwords needed to represent the cnt bits in a bitmap
#define BITMAP_N_DWORDS(cnt) \
    (((cnt) + BITMAP_DWORD_SIZE - 1) / BITMAP_DWORD_SIZE)

typedef struct Bitmap_t {
    BITMAP_DWORD dword;
} Bitmap;

#define BITMAP(name, cnt) \
    Bitmap (name)[BITMAP_N_DWORDS(cnt)]

#define BITMAP_DWORD(bitmap, bit) \
    ((bitmap)[(bit) / BITMAP_DWORD_SIZE].dword)

#define BITMAP_DWORD_BIT(bit) \
    (1UL << ((bit) & (BITMAP_DWORD_SIZE - 1)))

C_UTILS_INLINE static inline usize bitmap_sizeof(u32 cnt) {
    return BITMAP_N_DWORDS(cnt) * sizeof(BITMAP_DWORD);
}

C_UTILS_INLINE static inline void bitmap_set_bit(Bitmap *bitmap, u32 n) {
    BITMAP_DWORD(bitmap, n) |= BITMAP_DWORD_BIT(n);
}

C_UTILS_INLINE static inline void bitmap_set_bits(Bitmap *bitmap, u32 n, u32 cnt) {
    for (usize i = 0; i < cnt; ++i) {
        BITMAP_DWORD(bitmap, n + i) |= BITMAP_DWORD_BIT(n + i);
    }
}

C_UTILS_INLINE static inline void bitmap_clear_bit(Bitmap *bitmap, u32 n) {
    BITMAP_DWORD(bitmap, n) &= ~BITMAP_DWORD_BIT(n);
}

C_UTILS_INLINE static inline void bitmap_flip_bit(Bitmap *bitmap, u32 n) {
    BITMAP_DWORD(bitmap, n) ^= BITMAP_DWORD_BIT(n);
}

C_UTILS_INLINE static inline bool bitmap_test_bit(Bitmap *bitmap, u32 n) {
    return (bool) (BITMAP_DWORD(bitmap, n) & BITMAP_DWORD_BIT(n));
}

C_UTILS_INLINE static inline void bitmap_zero(Bitmap *bitmap, u32 cnt) {
    memset(bitmap, 0, bitmap_sizeof(cnt));
}

C_UTILS_INLINE static inline void bitmap_fill(Bitmap *bitmap, u32 cnt) {
    memset(bitmap, UINT8_MAX, bitmap_sizeof(cnt));
}

C_UTILS_INLINE static inline void bitmap_cpy(Bitmap *dst, Bitmap *src, u32 cnt) {
    memcpy(dst, src, bitmap_sizeof(cnt));
}

C_UTILS_INLINE static inline void bitmap_complement(Bitmap *dst, Bitmap *src, u32 cnt) {
    for (usize i = 0; i < BITMAP_N_DWORDS(cnt); ++i) {
            dst[i].dword = ~src[i].dword;
    }
}

C_UTILS_INLINE static inline u32 bitmap_popcount(Bitmap *bitmap, u32 cnt) {
    u32 popcnt = 0;

    // TODO: abstract popcount (no guarantee this builtin exists)
    for (usize i = 0; i < BITMAP_N_DWORDS(cnt); ++i) {
        popcnt += __builtin_popcount(bitmap[i].dword);
    }

    return popcnt;
}

C_UTILS_INLINE static inline u32 bitmap_clz(Bitmap *bitmap, u32 cnt) {
    u32 clz = 0;

    for (usize i = 0; i < BITMAP_N_DWORDS(cnt); ++i) {
        clz += CPU_RAYTRACING_CLZ(bitmap[i].dword);

        if (clz & (BITMAP_DWORD_SIZE - 1)) {
            return clz;
        }
    }

    return clz;
}

#define BITMAP_DEFINE_BINOP(_name, _op)                                                           \
    C_UTILS_INLINE static inline void bitmap_##_name(                                              \
        Bitmap *dst, const Bitmap *a, const Bitmap *b, u32 cnt) {                                 \
        for (usize i = 0; i < BITMAP_N_DWORDS(cnt); ++i) {                                        \
            dst[i].dword = a[i].dword _op b[i].dword;                                             \
        }                                                                                         \
    }

BITMAP_DEFINE_BINOP(and, &)
BITMAP_DEFINE_BINOP(nand, & ~)
BITMAP_DEFINE_BINOP(or, |)
BITMAP_DEFINE_BINOP(nor, | ~)
BITMAP_DEFINE_BINOP(xor, ^)
BITMAP_DEFINE_BINOP(nxor, ^ ~)

// comparison
BITMAP_DEFINE_BINOP(lt, <)
BITMAP_DEFINE_BINOP(eq, ==)
BITMAP_DEFINE_BINOP(neq, !=)
BITMAP_DEFINE_BINOP(gt, >)
BITMAP_DEFINE_BINOP(le, <=)
BITMAP_DEFINE_BINOP(ge, >=)

C_GUARD_END()

#endif //C_UTILS
