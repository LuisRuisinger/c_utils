//
// Created by Luis Ruisinger on 26.11.24.
//

#ifndef C_UTILS_AABB_H
#define C_UTILS_AABB_H

#include "../linalg/vec3.h"

C_GUARD_BEGINN()

typedef struct AABB_t {
    vec3f min;
    vec3f max;
} AABB;

#define AABB_1(v_max) \
    (AABB) { .min = VEC3(0.0F), .max = v_max }

#define AABB_2(v_min, v_max) \
    (AABB) { .min = v_min, .max = v_max }

#define AABB_GET_MACRO(_1, _2, NAME, ...) NAME

#define AABB(...) \
    AABB_GET_MACRO(__VA_ARGS__, AABB_2, AABB_1)(__VA_ARGS__)

C_UTILS_INLINE static inline void scale_vec(AABB *aabb, vec3f vec) {
    vec3f tmp = vec3_sub(aabb->max, aabb->min);
    aabb->max = vec3_mul(tmp, vec);
}

C_UTILS_INLINE static inline void scale_scalar(AABB *aabb, f32 scalar) {
    vec3f vec = VEC3(scalar, scalar, scalar);
    scale_vec(aabb, vec);
}

C_UTILS_INLINE static inline void scale_center_vec(AABB *aabb, vec3f vec) {
    vec3f tmp_0 = vec3_add(aabb->max, aabb->min);
    tmp_0 = vec3_muls(tmp_0, 0.5F);

    vec3f tmp_1 = vec3_sub(aabb->max, aabb->min);
    tmp_1 = vec3_muls(tmp_1, 0.5F);
    tmp_1 = vec3_mul(tmp_1, vec);

    aabb->min = vec3_sub(tmp_0, tmp_1);
    aabb->max = vec3_add(tmp_0, tmp_1);
}

C_UTILS_INLINE static inline void scale_center_scalar(AABB *aabb, f32 scalar) {
    vec3f vec = VEC3(scalar, scalar, scalar);
    scale_center_vec(aabb, vec);
}

C_UTILS_INLINE static inline void translate_vec(AABB *aabb, vec3f vec) {
    aabb->min = vec3_add(aabb->min, vec);
    aabb->max = vec3_add(aabb->max, vec);
}

C_UTILS_INLINE static inline void translate_scalar(AABB *aabb, f32 scalar) {
    vec3f vec = VEC3(scalar);
    translate_vec(aabb, vec);
}

C_UTILS_INLINE static inline bool aabb_aabb_intersection(const AABB *a, const AABB *b) {
    __m128 tmp_0 = _mm_cmple_ps(a->min.vec, a->max.vec);
    __m128 tmp_1 = _mm_cmpge_ps(a->max.vec, b->min.vec);
    tmp_1 = _mm_and_si128(tmp_0, tmp_1);

    return !!( _mm_test_all_ones(tmp_1));
}

C_UTILS_INLINE static inline void aabb_grow_vec(AABB *aabb, vec3f v) {
    aabb->min = vec3_min(aabb->min, v);
    aabb->max = vec3_max(aabb->max, v);
}

void aabb_grow_aabb(AABB *aabb, const AABB *other) {
    aabb->min = vec3_min(aabb->min, other->min);
    aabb->max = vec3_max(aabb->max, other->max);
}

C_UTILS_INLINE static inline f32 aabb_area(const AABB *aabb) {
    vec3f v = vec3_sub(aabb->max, aabb->min);

    __m128 _tmp_0 = _mm_shuffle_epi32(_mm_castps_si128(v.vec), _MM_SHUFFLE(3, 1, 2, 0));;
    vec3f w = (vec3f) { .vec = _mm_castsi128_ps(_tmp_0) };

    return vec3_hsum(vec3_mul(v, w));
}

C_GUARD_END()

#endif //SOFTWARE_RATYTRACING_AABB_H
