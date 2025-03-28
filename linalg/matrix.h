#ifndef C_UTILS_MATRIX_H
#define C_UTILS_MATRIX_H

#include <string.h>

#include "../defines.h"
#include "vec3.h"
#include "vec2.h"
#include "vec.h"

C_GUARD_BEGINN()

#if defined(C_UTILS_AVX512F)
    #define MAT_ALIGNMENT 64
#elif defined(C_UTILS_AVX)
    #define MAT_ALIGNMENT 32
#else
    #define MAT_ALIGNMENT 16
#endif





/**
 * @brief Definition of a 2x2 matrix. Struct itself stores the entire matrix. 
 *        Underlying storage is aligned to 128 bit SIMD register.  
 */
typedef struct Mat2x2_t {
    union {
        alignas(16) f32 val[4];
        __m128 row;
    };
} __attribute__((aligned(16))) Mat2x2;





/**
 * @brief Definition of a 4x4 matrix. Struct itself stores the entire matrix. 
 *        Underlying storage is aligned to the biggest available SIMD register of the ISA. 
 */
typedef struct Mat4x4_t {
    union {
        alignas(MAT_ALIGNMENT) f32 val[16];
        alignas(MAT_ALIGNMENT) __m128 rows[4];
    };
} __attribute__((aligned(MAT_ALIGNMENT))) Mat4x4;





#if defined(C_UTILS_AVX512F)
#define ROW_512(mat, i) \
    *((__m512*) &(mat)->val + i)
#endif

#if defined(C_UTILS_AVX)
#define ROW_256(mat, i) \
    *((__m256*) &(mat)->val + i)
#endif

#define ROW_128(mat, i) \
    *((__m128*) &(mat)->val + i)





C_UTILS_INLINE static inline Mat2x2 mat2x2_identity(void) {

    // this does not care about row-major / column-major ordering
    return (Mat2x2) { .row = _mm_set_ps(1.0F, 0.0F, 0.0F, 1.0F) };
}

C_UTILS_INLINE static inline Mat2x2 mat2x2_zeroed(void) {

    // this does not care about row-major / column-major ordering
    return (Mat2x2) { .row = _mm_set1_ps(0.0F) };
}

C_UTILS_INLINE static inline f32 mat2x2_det(const Mat2x2 *__restrict mat) {
    return (mat->val[0] * mat->val[3]) - (mat->val[1] * mat->val[2]);
}

/**
 * @brief  Function to emit a 4x4 identity matrix.
 *         Note: Using memcpy rather than explicit SIMD loads results in poor codegen. 
 * @return 4x4 identity matrix. 
 */
C_UTILS_INLINE static inline Mat4x4 mat4x4_identity(void) {
    const alignas(MAT_ALIGNMENT) f32 arr[16] = {
        0.0F, 0.0F, 1.0F, 0.0F, 0.0F, 0.0F, 0.0F, 1.0F,
        1.0F, 0.0F, 0.0F, 0.0F, 0.0F, 1.0F, 0.0F, 0.0F
    };

#if defined(C_UTILS_AVX512F)
    Mat4x4 mat;

    // this does not care about row-major / column-major ordering
    ROW_512(&mat, 0) = _mm512_load_ps(arr);
    return mat;

#elif defined(C_UTILS_AVX)
    Mat4x4 mat;

    // this does not care about row-major / column-major ordering
    //ROW_256(&mat, 0) = _mm256_set_ps(0.0F, 0.0F, 1.0F, 0.0F, 0.0F, 0.0F, 0.0F, 1.0F);
    //ROW_256(&mat, 1) = _mm256_set_ps(1.0F, 0.0F, 0.0F, 0.0F, 0.0F, 1.0F, 0.0F, 0.0F);

    ROW_256(&mat, 0) = _mm256_load_ps(arr);
    ROW_256(&mat, 1) = _mm256_load_ps(arr + 8);
    return mat;

#else
    Mat4x4 mat;

    // this does not care about row-major / column-major ordering
    ROW_128(&mat, 0) = _mm_load_ps(arr);
    ROW_128(&mat, 1) = _mm_load_ps(arr + 4);
    ROW_128(&mat, 2) = _mm_load_ps(arr + 8);
    ROW_128(&mat, 0) = _mm_load_ps(arr + 12);
    return mat;

#endif
}

/**
 * @brief  Function to emit a zero-ed 4x4 matrix.  
 * @return Zero-ed 4x4 matrix. 
 */
C_UTILS_INLINE static inline Mat4x4 mat4x4_zeroed(void) {
    Mat4x4 mat;
    memset(&mat.val, 0, sizeof(mat.val));

    return mat;
}

/** 
 * @brief  Determinant calculation for a arbitrary 4x4 matrices through the rule of Sarrus applied
 *         on 4x4 matrices. 
 * @param  mat The input matrix for calculating det(mat).    
 * @return det(mat)
 */
C_UTILS_INLINE static inline f32 mat4x4_det(const Mat4x4 *__restrict mat) {
    
    // src: https://stackoverflow.com/a/2980966
    return 
        (mat->val[3] * mat->val[6] * mat->val[9]  * mat->val[12]) - 
        (mat->val[2] * mat->val[7] * mat->val[9]  * mat->val[12]) - 
        (mat->val[3] * mat->val[5] * mat->val[10] * mat->val[12]) + 
        (mat->val[1] * mat->val[7] * mat->val[10] * mat->val[12]) + 
        (mat->val[2] * mat->val[5] * mat->val[11] * mat->val[12]) - 
        (mat->val[1] * mat->val[6] * mat->val[11] * mat->val[12]) - 
        (mat->val[3] * mat->val[6] * mat->val[8]  * mat->val[13]) + 
        (mat->val[2] * mat->val[7] * mat->val[8]  * mat->val[13]) + 
        (mat->val[3] * mat->val[4] * mat->val[10] * mat->val[13]) -
        (mat->val[0] * mat->val[7] * mat->val[10] * mat->val[13]) - 
        (mat->val[2] * mat->val[4] * mat->val[11] * mat->val[13]) + 
        (mat->val[0] * mat->val[6] * mat->val[11] * mat->val[13]) + 
        (mat->val[3] * mat->val[5] * mat->val[8]  * mat->val[14]) - 
        (mat->val[1] * mat->val[7] * mat->val[8]  * mat->val[14]) - 
        (mat->val[3] * mat->val[4] * mat->val[9]  * mat->val[14]) + 
        (mat->val[0] * mat->val[7] * mat->val[9]  * mat->val[14]) + 
        (mat->val[1] * mat->val[4] * mat->val[11] * mat->val[14]) - 
        (mat->val[0] * mat->val[5] * mat->val[11] * mat->val[14]) - 
        (mat->val[2] * mat->val[5] * mat->val[8]  * mat->val[15]) + 
        (mat->val[1] * mat->val[6] * mat->val[8]  * mat->val[15]) + 
        (mat->val[2] * mat->val[4] * mat->val[9]  * mat->val[15]) - 
        (mat->val[0] * mat->val[6] * mat->val[9]  * mat->val[15]) - 
        (mat->val[1] * mat->val[4] * mat->val[10] * mat->val[15]) + 
        (mat->val[0] * mat->val[5] * mat->val[10] * mat->val[15]);
}

/**
 * @brief  Determinant calculation for a 4x4 affine transformation matrices.
 *         Such matrix is composed of
 *         | A   0_v |
 *         | 0_v 1   | 
 *         which allows the laplace extension to det(M) = 1 * det(A) and just discards the other
 *         sub-matrices because of the 0_v's.
 *         The det(M) = det(A) can be calculated using the rule of Sarrus. 
 *         The order of the matrix does not matter as det(M^T) = det(M) and therefore 
 *         | A   0_v |^T        | A^T 0_v |
 *         | 0_v 1   |   equals | 0_v 1   | which yields det(A^T) = det(A)
 * @param  mat The input matrix for calculating det(mat).    
 * @return det(mat)
 */
C_UTILS_INLINE static inline f32 mat4x4_det_tm(const Mat4x4 *__restrict mat) {
    return 
        mat->val[0] * (mat->val[5] * mat->val[10] - mat->val[9] * mat->val[6]) - 
        mat->val[0] * (mat->val[1] * mat->val[10] - mat->val[9] * mat->val[2]) +
        mat->val[4] * (mat->val[1] * mat->val[6]  - mat->val[5] * mat->val[2]);
}

C_UTILS_INLINE static inline void mat2x2_muls(
        const Mat2x2 *__restrict src, Mat2x2 *__restrict dst, f32 s) {
    ROW_128(dst, 0) = _mm_mul_ps(ROW_128(src, 0), _mm_set1_ps(s));
}

C_UTILS_INLINE static inline vec2f mat2x2_mulv(const Mat2x2 *__restrict src, vec2f v) {
    __m128 _tmp_0 = _mm_shuffle_ps(ROW_128(src, 0), ROW_128(src, 0), _MM_SHUFFLE(3, 1, 2, 0));
    __m128 _tmp_1 = _mm_set_ps(GET_VEC2_X(v), GET_VEC2_Y(v), GET_VEC2_X(v), GET_VEC2_Y(v));
    _tmp_1 = _mm_mul_ps(_tmp_0, _tmp_1);
    _tmp_1 = _mm_hadd_ps(_tmp_1, _tmp_1);

    return (vec2f) {
        .x = _mm_cvtss_f32(_tmp_1),
        .y = _mm_cvtss_f32(_mm_movehl_ps(_tmp_1, _tmp_1))
    };
}

C_UTILS_INLINE static inline void mat2x2_mulm(
        const Mat2x2 *__restrict a, const Mat2x2 *__restrict b, Mat2x2 *__restrict c) {
    __m128 _tmp_0 = _mm_shuffle_epi32(_mm_castps_si128(ROW_128(b, 0)), _MM_SHUFFLE(3, 0, 3, 0));
    __m128 _tmp_1 = _mm_shuffle_epi32(_mm_castps_si128(ROW_128(a, 0)), _MM_SHUFFLE(2, 3, 0, 1));
    __m128 _tmp_2 = _mm_shuffle_epi32(_mm_castps_si128(ROW_128(b, 0)), _MM_SHUFFLE(1, 2, 1, 2));

    __m128 _tmp_3 = _mm_mul_ps(ROW_128(a, 0), _mm_castsi128_ps(_tmp_0));
    __m128 _tmp_4 = _mm_fmadd_ps(_mm_castsi128_ps(_tmp_1), _mm_castsi128_ps(_tmp_2), _tmp_3);

    ROW_128(c, 0) = _tmp_4;
}

C_UTILS_INLINE static inline void mat2x2_transpose(
        const Mat2x2 *__restrict src, Mat2x2 *__restrict dst) {
    ROW_128(dst, 0) = _mm_shuffle_ps(ROW_128(src, 0), ROW_128(src, 0), _MM_SHUFFLE(3, 1, 2, 0));
}

C_UTILS_INLINE static inline void mat2x2_adjmul(
        const Mat2x2 *__restrict a, const Mat2x2 *__restrict b, Mat2x2 *__restrict c) {
    __m128 _tmp_0 = _mm_shuffle_epi32(_mm_castps_si128(ROW_128(a, 0)), _MM_SHUFFLE(0, 3, 0, 3));
    __m128 _tmp_1 = _mm_shuffle_epi32(_mm_castps_si128(ROW_128(a, 0)), _MM_SHUFFLE(2, 2, 1, 1));
    __m128 _tmp_2 = _mm_shuffle_epi32(_mm_castps_si128(ROW_128(b, 0)), _MM_SHUFFLE(1, 0, 3, 2));

    __m128 _tmp_3 = _mm_mul_ps(_mm_castsi128_ps(_tmp_1), _mm_castsi128_ps(_tmp_2));
    __m128 _tmp_4 = _mm_fmsub_ps(_mm_castsi128_ps(_tmp_0), ROW_128(b, 0), _tmp_3);

    ROW_128(c, 0) = _tmp_4;
}

C_UTILS_INLINE static inline void mat2x2_muladj(
        const Mat2x2 *__restrict a, const Mat2x2 *__restrict b, Mat2x2 *__restrict c) {
    __m128 _tmp_0 = _mm_shuffle_epi32(_mm_castps_si128(ROW_128(b, 0)), _MM_SHUFFLE(0, 3, 0, 3));
    __m128 _tmp_1 = _mm_shuffle_epi32(_mm_castps_si128(ROW_128(a, 0)), _MM_SHUFFLE(2, 3, 0, 1));
    __m128 _tmp_2 = _mm_shuffle_epi32(_mm_castps_si128(ROW_128(b, 0)), _MM_SHUFFLE(1, 2, 1, 2));

    __m128 _tmp_3 = _mm_mul_ps(_mm_castsi128_ps(_tmp_1), _mm_castsi128_ps( _tmp_2));
    __m128 _tmp_4 = _mm_fmsub_ps(ROW_128(a, 0), _mm_castsi128_ps(_tmp_0), _tmp_3);

    ROW_128(c, 0) = _tmp_4;
}

C_UTILS_INLINE static inline void mat4x4_muls(
        const Mat4x4 *__restrict src, Mat4x4 *__restrict dst, f32 s) {
#if defined(C_UTILS_AVX512F)
    __m512 _tmp_0 = _mm512_set1_ps(s);

    // this does not care about row-major / column-major ordering
    ROW_512(dst, 0) = _mm512_mul_ps(ROW_512(src, 0), _tmp_0);

#elif defined(C_UTILS_AVX)
    __m256 _tmp_0 = _mm256_set1_ps(s);

    // this does not care about row-major / column-major ordering
    ROW_256(dst, 0) = _mm256_mul_ps(ROW_256(src, 0), _tmp_0);
    ROW_256(dst, 1) = _mm256_mul_ps(ROW_256(src, 1), _tmp_0);

#else
    __m128 _tmp_0 = _mm_set1_ps(s);

    // this does not care about row-major / column-major ordering
    ROW_128(dst, 0) = _mm_mul_ps(ROW_128(src, 0), _tmp_0);
    ROW_128(dst, 1) = _mm_mul_ps(ROW_128(src, 1), _tmp_0);
    ROW_128(dst, 2) = _mm_mul_ps(ROW_128(src, 2), _tmp_0);
    ROW_128(dst, 3) = _mm_mul_ps(ROW_128(src, 3), _tmp_0);
#endif
}

C_UTILS_INLINE static inline void mat4x4_transpose(
        const Mat4x4 *__restrict src, Mat4x4 *__restrict dst) {
    __m128 _tmp_0 = _mm_shuffle_ps(ROW_128(src, 0), ROW_128(src, 1), 0x44);
    __m128 _tmp_1 = _mm_shuffle_ps(ROW_128(src, 0), ROW_128(src, 1), 0xEE);
    __m128 _tmp_2 = _mm_shuffle_ps(ROW_128(src, 2), ROW_128(src, 3), 0x44);
    __m128 _tmp_3 = _mm_shuffle_ps(ROW_128(src, 2), ROW_128(src, 3), 0xEE);

    ROW_128(dst, 0) = _mm_shuffle_ps(_tmp_0, _tmp_2, 0x88);
    ROW_128(dst, 1) = _mm_shuffle_ps(_tmp_0, _tmp_2, 0xDD);
    ROW_128(dst, 2) = _mm_shuffle_ps(_tmp_1, _tmp_3, 0x88);
    ROW_128(dst, 3) = _mm_shuffle_ps(_tmp_1, _tmp_3, 0xDD);
}

C_UTILS_INLINE static inline vec4f mat4x4_mulv(const Mat4x4 *__restrict mat, vec4f v) {
    __m128 _tmp_0 = _mm_mul_ps(ROW_128(mat, 0), _mm_set1_ps(VEC4_GET(v, 0)));
    __m128 _tmp_1 = _mm_fmadd_ps(ROW_128(mat, 1), _mm_set1_ps(VEC4_GET(v, 1)), _tmp_0);
    __m128 _tmp_2 = _mm_fmadd_ps(ROW_128(mat, 2), _mm_set1_ps(VEC4_GET(v, 2)), _tmp_1);
    __m128 _tmp_3 = _mm_fmadd_ps(ROW_128(mat, 3), _mm_set1_ps(VEC4_GET(v, 3)), _tmp_2);

    return (vec4f) { .vec = _tmp_3 };
}

#if defined(C_UTILS_AVX)
C_UTILS_INLINE static inline __m256 mat4x4_mulv2(const Mat4x4 *__restrict mat, __m256 pv) {
    __m256 _tmp_0 = _mm256_shuffle_epi32(_mm256_castps_si256(pv), _MM_SHUFFLE(0, 0, 0, 0));
    __m256 _tmp_1 = _mm256_shuffle_epi32(_mm256_castps_si256(pv), _MM_SHUFFLE(1, 1, 1, 1));
    __m256 _tmp_2 = _mm256_shuffle_epi32(_mm256_castps_si256(pv), _MM_SHUFFLE(2, 2, 2, 2));
    __m256 _tmp_3 = _mm256_shuffle_epi32(_mm256_castps_si256(pv), _MM_SHUFFLE(3, 3, 3, 3));

    __m256 _tmp_4 = _mm256_mul_ps(ROW_256(mat, 0), _mm256_castsi256_ps(_tmp_0));
    __m256 _tmp_5 = _mm256_fmadd_ps(ROW_256(mat, 0), _mm256_castsi256_ps(_tmp_1), _tmp_4);
    __m256 _tmp_6 = _mm256_fmadd_ps(ROW_256(mat, 1), _mm256_castsi256_ps(_tmp_2), _tmp_5);
    __m256 _tmp_7 = _mm256_fmadd_ps(ROW_256(mat, 1), _mm256_castsi256_ps(_tmp_3), _tmp_6);

    return _tmp_7;
}
#endif

#if defined(C_UTILS_AVX512F)
C_UTILS_INLINE static inline __m512 mat4x4_mulv4(const Mat4x4 *__restrict mat, __m512 pv) {
    __m512 _tmp_0 = _mm512_shuffle_epi32(_mm512_castps_si512(pv), _MM_SHUFFLE(0, 0, 0, 0));
    __m512 _tmp_1 = _mm512_shuffle_epi32(_mm512_castps_si512(pv), _MM_SHUFFLE(1, 1, 1, 1));
    __m512 _tmp_2 = _mm512_shuffle_epi32(_mm512_castps_si512(pv), _MM_SHUFFLE(2, 2, 2, 2));
    __m512 _tmp_3 = _mm512_shuffle_epi32(_mm512_castps_si512(pv), _MM_SHUFFLE(3, 3, 3, 3));

    __m512 _tmp_4 = _mm512_castps256_ps512(ROW_256(mat, 0));
    _tmp_4 = _mm512_inserti64x4(_tmp_4, _mm256_castps_si256(ROW_256(mat, 0)), 1);

    __m512 _tmp_5 = _mm512_castps256_ps512(ROW_256(mat, 1));
    _tmp_5 = _mm512_inserti64x4(_tmp_5, _mm256_castps_si256(ROW_256(mat, 1)), 1);

    __m512 _tmp_6 = _mm512_mul_ps(_mm512_castsi512_ps(_tmp_4), _mm512_castsi512_ps(_tmp_0));
    __m512 _tmp_7 = _mm512_fmadd_ps(_mm512_castsi512_ps(_tmp_4), _mm512_castsi512_ps(_tmp_1), _tmp_6);
    __m512 _tmp_8 = _mm512_fmadd_ps(_mm512_castsi512_ps(_tmp_5), _mm512_castsi512_ps(_tmp_2), _tmp_7);
    __m512 _tmp_9 = _mm512_fmadd_ps(_mm512_castsi512_ps(_tmp_5), _mm512_castsi512_ps(_tmp_3), _tmp_8);

    return _tmp_9;
}
#endif

C_UTILS_INLINE static inline void mat4x4_mulm(
        const Mat4x4 *__restrict a, const Mat4x4 *__restrict b, Mat4x4 *__restrict c) {

    // 4x4 matrices are stored in column-major order rather than row-major to speed up
    // 4x4 matrix vector and matrix multiplication
    // this implies that b itself doesn't need to be transposed to linearize the multiplication
#if defined(C_UTILS_AVX512F)
    ROW_512(c, 0) = mat4x4_mulv4(a, ROW_512(b, 0));

#elif defined(C_UTILS_AVX)
    ROW_256(c, 0) = mat4x4_mulv2(a, ROW_256(b, 0));
    ROW_256(c, 1) = mat4x4_mulv2(a, ROW_256(b, 1));

#else
    ROW_128(c, 0) = mat4x4_mulv(a, *((vec4f *) b->val + 0)).vec;
    ROW_128(c, 1) = mat4x4_mulv(a, *((vec4f *) b->val + 1)).vec;
    ROW_128(c, 2) = mat4x4_mulv(a, *((vec4f *) b->val + 2)).vec;
    ROW_128(c, 3) = mat4x4_mulv(a, *((vec4f *) b->val + 3)).vec;

#endif
}

C_UTILS_INLINE static inline Mat4x4 mat4x4_from_cvec(const vec4f v[4]) {
    Mat4x4 mat; 
    
    // 4x4 matrices are stored in column-major order rather than row-major to speed up
    // 4x4 matrix vector and matrix multiplication
    memcpy(mat.val, v, sizeof(vec4f) * 4);
    return mat;
}

C_UTILS_INLINE static inline Mat4x4 mat4x4_from_rvec(const vec3f v[4]) {

    // transpose needed to ensure column-major order
    Mat4x4 mat = mat4x4_from_cvec(v);
    Mat4x4 dst;

    mat4x4_transpose(&mat, &dst);
    return dst;
}

C_UTILS_INLINE static inline Mat2x2 mat4x4_to_mat2x2(const Mat4x4 *__restrict src) {
    Mat2x2 mat = { .row = _mm_movelh_ps(ROW_128(src, 0), ROW_128(src, 1)) };
    Mat2x2 dst;

    // 2x2 matrices are stored in row-major order in contrast to
    // the column-major representation of 4x4 matrices
    mat2x2_transpose(&mat, &dst);
    return dst;
}

C_UTILS_INLINE static inline Mat4x4 mat2x2_to_mat4x4(const Mat2x2 *__restrict src) {

    // assumes row-major ordering in 2x2 matrices
    Mat2x2 mat;
    mat2x2_transpose(src, &mat);

    Mat4x4 dst = mat4x4_zeroed();
    ROW_128(&dst, 0) = _mm_move_sd(_mm_setzero_ps(), ROW_128(src, 0));
    ROW_128(&dst, 1) = _mm_movehl_ps(_mm_setzero_ps(), ROW_128(src, 0));

    return dst;
}

/**
 * @brief Computes the inverse of an arbitrary 2x2 matrix in row-major order
 * @param src The source transformation matrix.
 * @param dst The destination transformation matrix where we want to store the inverse of src.
 */
void mat2x2_inverse(const Mat2x2 *__restrict src, Mat2x2 *__restrict dst);

/**
 * @brief Computes the inverse of an arbitrary 4x4 matrix in column-major order
 * @param src The source transformation matrix.
 * @param dst The destination transformation matrix where we want to store the inverse of src.
 */
void mat4x4_inverse(const Mat4x4 *__restrict src, Mat4x4 *__restrict dst);

/**
 * @brief Computes the inverse of an affine transformation matrix in column-major order
 *        which does not apply any scaling.
 * @param src The source transformation matrix.
 * @param dst The destination transformation matrix where we want to store the inverse of src.
 */
void mat4x4_inverse_tns(const Mat4x4 *__restrict src, Mat4x4 *__restrict dst);

/**
 * @brief Computes the inverse of a general transformation matrix in column-major order
 * @param src The source transformation matrix.
 * @param dst The destination transformation matrix where we want to store the inverse of src.
 */
void mat4x4_inverse_t(const Mat4x4 *__restrict src, Mat4x4 *__restrict dst);





/**
 * @brief Definition of an arbitrary run-time sized m x n matrix.
 *        The storage of such allocation will always be aligned to MAT_ALIGNMENT.
 *        The storage will always be padded to a multiple of 
 */
typedef struct Mat_t {
    // convention for width, height naming scheme of matrices
    // n is the width and m the height of the matrix
    usize n, m;

    // compiler assumption for the pointed to location to be aligned to MAT_ALIGNMENT
    // the assumption will be applied through runtime overload of function using this type. 
    f32 *val;
} Mat;

#define MAT_GET(_mat, _x, _y) \
    ((_y) * RND_POW_2((_mat)->n, MAT_ALIGNMENT) + (_x))

void *linalg_malloc(usize align, usize n, usize m) ;

void mat_mulv(const Mat *__restrict mat, const Vec *__restrict v, Vec *__restrict w);

void mat_mulm(const Mat *__restrict a, const Mat *__restrict b, Mat *__restrict c);

C_GUARD_END()

#endif //C_UTILS
