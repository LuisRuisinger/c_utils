#include "matrix.h"
#include "assert/assert.h"

void mat2x2_inverse(const Mat2x2 *__restrict src, Mat2x2 *__restrict dst) {
    f32 s = 1.0F / mat2x2_det(src);

    // this is strictly speaking wrong because we store 2 rows inside of one __m128
    // compared to a 4x4 matrix where one __m128 is exactly one row
    ROW_128(dst, 0) = _mm_shuffle_ps(ROW_128(src, 0), ROW_128(src, 0), _MM_SHUFFLE(3, 1, 2, 0));
    ROW_128(dst, 0) = _mm_mul_ps(ROW_128(dst, 0), _mm_set_ps(1.0F, -1.0F, -1.0F, 1.0F));
    ROW_128(dst, 0) = _mm_mul_ps(ROW_128(dst, 0), _mm_set1_ps(s));
}

void mat4x4_inverse(const Mat4x4 *__restrict src, Mat4x4 *__restrict dst) {

    // 4x4 matrices are stored in row-major order rather than column-major to speed up
    // 4x4 matrix vector and matrix multiplication
    // the current inverse implementation assumes a column-major 4x4 matrix
    // this fact itself is subject to change but currently of less priority
    // because inverse itself is a costly function which isn'mat really used often
    Mat4x4 mat;
    mat4x4_transpose(src, &mat);

    Mat2x2 a = { .row = _mm_movelh_ps(ROW_128(&mat, 0), ROW_128(&mat, 1)) };
    Mat2x2 b = { .row = _mm_movehl_ps(ROW_128(&mat, 1), ROW_128(&mat, 0)) };
    Mat2x2 c = { .row = _mm_movelh_ps(ROW_128(&mat, 2), ROW_128(&mat, 3)) };
    Mat2x2 d = { .row = _mm_movehl_ps(ROW_128(&mat, 3), ROW_128(&mat, 2)) };

    // det
    __m128 _tmp_0 = _mm_shuffle_ps(ROW_128(&mat, 0), ROW_128(&mat, 2), _MM_SHUFFLE(2, 0, 2, 0));
    __m128 _tmp_1 = _mm_shuffle_ps(ROW_128(&mat, 1), ROW_128(&mat, 3), _MM_SHUFFLE(3, 1, 3, 1));
    __m128 _tmp_2 = _mm_shuffle_ps(ROW_128(&mat, 0), ROW_128(&mat, 2), _MM_SHUFFLE(3, 1, 3, 1));
    __m128 _tmp_3 = _mm_shuffle_ps(ROW_128(&mat, 1), ROW_128(&mat, 3), _MM_SHUFFLE(2, 0, 2, 0));
    __m128 _tmp_4 = _mm_fmsub_ps(_tmp_0, _tmp_1, _mm_mul_ps(_tmp_2, _tmp_3));

    __m128 _det_a = _mm_shuffle_ps(_tmp_4, _tmp_4, _MM_SHUFFLE(0, 0, 0, 0));
    __m128 _det_b = _mm_shuffle_ps(_tmp_4, _tmp_4, _MM_SHUFFLE(1, 1, 1, 1));
    __m128 _det_c = _mm_shuffle_ps(_tmp_4, _tmp_4, _MM_SHUFFLE(2, 2, 2, 2));
    __m128 _det_d = _mm_shuffle_ps(_tmp_4, _tmp_4, _MM_SHUFFLE(3, 3, 3, 3));

    Mat2x2 d_c, a_b;
    mat2x2_adjmul(&d, &c, &d_c);
    mat2x2_adjmul(&a, &b, &a_b);

    Mat2x2 tmp_0, tmp_1;
    mat2x2_mulm(&b, &d_c, &tmp_0);
    mat2x2_mulm(&c, &a_b, &tmp_1);

    __m128 _x = _mm_fmsub_ps(_det_d, ROW_128(&a, 0), ROW_128(&tmp_0, 0));
    __m128 _w = _mm_fmsub_ps(_det_a, ROW_128(&d, 0), ROW_128(&tmp_1, 0));

    Mat2x2 tmp_2, tmp_3;
    mat2x2_muladj(&d, &a_b, &tmp_2);
    mat2x2_muladj(&a, &d_c, &tmp_3);

    __m128 _y = _mm_fmsub_ps(_det_b, ROW_128(&c, 0), ROW_128(&tmp_2, 0));
    __m128 _z = _mm_fmsub_ps(_det_c, ROW_128(&b, 0), ROW_128(&tmp_3, 0));

    __m128 _tr = _mm_shuffle_epi32(_mm_castps_si128(ROW_128(&d_c, 0)), _MM_SHUFFLE(3, 1, 2, 0));
    _tr = _mm_mul_ps(ROW_128(&a_b, 0), _mm_castsi128_ps(_tr));
    _tr = _mm_hadd_ps(_tr, _tr);
    _tr = _mm_hadd_ps(_tr, _tr);

    // matrix det
    __m128 _det = _mm_fmadd_ps(_det_a, _det_d, _mm_mul_ps(_det_b, _det_c));
    _det = _mm_sub_ps(_det, _tr);
    _det = _mm_div_ps(_mm_setr_ps(1.0F, -1.0F, -1.0F, 1.0F), _det);

    _x = _mm_mul_ps(_x, _det);
    _y = _mm_mul_ps(_y, _det);
    _z = _mm_mul_ps(_z, _det);
    _w = _mm_mul_ps(_w, _det);

    ROW_128(dst, 0) = _mm_shuffle_ps(_x, _y, _MM_SHUFFLE(1, 3, 1, 3));
    ROW_128(dst, 1) = _mm_shuffle_ps(_x, _y, _MM_SHUFFLE(0, 2, 0, 2));
    ROW_128(dst, 2) = _mm_shuffle_ps(_z, _w, _MM_SHUFFLE(1, 3, 1, 3));
    ROW_128(dst, 3) = _mm_shuffle_ps(_z, _w, _MM_SHUFFLE(0, 2, 0, 2));
}

void mat4x4_inverse_tns(const Mat4x4 *__restrict src, Mat4x4 *__restrict dst) {

    // 4x4 matrices are stored in row-major order rather than column-major to speed up
    // 4x4 matrix vector and matrix multiplication
    // the current inverse implementation assumes a column-major 4x4 matrix
    // this fact itself is subject to change but currently of less priority
    // because inverse itself is a costly function which isn'mat really used often
    Mat4x4 mat;
    mat4x4_transpose(src, &mat);

    Mat2x2 a = { .row = _mm_movelh_ps(ROW_128(&mat, 0), ROW_128(&mat, 1)) };
    Mat2x2 b = { .row = _mm_movehl_ps(ROW_128(&mat, 1), ROW_128(&mat, 0)) };

    ROW_128(dst, 0) = _mm_shuffle_ps(ROW_128(&a, 0), ROW_128(&mat, 2), _MM_SHUFFLE(3, 0, 2, 0));
    ROW_128(dst, 1) = _mm_shuffle_ps(ROW_128(&a, 0), ROW_128(&mat, 2), _MM_SHUFFLE(3, 1, 3, 1));
    ROW_128(dst, 2) = _mm_shuffle_ps(ROW_128(&b, 0), ROW_128(&mat, 2), _MM_SHUFFLE(3, 2, 2, 0));

    __m128 _tmp_0 = _mm_shuffle_epi32(_mm_castps_si128(ROW_128(&mat, 0)), _MM_SHUFFLE(0, 0, 0, 0));
    __m128 _tmp_1 = _mm_shuffle_epi32(_mm_castps_si128(ROW_128(&mat, 0)), _MM_SHUFFLE(1, 1, 1, 1));
    __m128 _tmp_2 = _mm_shuffle_epi32(_mm_castps_si128(ROW_128(&mat, 0)), _MM_SHUFFLE(2, 2, 2, 2));

    _tmp_0 = _mm_mul_ps(ROW_128(dst, 0), _mm_castsi128_ps(_tmp_0));
    _tmp_1 = _mm_mul_ps(ROW_128(dst, 1), _mm_castsi128_ps(_tmp_1));
    _tmp_2 = _mm_mul_ps(ROW_128(dst, 2), _mm_castsi128_ps(_tmp_2));

    __m128 _tmp_3 = _mm_mul_ps(ROW_128(dst, 0), _mm_castsi128_ps(_tmp_0));
    __m128 _tmp_4 = _mm_fmadd_ps(ROW_128(dst, 1), _mm_castsi128_ps(_tmp_1), _tmp_3);
    __m128 _tmp_5 = _mm_fmadd_ps(ROW_128(dst, 2), _mm_castsi128_ps(_tmp_2), _tmp_4);

    ROW_128(dst, 3) = _mm_sub_ps(_mm_setr_ps(0.0F, 0.0F, 0.0F, 1.0F), _tmp_5);
}

void mat4x4_inverse_t(const Mat4x4 *__restrict src, Mat4x4 *__restrict dst) {

    // 4x4 matrices are stored in row-major order rather than column-major to speed up
    // 4x4 matrix vector and matrix multiplication
    // the current inverse implementation assumes a column-major 4x4 matrix
    // this fact itself is subject to change but currently of less priority
    // because inverse itself is a costly function which isn'mat really used often
    Mat4x4 mat;
    mat4x4_transpose(src, &mat);

    /* TODO */
    // maybe 2 2x2 transpose are faster; however I think this needs more complex extraction logic
    // and is therefore less feasible than the current single 4x4 transpose approach
    Mat2x2 a = { .row = _mm_movelh_ps(ROW_128(&mat, 0), ROW_128(&mat, 1)) };
    Mat2x2 b = { .row = _mm_movehl_ps(ROW_128(&mat, 1), ROW_128(&mat, 0)) };

    ROW_128(dst, 0) = _mm_shuffle_ps(ROW_128(&a, 0), ROW_128(&mat, 2), _MM_SHUFFLE(3, 0, 2, 0));
    ROW_128(dst, 1) = _mm_shuffle_ps(ROW_128(&a, 0), ROW_128(&mat, 2), _MM_SHUFFLE(3, 1, 3, 1));
    ROW_128(dst, 2) = _mm_shuffle_ps(ROW_128(&b, 0), ROW_128(&mat, 2), _MM_SHUFFLE(3, 2, 2, 0));

    __m128 _tmp_0 = _mm_mul_ps(ROW_128(dst, 0), ROW_128(dst, 0));
    _tmp_0 = _mm_fmadd_ps(ROW_128(dst, 1), ROW_128(dst, 1), _tmp_0);
    _tmp_0 = _mm_fmadd_ps(ROW_128(dst, 2), ROW_128(dst, 2), _tmp_0);

    __m128 _tmp_1 = _mm_set1_ps(1.0F);
    __m128 _tmp_2 = _mm_div_ps(_tmp_1, _tmp_0);
    _tmp_2 = _mm_blendv_ps(_tmp_2, _tmp_1, _mm_cmplt_ps(_tmp_0, _mm_set1_ps(1.0E-8F)));

    ROW_128(dst, 0) = _mm_mul_ps(ROW_128(dst, 0), _tmp_2);
    ROW_128(dst, 1) = _mm_mul_ps(ROW_128(dst, 1), _tmp_2);
    ROW_128(dst, 2) = _mm_mul_ps(ROW_128(dst, 2), _tmp_2);

    __m128 _tmp_3 = _mm_shuffle_epi32(_mm_castps_si128(ROW_128(&mat, 0)), _MM_SHUFFLE(0, 0, 0, 0));
    __m128 _tmp_4 = _mm_shuffle_epi32(_mm_castps_si128(ROW_128(&mat, 0)), _MM_SHUFFLE(1, 1, 1, 1));
    __m128 _tmp_5 = _mm_shuffle_epi32(_mm_castps_si128(ROW_128(&mat, 0)), _MM_SHUFFLE(2, 2, 2, 2));

    __m128 _tmp_6 = _mm_mul_ps(ROW_128(dst, 0), _mm_castsi128_ps(_tmp_3));
    _tmp_6 = _mm_fmadd_ps(ROW_128(dst, 1), _mm_castsi128_ps(_tmp_4), _tmp_6);
    _tmp_6 = _mm_fmadd_ps(ROW_128(dst, 2), _mm_castsi128_ps(_tmp_5), _tmp_6);

    ROW_128(dst, 3) = _mm_sub_ps(_mm_setr_ps(0.0F, 0.0F, 0.0F, 1.0F), _tmp_6);
}

void mat_mulv(const Mat *__restrict m, const Vec *__restrict v, Vec *__restrict w) {
    if (m->n != v->m || m->m != w->m) {
        // ASSERT(false);
        exit(EXIT_FAILURE);
    }
    
    f32 *m_val = C_UTILS_ASSUME_ALIGNED(m->val, sizeof(MAT_ALIGNMENT)); 
    f32 *v_val = C_UTILS_ASSUME_ALIGNED(v->val, sizeof(MAT_ALIGNMENT)); 
    f32 *w_val = C_UTILS_ASSUME_ALIGNED(w->val, sizeof(MAT_ALIGNMENT)); 
}

void *linalg_malloc(usize align, usize n, usize m) {
    n = RND_POW_2(n, MAT_ALIGNMENT);
    m = RND_POW_2(m, MAT_ALIGNMENT);

    return aligned_alloc(align, n * m * sizeof(f32));
}

usize l1_cache_size;
usize l2_cache_size;
usize l3_cache_size;

#define L1_CACHE_SIZE l1_cache_size
#define L2_CACHE_SIZE l2_cache_size
#define L3_CACHE_SiZE l3_cache_size

C_UTILS_ON_STARTUP void linalg_setup(void) {
    LOG("hello");

    /* TODO */
    // call cpuid and fetch information about the cache layers 
    // to build divide-and-conquer matrices based on cache information
}

#define DEFAULT_MAT_MAT_MUL

void mat_mulm(const Mat *__restrict a, const Mat *__restrict b, Mat *__restrict c) {
    if (a->n != b->m || a->m != c->m || b->n != c->n) {
        // ASSERT(false);
        LOG("ERROR");
        exit(EXIT_FAILURE);
    }
    
    f32 *a_val = C_UTILS_ASSUME_ALIGNED(a->val, sizeof(MAT_ALIGNMENT)); 
    f32 *b_val = C_UTILS_ASSUME_ALIGNED(b->val, sizeof(MAT_ALIGNMENT)); 
    f32 *c_val = C_UTILS_ASSUME_ALIGNED(c->val, sizeof(MAT_ALIGNMENT)); 

#ifdef DEFAULT_MAT_MAT_MUL
    for (usize i = 0; i < a->m; ++i) {
        for (usize j = 0; j < b->n; ++j) {
            for (usize k = 0; k < a->n; ++k) {
                c_val[i * c->n + j] += a_val[i * a->n + k] * b_val[k * b->n + j];
            }
        }
    }
#endif
}



