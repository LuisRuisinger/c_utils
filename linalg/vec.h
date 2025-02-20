#ifndef C_UTILS_VEC_H
#define C_UTILS_VEC_H

#include "../defines.h"

typedef struct Vec_t {
    
    // continous storage of vector values
    // aligned to the biggest available SIMD register
    // alignment resolved at runtime
    f32 *val;
    
    // height of the row vector
    usize m;
} Vec;

#endif //C_UTILS_VEC_H