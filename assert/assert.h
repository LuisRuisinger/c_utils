#ifndef C_UTILS_ASSERT_H
#define C_UTILS_ASSERT_H

#include "../build.h"
#include "../defines.h"
#include "../fmt/fmt.h"

#ifdef C_UTILS_CPP11_LEAST 
    #define STATIC_ASSERT(...) \
	    static_assert(...)
#else
    #define STATIC_ASSERT(_c) \
	    do { (void) sizeof(u8 [1 - 2 * !(_c)]); } while(0)
#endif

#define __FMT_ASSERT(_m, ...) \
    do { fmt_log(LEVEL_ERROR, __FILENAME__, __LINE__, __func__, _m, ##__VA_ARGS__);  } while(0)

#define ASSERT(...) VFUNC(ASSERT, __VA_ARGS__)

#define ASSERT_1(_c) \
    do { if(_c) { __FMT_ASSERT(""); TRAP(); } } while (0)
#define ASSERT_2(_c, _m_0) \
    do { if(_c) { __FMT_ASSERT(_m_0); TRAP(); } } while (0)
#define ASSERT_3(_c, _m_0, _m_1) \
    do { if(_c) { __FMT_ASSERT(_m_0, _m_1); TRAP(); } } while (0)
#define ASSERT_4(_c, _m_0, _m_1, _m_2) \
    do { if(_c) { __FMT_ASSERT(_m_0, _m_1, _m_2); TRAP(); } } while (0)
#define ASSERT_5(_c, _m_0, _m_1, _m_2, _m_3) \
    do { if(_c) { __FMT_ASSERT(_m_0, _m_1, _m_2, _m_3); TRAP(); } } while (0)
#define ASSERT_6(_c, _m_0, _m_1, _m_2, _m_3, _m_4) \
    do { if(_c) { __FMT_ASSERT(_m_0, _m_1, _m_2, _m_3, _m_4); TRAP(); } } while (0)
#define ASSERT_7(_c, _m_0, _m_1, _m_2, _m_3, _m_4, _m_5) \
    do { if(_c) { __FMT_ASSERT(_m_0, _m_1, _m_2, _m_3, _m_4, _m_5); TRAP(); } } while (0)
#define ASSERT_8(_c, _m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6) \
    do { if(_c) { __FMT_ASSERT(_m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6); TRAP(); } } while (0)
#define ASSERT_9(_c, _m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6, _m_7) \
    do { if(_c) { __FMT_ASSERT(_m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6, _m_7); TRAP(); } } while (0)
#define ASSERT_10(_c, _m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6, _m_7, _m_8) \
    do { if(_c) { __FMT_ASSERT(_m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6, _m_7, _m_8); TRAP(); } } while (0)
#define ASSERT_11(_c, _m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6, _m_7, _m_8, _m_9) \
    do { if(_c) { __FMT_ASSERT(_m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6, _m_7, _m_8, _m_9); TRAP(); } } while (0)
#define ASSERT_12(_c, _m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6, _m_7, _m_8, _m_9, _m_10) \
    do { if(_c) { __FMT_ASSERT(_m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6, _m_7, _m_8, _m_9, _m_10); TRAP(); } } while (0)
#define ASSERT_13(_c, _m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6, _m_7, _m_8, _m_9, _m_10, _m_11) \
    do { if(_c) { __FMT_ASSERT(_m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6, _m_7, _m_8, _m_9, _m_10, _m_11); TRAP(); } } while (0)
#define ASSERT_14(_c, _m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6, _m_7, _m_8, _m_9, _m_10, _m_11, _m_12) \
    do { if(_c) { __FMT_ASSERT(_m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6, _m_7, _m_8, _m_9, _m_10, _m_11, _m_12); TRAP(); } } while (0)
#define ASSERT_15(_c, _m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6, _m_7, _m_8, _m_9, _m_10, _m_11, _m_12, _m_13) \
    do { if(_c) { __FMT_ASSERT(_m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6, _m_7, _m_8, _m_9, _m_10, _m_11, _m_12, _m_13); TRAP(); } } while (0)
#define ASSERT_16(_c, _m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6, _m_7, _m_8, _m_9, _m_10, _m_11, _m_12, _m_13, _m_14) \
    do { if(_c) { __FMT_ASSERT(_m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6, _m_7, _m_8, _m_9, _m_10, _m_11, _m_12, _m_13, _m_14); TRAP(); } } while (0)
#define ASSERT_17(_c, _m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6, _m_7, _m_8, _m_9, _m_10, _m_11, _m_12, _m_13, _m_14, _m_15) \
    do { if(_c) { __FMT_ASSERT(_m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6, _m_7, _m_8, _m_9, _m_10, _m_11, _m_12, _m_13, _m_14, _m_15); TRAP(); } } while (0)
#define ASSERT_18(_c, _m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6, _m_7, _m_8, _m_9, _m_10, _m_11, _m_12, _m_13, _m_14, _m_15, _m_16) \
    do { if(_c) { __FMT_ASSERT(_m_0, _m_1, _m_2, _m_3, _m_4, _m_5, _m_6, _m_7, _m_8, _m_9, _m_10, _m_11, _m_12, _m_13, _m_14, _m_15, _m_16); TRAP(); } } while (0)

#endif // C_UTILS_ASSERT_H