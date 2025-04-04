//
// Created by Luis Ruisinger on 17.12.24.
//

#ifndef C_UTILS_DICT_H
#define C_UTILS_DICT_H

#include "array/dynamic_array.h"
#include "defines.h"

C_GUARD_BEGINN()

#define ENTRY(_key, _value) \
    struct { _key key; _value value; }

#define ENTRY_FUNC(_key, _value) \
    struct { _key key; _value; }

#define MAP(_key, _value)                                                                           \
    struct {                                                                                        \
        ARRAY(ENTRY(_key, _value)) entries;                                                         \
        i32 (*cmp) (const void *a, const void *b);                                                  \
    }

#define MAP_FUNC(_key, _value)                                                                      \
    struct {                                                                                        \
        ARRAY(ENTRY_FUNC(_key, _value)) entries;                                                    \
        i32 (*cmp) (const void *a, const void *b);                                                  \
    }

#define MAP_TYPEOF(_map) \
    ARRAY_TYPEOF((_map).entries)

#define MAP_INIT_1(_map) \
    do { ARRAY_INIT((_map).entries, 16); } while (0)

#define MAP_INIT_2(_map, _arr)  \
    do {                    \
        MAP_INIT_1(_map);    \
        ARRAY_APPEND((_map).entries, _arr); \
    } while (0)

#define GET_MAP_INIT(_1, _2, NAME, ...) NAME

#define MAP_INIT(...) \
    GET_MAP_INIT(__VA_ARGS__, MAP_INIT_2, MAP_INIT_1)(__VA_ARGS__)

#define MAP_CONTAINS(_map, _key, _res)                                                              \
    do {                                                                                            \
        __typeof(_key) *__key = &_key;                                                              \
                                                                                                    \
        ARRAY_TYPEOF((_map).entries) *entry = NULL;                                                 \
        ARRAY_FOREACH((_map).entries, entry) {                                                      \
            if ((_map).cmp_fun(__key, &(entry)->key) == 0) {                                        \
                *(_res) = true;                                                                     \
                break;                                                                              \
            }                                                                                       \
        }                                                                                           \
    }                                                                                               \
    while (0)

#define MAP_INSERT(_map, _e, _res)                                                                  \
    do {                                                                                            \
        MAP_CONTAINS(_map, (_e).key, _res);                                                         \
                                                                                                    \
        if (!*(_res)) {                                                                             \
            ARRAY_PUSH_BACK((_map).entries, _e);                                                    \
            *(_res) = true;                                                                         \
        }                                                                                           \
    }                                                                                               \
    while (0)

#define MAP_GET(_map, _key, _res)                                                                   \
    do {                                                                                            \
        __typeof(_key) *__key = &_key;                                                              \
                                                                                                    \
        ARRAY_TYPEOF((_map).entries) *_entry = NULL;                                                \
        ARRAY_FOREACH((_map).entries, _entry) {                                                     \
            if ((_map).cmp(__key, &(_entry)->key) == 0) {                                           \
                *((ARRAY_TYPEOF((_map).entries) *) _res) = *_entry;                                 \
                break;                                                                              \
            }                                                                                       \
        }                                                                                           \
    }                                                                                               \
    while (0)

#define MAP_SORT_1(_map) \
    do { ARRAY_SORT((_map).entries, (_map).cmp); } while (0)

#define MAP_SORT_2(_map, _cmp) \
    do { ARRAY_SORT((_map).entries, (_cmp)); } while (0)

#define GET_MAP_SORT(_1, _2, NAME, ...) NAME

#define MAP_SORT(...) \
    GET_MAP_SORT(__VA_ARGS__, MAP_SORT_2, MAP_SORT_1)(__VA_ARGS__)

C_GUARD_END()

#endif //C_UTILS