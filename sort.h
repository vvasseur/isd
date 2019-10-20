/*
   Copyright (c) 2010-2017 Christopher Swenson.
   Copyright (c) 2012 Vojtech Fried.
   Copyright (c) 2012 Google Inc. All Rights Reserved.
   Copyright (c) 2019 Valentin Vasseur

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to
   deal in the Software without restriction, including without limitation the
   rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
   sell copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
   IN THE SOFTWARE
*/

/*
 * This code was taken from: https://github.com/swenson/sort
 *
 * Everything but the functions needed for the quick sort have been removed.
 * The quick sort function has been modified to sort two arrays at the same
 * time, the first array 'idx' is supposed to contain the index of the elements
 * being sorted.
 */

#include <stdint.h>
#include <stdlib.h>

#ifndef SORT_NAME
#error "Must declare SORT_NAME"
#endif

#ifndef SORT_TYPE
#error "Must declare SORT_TYPE"
#endif

#ifndef SORT_CMP
#define SORT_CMP(x, y) ((x) < (y) ? -1 : ((x) == (y) ? 0 : 1))
#endif

#ifndef IDX_SWAP
#define IDX_SWAP(x, y)            \
  {                               \
    size_t _sort_swap_temp = (x); \
    (x) = (y);                    \
    (y) = _sort_swap_temp;        \
  }
#endif

#ifndef SORT_SWAP
#define SORT_SWAP(x, y)              \
  {                                  \
    SORT_TYPE _sort_swap_temp = (x); \
    (x) = (y);                       \
    (y) = _sort_swap_temp;           \
  }
#endif

#define SORT_CONCAT(x, y) x##_##y
#define SORT_MAKE_STR1(x, y) SORT_CONCAT(x, y)
#define SORT_MAKE_STR(x) SORT_MAKE_STR1(SORT_NAME, x)

#ifndef SMALL_SORT_BND
#define SMALL_SORT_BND 16
#endif
#ifndef SMALL_SORT_PAIR
#define SMALL_SORT_PAIR BITONIC_SORT_PAIR
#endif

#define BITONIC_SORT_PAIR SORT_MAKE_STR(bitonic_sort_pair)
#define QUICK_SORT_PAIR SORT_MAKE_STR(quick_sort_pair)
#define QUICK_SORT_PARTITION_PAIR SORT_MAKE_STR(quick_sort_partition_pair)
#define QUICK_SORT_RECURSIVE_PAIR SORT_MAKE_STR(quick_sort_recursive_pair)
#define BINARY_INSERTION_SORT_START_PAIR \
  SORT_MAKE_STR(binary_insertion_sort_start_pair)
#define BINARY_INSERTION_SORT_PAIR SORT_MAKE_STR(binary_insertion_sort_pair)

#ifndef SORT_CSWAP_PAIR
#define SORT_CSWAP_PAIR(x, y, idxx, idxy) \
  {                                       \
    if (SORT_CMP((x), (y)) > 0) {         \
      SORT_SWAP((x), (y));                \
      IDX_SWAP((idxx), (idxy));           \
    }                                     \
  }
#endif

void BINARY_INSERTION_SORT_PAIR(size_t *idx, SORT_TYPE *dst, const size_t size);
void QUICK_SORT_PAIR(size_t *idx, SORT_TYPE *dst, const size_t size);
void BITONIC_SORT_PAIR(size_t *idx, SORT_TYPE *dst, const size_t size);

/* The full implementation of a bitonic sort is not here. Since we only want to
   use sorting networks for small length lists we create optimal sorting
   networks for lists of length <= 16 and call out to BINARY_INSERTION_SORT for
   anything larger than 16. Optimal sorting networks for small length lists.
   Taken from https://pages.ripco.net/~jgamble/nw.html */
#define BITONIC_SORT_PAIR_2_PAIR SORT_MAKE_STR(bitonic_sort_2_pair)
static __inline void BITONIC_SORT_PAIR_2(size_t *idx, SORT_TYPE *dst) {
  SORT_CSWAP_PAIR(dst[0], dst[1], idx[0], idx[1]);
}

#define BITONIC_SORT_PAIR_3_PAIR SORT_MAKE_STR(bitonic_sort_3_pair)
static __inline void BITONIC_SORT_PAIR_3(size_t *idx, SORT_TYPE *dst) {
  SORT_CSWAP_PAIR(dst[1], dst[2], idx[1], idx[2]);
  SORT_CSWAP_PAIR(dst[0], dst[2], idx[0], idx[2]);
  SORT_CSWAP_PAIR(dst[0], dst[1], idx[0], idx[1]);
}

#define BITONIC_SORT_PAIR_4_PAIR SORT_MAKE_STR(bitonic_sort_4_pair)
static __inline void BITONIC_SORT_PAIR_4(size_t *idx, SORT_TYPE *dst) {
  SORT_CSWAP_PAIR(dst[0], dst[1], idx[0], idx[1]);
  SORT_CSWAP_PAIR(dst[2], dst[3], idx[2], idx[3]);
  SORT_CSWAP_PAIR(dst[0], dst[2], idx[0], idx[2]);
  SORT_CSWAP_PAIR(dst[1], dst[3], idx[1], idx[3]);
  SORT_CSWAP_PAIR(dst[1], dst[2], idx[1], idx[2]);
}

#define BITONIC_SORT_PAIR_5_PAIR SORT_MAKE_STR(bitonic_sort_5_pair)
static __inline void BITONIC_SORT_PAIR_5(size_t *idx, SORT_TYPE *dst) {
  SORT_CSWAP_PAIR(dst[0], dst[1], idx[0], idx[1]);
  SORT_CSWAP_PAIR(dst[3], dst[4], idx[3], idx[4]);
  SORT_CSWAP_PAIR(dst[2], dst[4], idx[2], idx[4]);
  SORT_CSWAP_PAIR(dst[2], dst[3], idx[2], idx[3]);
  SORT_CSWAP_PAIR(dst[1], dst[4], idx[1], idx[4]);
  SORT_CSWAP_PAIR(dst[0], dst[3], idx[0], idx[3]);
  SORT_CSWAP_PAIR(dst[0], dst[2], idx[0], idx[2]);
  SORT_CSWAP_PAIR(dst[1], dst[3], idx[1], idx[3]);
  SORT_CSWAP_PAIR(dst[1], dst[2], idx[1], idx[2]);
}

#define BITONIC_SORT_PAIR_6_PAIR SORT_MAKE_STR(bitonic_sort_6_pair)
static __inline void BITONIC_SORT_PAIR_6(size_t *idx, SORT_TYPE *dst) {
  SORT_CSWAP_PAIR(dst[1], dst[2], idx[1], idx[2]);
  SORT_CSWAP_PAIR(dst[4], dst[5], idx[4], idx[5]);
  SORT_CSWAP_PAIR(dst[0], dst[2], idx[0], idx[2]);
  SORT_CSWAP_PAIR(dst[3], dst[5], idx[3], idx[5]);
  SORT_CSWAP_PAIR(dst[0], dst[1], idx[0], idx[1]);
  SORT_CSWAP_PAIR(dst[3], dst[4], idx[3], idx[4]);
  SORT_CSWAP_PAIR(dst[2], dst[5], idx[2], idx[5]);
  SORT_CSWAP_PAIR(dst[0], dst[3], idx[0], idx[3]);
  SORT_CSWAP_PAIR(dst[1], dst[4], idx[1], idx[4]);
  SORT_CSWAP_PAIR(dst[2], dst[4], idx[2], idx[4]);
  SORT_CSWAP_PAIR(dst[1], dst[3], idx[1], idx[3]);
  SORT_CSWAP_PAIR(dst[2], dst[3], idx[2], idx[3]);
}

#define BITONIC_SORT_PAIR_7_PAIR SORT_MAKE_STR(bitonic_sort_7_pair)
static __inline void BITONIC_SORT_PAIR_7(size_t *idx, SORT_TYPE *dst) {
  SORT_CSWAP_PAIR(dst[1], dst[2], idx[1], idx[2]);
  SORT_CSWAP_PAIR(dst[3], dst[4], idx[3], idx[4]);
  SORT_CSWAP_PAIR(dst[5], dst[6], idx[5], idx[6]);
  SORT_CSWAP_PAIR(dst[0], dst[2], idx[0], idx[2]);
  SORT_CSWAP_PAIR(dst[3], dst[5], idx[3], idx[5]);
  SORT_CSWAP_PAIR(dst[4], dst[6], idx[4], idx[6]);
  SORT_CSWAP_PAIR(dst[0], dst[1], idx[0], idx[1]);
  SORT_CSWAP_PAIR(dst[4], dst[5], idx[4], idx[5]);
  SORT_CSWAP_PAIR(dst[2], dst[6], idx[2], idx[6]);
  SORT_CSWAP_PAIR(dst[0], dst[4], idx[0], idx[4]);
  SORT_CSWAP_PAIR(dst[1], dst[5], idx[1], idx[5]);
  SORT_CSWAP_PAIR(dst[0], dst[3], idx[0], idx[3]);
  SORT_CSWAP_PAIR(dst[2], dst[5], idx[2], idx[5]);
  SORT_CSWAP_PAIR(dst[1], dst[3], idx[1], idx[3]);
  SORT_CSWAP_PAIR(dst[2], dst[4], idx[2], idx[4]);
  SORT_CSWAP_PAIR(dst[2], dst[3], idx[2], idx[3]);
}

#define BITONIC_SORT_PAIR_8_PAIR SORT_MAKE_STR(bitonic_sort_8_pair)
static __inline void BITONIC_SORT_PAIR_8(size_t *idx, SORT_TYPE *dst) {
  SORT_CSWAP_PAIR(dst[0], dst[1], idx[0], idx[1]);
  SORT_CSWAP_PAIR(dst[2], dst[3], idx[2], idx[3]);
  SORT_CSWAP_PAIR(dst[4], dst[5], idx[4], idx[5]);
  SORT_CSWAP_PAIR(dst[6], dst[7], idx[6], idx[7]);
  SORT_CSWAP_PAIR(dst[0], dst[2], idx[0], idx[2]);
  SORT_CSWAP_PAIR(dst[1], dst[3], idx[1], idx[3]);
  SORT_CSWAP_PAIR(dst[4], dst[6], idx[4], idx[6]);
  SORT_CSWAP_PAIR(dst[5], dst[7], idx[5], idx[7]);
  SORT_CSWAP_PAIR(dst[1], dst[2], idx[1], idx[2]);
  SORT_CSWAP_PAIR(dst[5], dst[6], idx[5], idx[6]);
  SORT_CSWAP_PAIR(dst[0], dst[4], idx[0], idx[4]);
  SORT_CSWAP_PAIR(dst[3], dst[7], idx[3], idx[7]);
  SORT_CSWAP_PAIR(dst[1], dst[5], idx[1], idx[5]);
  SORT_CSWAP_PAIR(dst[2], dst[6], idx[2], idx[6]);
  SORT_CSWAP_PAIR(dst[1], dst[4], idx[1], idx[4]);
  SORT_CSWAP_PAIR(dst[3], dst[6], idx[3], idx[6]);
  SORT_CSWAP_PAIR(dst[2], dst[4], idx[2], idx[4]);
  SORT_CSWAP_PAIR(dst[3], dst[5], idx[3], idx[5]);
  SORT_CSWAP_PAIR(dst[3], dst[4], idx[3], idx[4]);
}

#define BITONIC_SORT_PAIR_9_PAIR SORT_MAKE_STR(bitonic_sort_9_pair)
static __inline void BITONIC_SORT_PAIR_9(size_t *idx, SORT_TYPE *dst) {
  SORT_CSWAP_PAIR(dst[0], dst[1], idx[0], idx[1]);
  SORT_CSWAP_PAIR(dst[3], dst[4], idx[3], idx[4]);
  SORT_CSWAP_PAIR(dst[6], dst[7], idx[6], idx[7]);
  SORT_CSWAP_PAIR(dst[1], dst[2], idx[1], idx[2]);
  SORT_CSWAP_PAIR(dst[4], dst[5], idx[4], idx[5]);
  SORT_CSWAP_PAIR(dst[7], dst[8], idx[7], idx[8]);
  SORT_CSWAP_PAIR(dst[0], dst[1], idx[0], idx[1]);
  SORT_CSWAP_PAIR(dst[3], dst[4], idx[3], idx[4]);
  SORT_CSWAP_PAIR(dst[6], dst[7], idx[6], idx[7]);
  SORT_CSWAP_PAIR(dst[2], dst[5], idx[2], idx[5]);
  SORT_CSWAP_PAIR(dst[0], dst[3], idx[0], idx[3]);
  SORT_CSWAP_PAIR(dst[1], dst[4], idx[1], idx[4]);
  SORT_CSWAP_PAIR(dst[5], dst[8], idx[5], idx[8]);
  SORT_CSWAP_PAIR(dst[3], dst[6], idx[3], idx[6]);
  SORT_CSWAP_PAIR(dst[4], dst[7], idx[4], idx[7]);
  SORT_CSWAP_PAIR(dst[2], dst[5], idx[2], idx[5]);
  SORT_CSWAP_PAIR(dst[0], dst[3], idx[0], idx[3]);
  SORT_CSWAP_PAIR(dst[1], dst[4], idx[1], idx[4]);
  SORT_CSWAP_PAIR(dst[5], dst[7], idx[5], idx[7]);
  SORT_CSWAP_PAIR(dst[2], dst[6], idx[2], idx[6]);
  SORT_CSWAP_PAIR(dst[1], dst[3], idx[1], idx[3]);
  SORT_CSWAP_PAIR(dst[4], dst[6], idx[4], idx[6]);
  SORT_CSWAP_PAIR(dst[2], dst[4], idx[2], idx[4]);
  SORT_CSWAP_PAIR(dst[5], dst[6], idx[5], idx[6]);
  SORT_CSWAP_PAIR(dst[2], dst[3], idx[2], idx[3]);
}

#define BITONIC_SORT_PAIR_10_PAIR SORT_MAKE_STR(bitonic_sort_10_pair)
static __inline void BITONIC_SORT_PAIR_10(size_t *idx, SORT_TYPE *dst) {
  SORT_CSWAP_PAIR(dst[4], dst[9], idx[4], idx[9]);
  SORT_CSWAP_PAIR(dst[3], dst[8], idx[3], idx[8]);
  SORT_CSWAP_PAIR(dst[2], dst[7], idx[2], idx[7]);
  SORT_CSWAP_PAIR(dst[1], dst[6], idx[1], idx[6]);
  SORT_CSWAP_PAIR(dst[0], dst[5], idx[0], idx[5]);
  SORT_CSWAP_PAIR(dst[1], dst[4], idx[1], idx[4]);
  SORT_CSWAP_PAIR(dst[6], dst[9], idx[6], idx[9]);
  SORT_CSWAP_PAIR(dst[0], dst[3], idx[0], idx[3]);
  SORT_CSWAP_PAIR(dst[5], dst[8], idx[5], idx[8]);
  SORT_CSWAP_PAIR(dst[0], dst[2], idx[0], idx[2]);
  SORT_CSWAP_PAIR(dst[3], dst[6], idx[3], idx[6]);
  SORT_CSWAP_PAIR(dst[7], dst[9], idx[7], idx[9]);
  SORT_CSWAP_PAIR(dst[0], dst[1], idx[0], idx[1]);
  SORT_CSWAP_PAIR(dst[2], dst[4], idx[2], idx[4]);
  SORT_CSWAP_PAIR(dst[5], dst[7], idx[5], idx[7]);
  SORT_CSWAP_PAIR(dst[8], dst[9], idx[8], idx[9]);
  SORT_CSWAP_PAIR(dst[1], dst[2], idx[1], idx[2]);
  SORT_CSWAP_PAIR(dst[4], dst[6], idx[4], idx[6]);
  SORT_CSWAP_PAIR(dst[7], dst[8], idx[7], idx[8]);
  SORT_CSWAP_PAIR(dst[3], dst[5], idx[3], idx[5]);
  SORT_CSWAP_PAIR(dst[2], dst[5], idx[2], idx[5]);
  SORT_CSWAP_PAIR(dst[6], dst[8], idx[6], idx[8]);
  SORT_CSWAP_PAIR(dst[1], dst[3], idx[1], idx[3]);
  SORT_CSWAP_PAIR(dst[4], dst[7], idx[4], idx[7]);
  SORT_CSWAP_PAIR(dst[2], dst[3], idx[2], idx[3]);
  SORT_CSWAP_PAIR(dst[6], dst[7], idx[6], idx[7]);
  SORT_CSWAP_PAIR(dst[3], dst[4], idx[3], idx[4]);
  SORT_CSWAP_PAIR(dst[5], dst[6], idx[5], idx[6]);
  SORT_CSWAP_PAIR(dst[4], dst[5], idx[4], idx[5]);
}

#define BITONIC_SORT_PAIR_11_PAIR SORT_MAKE_STR(bitonic_sort_11_pair)
static __inline void BITONIC_SORT_PAIR_11(size_t *idx, SORT_TYPE *dst) {
  SORT_CSWAP_PAIR(dst[0], dst[1], idx[0], idx[1]);
  SORT_CSWAP_PAIR(dst[2], dst[3], idx[2], idx[3]);
  SORT_CSWAP_PAIR(dst[4], dst[5], idx[4], idx[5]);
  SORT_CSWAP_PAIR(dst[6], dst[7], idx[6], idx[7]);
  SORT_CSWAP_PAIR(dst[8], dst[9], idx[8], idx[9]);
  SORT_CSWAP_PAIR(dst[1], dst[3], idx[1], idx[3]);
  SORT_CSWAP_PAIR(dst[5], dst[7], idx[5], idx[7]);
  SORT_CSWAP_PAIR(dst[0], dst[2], idx[0], idx[2]);
  SORT_CSWAP_PAIR(dst[4], dst[6], idx[4], idx[6]);
  SORT_CSWAP_PAIR(dst[8], dst[10], idx[8], idx[10]);
  SORT_CSWAP_PAIR(dst[1], dst[2], idx[1], idx[2]);
  SORT_CSWAP_PAIR(dst[5], dst[6], idx[5], idx[6]);
  SORT_CSWAP_PAIR(dst[9], dst[10], idx[9], idx[10]);
  SORT_CSWAP_PAIR(dst[0], dst[4], idx[0], idx[4]);
  SORT_CSWAP_PAIR(dst[3], dst[7], idx[3], idx[7]);
  SORT_CSWAP_PAIR(dst[1], dst[5], idx[1], idx[5]);
  SORT_CSWAP_PAIR(dst[6], dst[10], idx[6], idx[10]);
  SORT_CSWAP_PAIR(dst[4], dst[8], idx[4], idx[8]);
  SORT_CSWAP_PAIR(dst[5], dst[9], idx[5], idx[9]);
  SORT_CSWAP_PAIR(dst[2], dst[6], idx[2], idx[6]);
  SORT_CSWAP_PAIR(dst[0], dst[4], idx[0], idx[4]);
  SORT_CSWAP_PAIR(dst[3], dst[8], idx[3], idx[8]);
  SORT_CSWAP_PAIR(dst[1], dst[5], idx[1], idx[5]);
  SORT_CSWAP_PAIR(dst[6], dst[10], idx[6], idx[10]);
  SORT_CSWAP_PAIR(dst[2], dst[3], idx[2], idx[3]);
  SORT_CSWAP_PAIR(dst[8], dst[9], idx[8], idx[9]);
  SORT_CSWAP_PAIR(dst[1], dst[4], idx[1], idx[4]);
  SORT_CSWAP_PAIR(dst[7], dst[10], idx[7], idx[10]);
  SORT_CSWAP_PAIR(dst[3], dst[5], idx[3], idx[5]);
  SORT_CSWAP_PAIR(dst[6], dst[8], idx[6], idx[8]);
  SORT_CSWAP_PAIR(dst[2], dst[4], idx[2], idx[4]);
  SORT_CSWAP_PAIR(dst[7], dst[9], idx[7], idx[9]);
  SORT_CSWAP_PAIR(dst[5], dst[6], idx[5], idx[6]);
  SORT_CSWAP_PAIR(dst[3], dst[4], idx[3], idx[4]);
  SORT_CSWAP_PAIR(dst[7], dst[8], idx[7], idx[8]);
}

#define BITONIC_SORT_PAIR_12_PAIR SORT_MAKE_STR(bitonic_sort_12_pair)
static __inline void BITONIC_SORT_PAIR_12(size_t *idx, SORT_TYPE *dst) {
  SORT_CSWAP_PAIR(dst[0], dst[1], idx[0], idx[1]);
  SORT_CSWAP_PAIR(dst[2], dst[3], idx[2], idx[3]);
  SORT_CSWAP_PAIR(dst[4], dst[5], idx[4], idx[5]);
  SORT_CSWAP_PAIR(dst[6], dst[7], idx[6], idx[7]);
  SORT_CSWAP_PAIR(dst[8], dst[9], idx[8], idx[9]);
  SORT_CSWAP_PAIR(dst[10], dst[11], idx[10], idx[11]);
  SORT_CSWAP_PAIR(dst[1], dst[3], idx[1], idx[3]);
  SORT_CSWAP_PAIR(dst[5], dst[7], idx[5], idx[7]);
  SORT_CSWAP_PAIR(dst[9], dst[11], idx[9], idx[11]);
  SORT_CSWAP_PAIR(dst[0], dst[2], idx[0], idx[2]);
  SORT_CSWAP_PAIR(dst[4], dst[6], idx[4], idx[6]);
  SORT_CSWAP_PAIR(dst[8], dst[10], idx[8], idx[10]);
  SORT_CSWAP_PAIR(dst[1], dst[2], idx[1], idx[2]);
  SORT_CSWAP_PAIR(dst[5], dst[6], idx[5], idx[6]);
  SORT_CSWAP_PAIR(dst[9], dst[10], idx[9], idx[10]);
  SORT_CSWAP_PAIR(dst[0], dst[4], idx[0], idx[4]);
  SORT_CSWAP_PAIR(dst[7], dst[11], idx[7], idx[11]);
  SORT_CSWAP_PAIR(dst[1], dst[5], idx[1], idx[5]);
  SORT_CSWAP_PAIR(dst[6], dst[10], idx[6], idx[10]);
  SORT_CSWAP_PAIR(dst[3], dst[7], idx[3], idx[7]);
  SORT_CSWAP_PAIR(dst[4], dst[8], idx[4], idx[8]);
  SORT_CSWAP_PAIR(dst[5], dst[9], idx[5], idx[9]);
  SORT_CSWAP_PAIR(dst[2], dst[6], idx[2], idx[6]);
  SORT_CSWAP_PAIR(dst[0], dst[4], idx[0], idx[4]);
  SORT_CSWAP_PAIR(dst[7], dst[11], idx[7], idx[11]);
  SORT_CSWAP_PAIR(dst[3], dst[8], idx[3], idx[8]);
  SORT_CSWAP_PAIR(dst[1], dst[5], idx[1], idx[5]);
  SORT_CSWAP_PAIR(dst[6], dst[10], idx[6], idx[10]);
  SORT_CSWAP_PAIR(dst[2], dst[3], idx[2], idx[3]);
  SORT_CSWAP_PAIR(dst[8], dst[9], idx[8], idx[9]);
  SORT_CSWAP_PAIR(dst[1], dst[4], idx[1], idx[4]);
  SORT_CSWAP_PAIR(dst[7], dst[10], idx[7], idx[10]);
  SORT_CSWAP_PAIR(dst[3], dst[5], idx[3], idx[5]);
  SORT_CSWAP_PAIR(dst[6], dst[8], idx[6], idx[8]);
  SORT_CSWAP_PAIR(dst[2], dst[4], idx[2], idx[4]);
  SORT_CSWAP_PAIR(dst[7], dst[9], idx[7], idx[9]);
  SORT_CSWAP_PAIR(dst[5], dst[6], idx[5], idx[6]);
  SORT_CSWAP_PAIR(dst[3], dst[4], idx[3], idx[4]);
  SORT_CSWAP_PAIR(dst[7], dst[8], idx[7], idx[8]);
}

#define BITONIC_SORT_PAIR_13_PAIR SORT_MAKE_STR(bitonic_sort_13_pair)
static __inline void BITONIC_SORT_PAIR_13(size_t *idx, SORT_TYPE *dst) {
  SORT_CSWAP_PAIR(dst[1], dst[7], idx[1], idx[7]);
  SORT_CSWAP_PAIR(dst[9], dst[11], idx[9], idx[11]);
  SORT_CSWAP_PAIR(dst[3], dst[4], idx[3], idx[4]);
  SORT_CSWAP_PAIR(dst[5], dst[8], idx[5], idx[8]);
  SORT_CSWAP_PAIR(dst[0], dst[12], idx[0], idx[12]);
  SORT_CSWAP_PAIR(dst[2], dst[6], idx[2], idx[6]);
  SORT_CSWAP_PAIR(dst[0], dst[1], idx[0], idx[1]);
  SORT_CSWAP_PAIR(dst[2], dst[3], idx[2], idx[3]);
  SORT_CSWAP_PAIR(dst[4], dst[6], idx[4], idx[6]);
  SORT_CSWAP_PAIR(dst[8], dst[11], idx[8], idx[11]);
  SORT_CSWAP_PAIR(dst[7], dst[12], idx[7], idx[12]);
  SORT_CSWAP_PAIR(dst[5], dst[9], idx[5], idx[9]);
  SORT_CSWAP_PAIR(dst[0], dst[2], idx[0], idx[2]);
  SORT_CSWAP_PAIR(dst[3], dst[7], idx[3], idx[7]);
  SORT_CSWAP_PAIR(dst[10], dst[11], idx[10], idx[11]);
  SORT_CSWAP_PAIR(dst[1], dst[4], idx[1], idx[4]);
  SORT_CSWAP_PAIR(dst[6], dst[12], idx[6], idx[12]);
  SORT_CSWAP_PAIR(dst[7], dst[8], idx[7], idx[8]);
  SORT_CSWAP_PAIR(dst[11], dst[12], idx[11], idx[12]);
  SORT_CSWAP_PAIR(dst[4], dst[9], idx[4], idx[9]);
  SORT_CSWAP_PAIR(dst[6], dst[10], idx[6], idx[10]);
  SORT_CSWAP_PAIR(dst[3], dst[4], idx[3], idx[4]);
  SORT_CSWAP_PAIR(dst[5], dst[6], idx[5], idx[6]);
  SORT_CSWAP_PAIR(dst[8], dst[9], idx[8], idx[9]);
  SORT_CSWAP_PAIR(dst[10], dst[11], idx[10], idx[11]);
  SORT_CSWAP_PAIR(dst[1], dst[7], idx[1], idx[7]);
  SORT_CSWAP_PAIR(dst[2], dst[6], idx[2], idx[6]);
  SORT_CSWAP_PAIR(dst[9], dst[11], idx[9], idx[11]);
  SORT_CSWAP_PAIR(dst[1], dst[3], idx[1], idx[3]);
  SORT_CSWAP_PAIR(dst[4], dst[7], idx[4], idx[7]);
  SORT_CSWAP_PAIR(dst[8], dst[10], idx[8], idx[10]);
  SORT_CSWAP_PAIR(dst[0], dst[5], idx[0], idx[5]);
  SORT_CSWAP_PAIR(dst[2], dst[5], idx[2], idx[5]);
  SORT_CSWAP_PAIR(dst[6], dst[8], idx[6], idx[8]);
  SORT_CSWAP_PAIR(dst[9], dst[10], idx[9], idx[10]);
  SORT_CSWAP_PAIR(dst[1], dst[2], idx[1], idx[2]);
  SORT_CSWAP_PAIR(dst[3], dst[5], idx[3], idx[5]);
  SORT_CSWAP_PAIR(dst[7], dst[8], idx[7], idx[8]);
  SORT_CSWAP_PAIR(dst[4], dst[6], idx[4], idx[6]);
  SORT_CSWAP_PAIR(dst[2], dst[3], idx[2], idx[3]);
  SORT_CSWAP_PAIR(dst[4], dst[5], idx[4], idx[5]);
  SORT_CSWAP_PAIR(dst[6], dst[7], idx[6], idx[7]);
  SORT_CSWAP_PAIR(dst[8], dst[9], idx[8], idx[9]);
  SORT_CSWAP_PAIR(dst[3], dst[4], idx[3], idx[4]);
  SORT_CSWAP_PAIR(dst[5], dst[6], idx[5], idx[6]);
}

#define BITONIC_SORT_PAIR_14_PAIR SORT_MAKE_STR(bitonic_sort_14_pair)
static __inline void BITONIC_SORT_PAIR_14(size_t *idx, SORT_TYPE *dst) {
  SORT_CSWAP_PAIR(dst[0], dst[1], idx[0], idx[1]);
  SORT_CSWAP_PAIR(dst[2], dst[3], idx[2], idx[3]);
  SORT_CSWAP_PAIR(dst[4], dst[5], idx[4], idx[5]);
  SORT_CSWAP_PAIR(dst[6], dst[7], idx[6], idx[7]);
  SORT_CSWAP_PAIR(dst[8], dst[9], idx[8], idx[9]);
  SORT_CSWAP_PAIR(dst[10], dst[11], idx[10], idx[11]);
  SORT_CSWAP_PAIR(dst[12], dst[13], idx[12], idx[13]);
  SORT_CSWAP_PAIR(dst[0], dst[2], idx[0], idx[2]);
  SORT_CSWAP_PAIR(dst[4], dst[6], idx[4], idx[6]);
  SORT_CSWAP_PAIR(dst[8], dst[10], idx[8], idx[10]);
  SORT_CSWAP_PAIR(dst[1], dst[3], idx[1], idx[3]);
  SORT_CSWAP_PAIR(dst[5], dst[7], idx[5], idx[7]);
  SORT_CSWAP_PAIR(dst[9], dst[11], idx[9], idx[11]);
  SORT_CSWAP_PAIR(dst[0], dst[4], idx[0], idx[4]);
  SORT_CSWAP_PAIR(dst[8], dst[12], idx[8], idx[12]);
  SORT_CSWAP_PAIR(dst[1], dst[5], idx[1], idx[5]);
  SORT_CSWAP_PAIR(dst[9], dst[13], idx[9], idx[13]);
  SORT_CSWAP_PAIR(dst[2], dst[6], idx[2], idx[6]);
  SORT_CSWAP_PAIR(dst[3], dst[7], idx[3], idx[7]);
  SORT_CSWAP_PAIR(dst[0], dst[8], idx[0], idx[8]);
  SORT_CSWAP_PAIR(dst[1], dst[9], idx[1], idx[9]);
  SORT_CSWAP_PAIR(dst[2], dst[10], idx[2], idx[10]);
  SORT_CSWAP_PAIR(dst[3], dst[11], idx[3], idx[11]);
  SORT_CSWAP_PAIR(dst[4], dst[12], idx[4], idx[12]);
  SORT_CSWAP_PAIR(dst[5], dst[13], idx[5], idx[13]);
  SORT_CSWAP_PAIR(dst[5], dst[10], idx[5], idx[10]);
  SORT_CSWAP_PAIR(dst[6], dst[9], idx[6], idx[9]);
  SORT_CSWAP_PAIR(dst[3], dst[12], idx[3], idx[12]);
  SORT_CSWAP_PAIR(dst[7], dst[11], idx[7], idx[11]);
  SORT_CSWAP_PAIR(dst[1], dst[2], idx[1], idx[2]);
  SORT_CSWAP_PAIR(dst[4], dst[8], idx[4], idx[8]);
  SORT_CSWAP_PAIR(dst[1], dst[4], idx[1], idx[4]);
  SORT_CSWAP_PAIR(dst[7], dst[13], idx[7], idx[13]);
  SORT_CSWAP_PAIR(dst[2], dst[8], idx[2], idx[8]);
  SORT_CSWAP_PAIR(dst[5], dst[6], idx[5], idx[6]);
  SORT_CSWAP_PAIR(dst[9], dst[10], idx[9], idx[10]);
  SORT_CSWAP_PAIR(dst[2], dst[4], idx[2], idx[4]);
  SORT_CSWAP_PAIR(dst[11], dst[13], idx[11], idx[13]);
  SORT_CSWAP_PAIR(dst[3], dst[8], idx[3], idx[8]);
  SORT_CSWAP_PAIR(dst[7], dst[12], idx[7], idx[12]);
  SORT_CSWAP_PAIR(dst[6], dst[8], idx[6], idx[8]);
  SORT_CSWAP_PAIR(dst[10], dst[12], idx[10], idx[12]);
  SORT_CSWAP_PAIR(dst[3], dst[5], idx[3], idx[5]);
  SORT_CSWAP_PAIR(dst[7], dst[9], idx[7], idx[9]);
  SORT_CSWAP_PAIR(dst[3], dst[4], idx[3], idx[4]);
  SORT_CSWAP_PAIR(dst[5], dst[6], idx[5], idx[6]);
  SORT_CSWAP_PAIR(dst[7], dst[8], idx[7], idx[8]);
  SORT_CSWAP_PAIR(dst[9], dst[10], idx[9], idx[10]);
  SORT_CSWAP_PAIR(dst[11], dst[12], idx[11], idx[12]);
  SORT_CSWAP_PAIR(dst[6], dst[7], idx[6], idx[7]);
  SORT_CSWAP_PAIR(dst[8], dst[9], idx[8], idx[9]);
}

#define BITONIC_SORT_PAIR_15_PAIR SORT_MAKE_STR(bitonic_sort_15_pair)
static __inline void BITONIC_SORT_PAIR_15(size_t *idx, SORT_TYPE *dst) {
  SORT_CSWAP_PAIR(dst[0], dst[1], idx[0], idx[1]);
  SORT_CSWAP_PAIR(dst[2], dst[3], idx[2], idx[3]);
  SORT_CSWAP_PAIR(dst[4], dst[5], idx[4], idx[5]);
  SORT_CSWAP_PAIR(dst[6], dst[7], idx[6], idx[7]);
  SORT_CSWAP_PAIR(dst[8], dst[9], idx[8], idx[9]);
  SORT_CSWAP_PAIR(dst[10], dst[11], idx[10], idx[11]);
  SORT_CSWAP_PAIR(dst[12], dst[13], idx[12], idx[13]);
  SORT_CSWAP_PAIR(dst[0], dst[2], idx[0], idx[2]);
  SORT_CSWAP_PAIR(dst[4], dst[6], idx[4], idx[6]);
  SORT_CSWAP_PAIR(dst[8], dst[10], idx[8], idx[10]);
  SORT_CSWAP_PAIR(dst[12], dst[14], idx[12], idx[14]);
  SORT_CSWAP_PAIR(dst[1], dst[3], idx[1], idx[3]);
  SORT_CSWAP_PAIR(dst[5], dst[7], idx[5], idx[7]);
  SORT_CSWAP_PAIR(dst[9], dst[11], idx[9], idx[11]);
  SORT_CSWAP_PAIR(dst[0], dst[4], idx[0], idx[4]);
  SORT_CSWAP_PAIR(dst[8], dst[12], idx[8], idx[12]);
  SORT_CSWAP_PAIR(dst[1], dst[5], idx[1], idx[5]);
  SORT_CSWAP_PAIR(dst[9], dst[13], idx[9], idx[13]);
  SORT_CSWAP_PAIR(dst[2], dst[6], idx[2], idx[6]);
  SORT_CSWAP_PAIR(dst[10], dst[14], idx[10], idx[14]);
  SORT_CSWAP_PAIR(dst[3], dst[7], idx[3], idx[7]);
  SORT_CSWAP_PAIR(dst[0], dst[8], idx[0], idx[8]);
  SORT_CSWAP_PAIR(dst[1], dst[9], idx[1], idx[9]);
  SORT_CSWAP_PAIR(dst[2], dst[10], idx[2], idx[10]);
  SORT_CSWAP_PAIR(dst[3], dst[11], idx[3], idx[11]);
  SORT_CSWAP_PAIR(dst[4], dst[12], idx[4], idx[12]);
  SORT_CSWAP_PAIR(dst[5], dst[13], idx[5], idx[13]);
  SORT_CSWAP_PAIR(dst[6], dst[14], idx[6], idx[14]);
  SORT_CSWAP_PAIR(dst[5], dst[10], idx[5], idx[10]);
  SORT_CSWAP_PAIR(dst[6], dst[9], idx[6], idx[9]);
  SORT_CSWAP_PAIR(dst[3], dst[12], idx[3], idx[12]);
  SORT_CSWAP_PAIR(dst[13], dst[14], idx[13], idx[14]);
  SORT_CSWAP_PAIR(dst[7], dst[11], idx[7], idx[11]);
  SORT_CSWAP_PAIR(dst[1], dst[2], idx[1], idx[2]);
  SORT_CSWAP_PAIR(dst[4], dst[8], idx[4], idx[8]);
  SORT_CSWAP_PAIR(dst[1], dst[4], idx[1], idx[4]);
  SORT_CSWAP_PAIR(dst[7], dst[13], idx[7], idx[13]);
  SORT_CSWAP_PAIR(dst[2], dst[8], idx[2], idx[8]);
  SORT_CSWAP_PAIR(dst[11], dst[14], idx[11], idx[14]);
  SORT_CSWAP_PAIR(dst[5], dst[6], idx[5], idx[6]);
  SORT_CSWAP_PAIR(dst[9], dst[10], idx[9], idx[10]);
  SORT_CSWAP_PAIR(dst[2], dst[4], idx[2], idx[4]);
  SORT_CSWAP_PAIR(dst[11], dst[13], idx[11], idx[13]);
  SORT_CSWAP_PAIR(dst[3], dst[8], idx[3], idx[8]);
  SORT_CSWAP_PAIR(dst[7], dst[12], idx[7], idx[12]);
  SORT_CSWAP_PAIR(dst[6], dst[8], idx[6], idx[8]);
  SORT_CSWAP_PAIR(dst[10], dst[12], idx[10], idx[12]);
  SORT_CSWAP_PAIR(dst[3], dst[5], idx[3], idx[5]);
  SORT_CSWAP_PAIR(dst[7], dst[9], idx[7], idx[9]);
  SORT_CSWAP_PAIR(dst[3], dst[4], idx[3], idx[4]);
  SORT_CSWAP_PAIR(dst[5], dst[6], idx[5], idx[6]);
  SORT_CSWAP_PAIR(dst[7], dst[8], idx[7], idx[8]);
  SORT_CSWAP_PAIR(dst[9], dst[10], idx[9], idx[10]);
  SORT_CSWAP_PAIR(dst[11], dst[12], idx[11], idx[12]);
  SORT_CSWAP_PAIR(dst[6], dst[7], idx[6], idx[7]);
  SORT_CSWAP_PAIR(dst[8], dst[9], idx[8], idx[9]);
}

#define BITONIC_SORT_PAIR_16 SORT_MAKE_STR(bitonic_sort_16_pair)
static __inline void BITONIC_SORT_PAIR_16(size_t *idx, SORT_TYPE *dst) {
  SORT_CSWAP_PAIR(dst[0], dst[1], idx[0], idx[1]);
  SORT_CSWAP_PAIR(dst[2], dst[3], idx[2], idx[3]);
  SORT_CSWAP_PAIR(dst[4], dst[5], idx[4], idx[5]);
  SORT_CSWAP_PAIR(dst[6], dst[7], idx[6], idx[7]);
  SORT_CSWAP_PAIR(dst[8], dst[9], idx[8], idx[9]);
  SORT_CSWAP_PAIR(dst[10], dst[11], idx[10], idx[11]);
  SORT_CSWAP_PAIR(dst[12], dst[13], idx[12], idx[13]);
  SORT_CSWAP_PAIR(dst[14], dst[15], idx[14], idx[15]);
  SORT_CSWAP_PAIR(dst[0], dst[2], idx[0], idx[2]);
  SORT_CSWAP_PAIR(dst[4], dst[6], idx[4], idx[6]);
  SORT_CSWAP_PAIR(dst[8], dst[10], idx[8], idx[10]);
  SORT_CSWAP_PAIR(dst[12], dst[14], idx[12], idx[14]);
  SORT_CSWAP_PAIR(dst[1], dst[3], idx[1], idx[3]);
  SORT_CSWAP_PAIR(dst[5], dst[7], idx[5], idx[7]);
  SORT_CSWAP_PAIR(dst[9], dst[11], idx[9], idx[11]);
  SORT_CSWAP_PAIR(dst[13], dst[15], idx[13], idx[15]);
  SORT_CSWAP_PAIR(dst[0], dst[4], idx[0], idx[4]);
  SORT_CSWAP_PAIR(dst[8], dst[12], idx[8], idx[12]);
  SORT_CSWAP_PAIR(dst[1], dst[5], idx[1], idx[5]);
  SORT_CSWAP_PAIR(dst[9], dst[13], idx[9], idx[13]);
  SORT_CSWAP_PAIR(dst[2], dst[6], idx[2], idx[6]);
  SORT_CSWAP_PAIR(dst[10], dst[14], idx[10], idx[14]);
  SORT_CSWAP_PAIR(dst[3], dst[7], idx[3], idx[7]);
  SORT_CSWAP_PAIR(dst[11], dst[15], idx[11], idx[15]);
  SORT_CSWAP_PAIR(dst[0], dst[8], idx[0], idx[8]);
  SORT_CSWAP_PAIR(dst[1], dst[9], idx[1], idx[9]);
  SORT_CSWAP_PAIR(dst[2], dst[10], idx[2], idx[10]);
  SORT_CSWAP_PAIR(dst[3], dst[11], idx[3], idx[11]);
  SORT_CSWAP_PAIR(dst[4], dst[12], idx[4], idx[12]);
  SORT_CSWAP_PAIR(dst[5], dst[13], idx[5], idx[13]);
  SORT_CSWAP_PAIR(dst[6], dst[14], idx[6], idx[14]);
  SORT_CSWAP_PAIR(dst[7], dst[15], idx[7], idx[15]);
  SORT_CSWAP_PAIR(dst[5], dst[10], idx[5], idx[10]);
  SORT_CSWAP_PAIR(dst[6], dst[9], idx[6], idx[9]);
  SORT_CSWAP_PAIR(dst[3], dst[12], idx[3], idx[12]);
  SORT_CSWAP_PAIR(dst[13], dst[14], idx[13], idx[14]);
  SORT_CSWAP_PAIR(dst[7], dst[11], idx[7], idx[11]);
  SORT_CSWAP_PAIR(dst[1], dst[2], idx[1], idx[2]);
  SORT_CSWAP_PAIR(dst[4], dst[8], idx[4], idx[8]);
  SORT_CSWAP_PAIR(dst[1], dst[4], idx[1], idx[4]);
  SORT_CSWAP_PAIR(dst[7], dst[13], idx[7], idx[13]);
  SORT_CSWAP_PAIR(dst[2], dst[8], idx[2], idx[8]);
  SORT_CSWAP_PAIR(dst[11], dst[14], idx[11], idx[14]);
  SORT_CSWAP_PAIR(dst[5], dst[6], idx[5], idx[6]);
  SORT_CSWAP_PAIR(dst[9], dst[10], idx[9], idx[10]);
  SORT_CSWAP_PAIR(dst[2], dst[4], idx[2], idx[4]);
  SORT_CSWAP_PAIR(dst[11], dst[13], idx[11], idx[13]);
  SORT_CSWAP_PAIR(dst[3], dst[8], idx[3], idx[8]);
  SORT_CSWAP_PAIR(dst[7], dst[12], idx[7], idx[12]);
  SORT_CSWAP_PAIR(dst[6], dst[8], idx[6], idx[8]);
  SORT_CSWAP_PAIR(dst[10], dst[12], idx[10], idx[12]);
  SORT_CSWAP_PAIR(dst[3], dst[5], idx[3], idx[5]);
  SORT_CSWAP_PAIR(dst[7], dst[9], idx[7], idx[9]);
  SORT_CSWAP_PAIR(dst[3], dst[4], idx[3], idx[4]);
  SORT_CSWAP_PAIR(dst[5], dst[6], idx[5], idx[6]);
  SORT_CSWAP_PAIR(dst[7], dst[8], idx[7], idx[8]);
  SORT_CSWAP_PAIR(dst[9], dst[10], idx[9], idx[10]);
  SORT_CSWAP_PAIR(dst[11], dst[12], idx[11], idx[12]);
  SORT_CSWAP_PAIR(dst[6], dst[7], idx[6], idx[7]);
  SORT_CSWAP_PAIR(dst[8], dst[9], idx[8], idx[9]);
}

void BITONIC_SORT_PAIR(size_t *idx, SORT_TYPE *dst, const size_t size) {
  switch (size) {
    case 0:
    case 1:
      break;

    case 2:
      BITONIC_SORT_PAIR_2(idx, dst);
      break;

    case 3:
      BITONIC_SORT_PAIR_3(idx, dst);
      break;

    case 4:
      BITONIC_SORT_PAIR_4(idx, dst);
      break;

    case 5:
      BITONIC_SORT_PAIR_5(idx, dst);
      break;

    case 6:
      BITONIC_SORT_PAIR_6(idx, dst);
      break;

    case 7:
      BITONIC_SORT_PAIR_7(idx, dst);
      break;

    case 8:
      BITONIC_SORT_PAIR_8(idx, dst);
      break;

    case 9:
      BITONIC_SORT_PAIR_9(idx, dst);
      break;

    case 10:
      BITONIC_SORT_PAIR_10(idx, dst);
      break;

    case 11:
      BITONIC_SORT_PAIR_11(idx, dst);
      break;

    case 12:
      BITONIC_SORT_PAIR_12(idx, dst);
      break;

    case 13:
      BITONIC_SORT_PAIR_13(idx, dst);
      break;

    case 14:
      BITONIC_SORT_PAIR_14(idx, dst);
      break;

    case 15:
      BITONIC_SORT_PAIR_15(idx, dst);
      break;

    case 16:
      BITONIC_SORT_PAIR_16(idx, dst);
      break;

    default:
      BINARY_INSERTION_SORT_PAIR(idx, dst, size);
  }
}

/* Function used to do a binary search for binary insertion sort */
static __inline size_t BINARY_INSERTION_FIND(SORT_TYPE *dst, const SORT_TYPE x,
                                             const size_t size) {
  size_t l, c, r;
  SORT_TYPE cx;
  l = 0;
  r = size - 1;
  c = r >> 1;

  /* check for out of bounds at the beginning. */
  if (SORT_CMP(x, dst[0]) < 0) {
    return 0;
  } else if (SORT_CMP(x, dst[r]) > 0) {
    return r;
  }

  cx = dst[c];

  while (1) {
    const int val = SORT_CMP(x, cx);

    if (val < 0) {
      if (c - l <= 1) {
        return c;
      }

      r = c;
    } else { /* allow = for stability. The binary search favors the right. */
      if (r - c <= 1) {
        return c + 1;
      }

      l = c;
    }

    c = l + ((r - l) >> 1);
    cx = dst[c];
  }
}

/* Binary insertion sort, but knowing that the first "start" entries are sorted.
 * Used in timsort. */
static void BINARY_INSERTION_SORT_START_PAIR(size_t *idx, SORT_TYPE *dst,
                                             const size_t start,
                                             const size_t size) {
  size_t i;

  for (i = start; i < size; i++) {
    size_t j;
    SORT_TYPE x;
    size_t idxi;
    size_t location;

    /* If this entry is already correct, just move along */
    if (SORT_CMP(dst[i - 1], dst[i]) <= 0) {
      continue;
    }

    /* Else we need to find the right place, shift everything over, and
     * squeeze in */
    x = dst[i];
    idxi = idx[i];
    location = BINARY_INSERTION_FIND(dst, x, i);

    for (j = i - 1; j >= location; j--) {
      dst[j + 1] = dst[j];
      idx[j + 1] = idx[j];

      if (j == 0) { /* check edge case because j is unsigned */
        break;
      }
    }

    dst[location] = x;
    idx[location] = idxi;
  }
}

/* Binary insertion sort */
void BINARY_INSERTION_SORT_PAIR(size_t *idx, SORT_TYPE *dst,
                                const size_t size) {
  /* don't bother sorting an array of size <= 1 */
  if (size <= 1) {
    return;
  }

  BINARY_INSERTION_SORT_START_PAIR(idx, dst, 1, size);
}

/* Quick sort: based on wikipedia */

static __inline size_t QUICK_SORT_PARTITION_PAIR(size_t *idx, SORT_TYPE *dst,
                                                 const size_t left,
                                                 const size_t right,
                                                 const size_t pivot) {
  SORT_TYPE value = dst[pivot];
  size_t index = left;
  size_t i;
  int not_all_same = 0;
  /* move the pivot to the right */
  SORT_SWAP(dst[pivot], dst[right]);
  IDX_SWAP(idx[pivot], idx[right]);

  for (i = left; i < right; i++) {
    int cmp = SORT_CMP(dst[i], value);
    /* check if everything is all the same */
    not_all_same |= cmp;

    if (cmp < 0) {
      SORT_SWAP(dst[i], dst[index]);
      IDX_SWAP(idx[i], idx[index]);
      index++;
    }
  }

  SORT_SWAP(dst[right], dst[index]);
  IDX_SWAP(idx[right], idx[index]);

  /* avoid degenerate case */
  if (not_all_same == 0) {
    return SIZE_MAX;
  }

  return index;
}

static void QUICK_SORT_RECURSIVE_PAIR(size_t *idx, SORT_TYPE *dst,
                                      const size_t left, const size_t right) {
  size_t pivot;
  size_t new_pivot;

  if (right <= left) {
    return;
  }

  if ((right - left + 1U) <= SMALL_SORT_BND) {
    SMALL_SORT_PAIR(&idx[left], &dst[left], right - left + 1U);
    return;
  }

  pivot = left + ((right - left) >> 1);
  /* this seems to perform worse by a small amount... ? */
  /* pivot = MEDIAN(dst, left, pivot, right); */
  new_pivot = QUICK_SORT_PARTITION_PAIR(idx, dst, left, right, pivot);

  /* check for partition all equal */
  if (new_pivot == SIZE_MAX) {
    return;
  }

  QUICK_SORT_RECURSIVE_PAIR(idx, dst, left, new_pivot - 1U);
  QUICK_SORT_RECURSIVE_PAIR(idx, dst, new_pivot + 1U, right);
}

void QUICK_SORT_PAIR(size_t *idx, SORT_TYPE *dst, const size_t size) {
  /* don't bother sorting an array of size 1 */
  if (size <= 1) {
    return;
  }

  QUICK_SORT_RECURSIVE_PAIR(idx, dst, 0U, size - 1U);
}

#undef BITONIC_SORT_PAIR
#undef QUICK_SORT_PAIR
#undef QUICK_SORT_PARTITION_PAIR
#undef QUICK_SORT_RECURSIVE_PAIR
#undef BINARY_INSERTION_SORT_START_PAIR
#undef BINARY_INSERTION_SORT_PAIR
