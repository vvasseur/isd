/*
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
#ifndef DUMER_H
#define DUMER_H
#include <omp.h>

#define TOKEN_CAT(x, y) x##y
#define XTOKEN_CAT(x, y) TOKEN_CAT(x, y)

#include "matrix.h"

#ifndef DUMER_L
#define DUMER_L 16L
#endif
#ifndef DUMER_P
#define DUMER_P 4L
#endif
#define DUMER_P1 (DUMER_P / 2)
#define DUMER_P2 (DUMER_P - DUMER_P1)
#ifndef DUMER_EPS
#define DUMER_EPS 40L
#endif
#ifndef DUMER_DOOM
#define DUMER_DOOM 0
#endif
#ifndef DUMER_LW
#define DUMER_LW 0
#endif
#ifndef DUMER_LUT
#define DUMER_LUT 11L
#endif
#if (DUMER_LUT) > (DUMER_L)
#undef DUMER_LUT
#define DUMER_LUT (DUMER_L)
#endif
#ifndef DUMER_LUT_SHIFT
#define DUMER_LUT_SHIFT ((DUMER_L) - (DUMER_LUT))
#endif
#if DUMER_P < 4 || DUMER_P > 8
#error "No implementation for this value of DUMER_P"
#endif
#if DUMER_L < 1
#error "DUMER_L should be greater than 0"
#endif
#if DUMER_L > 64
#error "No implementation for this value of DUMER_L"
#endif

#if DUMER_L <= 8
#define LIST_WIDTH 8
#define SORT_WIDTH 8
#define LIST_TYPE uint8_t
#define SORT_TYPE uint8_t
#elif DUMER_L <= 16
#define LIST_WIDTH 16
#define SORT_WIDTH 16
#define LIST_TYPE uint16_t
#define SORT_TYPE uint16_t
#elif DUMER_L <= 32
#define LIST_WIDTH 32
#define SORT_WIDTH 32
#define LIST_TYPE uint32_t
#define SORT_TYPE uint32_t
#elif DUMER_L <= 64
#define LIST_WIDTH 64
#define SORT_WIDTH 64
#define LIST_TYPE uint64_t
#define SORT_TYPE uint64_t
#endif
#if DUMER_L == 64
#define DUMER_L_MASK (~0UL)
#else
#define DUMER_L_MASK ((uint64_t)((1UL << DUMER_L) - 1))
#endif
#define xor_bcast XTOKEN_CAT(xor_bcast_, LIST_WIDTH)

enum type { QC, SD, LW, GO };

struct shared {
  uint16_t *list1_pos;

  uint16_t *combinations2;
  uint16_t *combinations2_diff;

  uint64_t nb_combinations1;
  uint64_t nb_combinations2;

#if DUMER_LW
  omp_lock_t w_best_lock;
  size_t w_best;
#endif
  int **gray_rev;
  int **gray_diff;
  size_t k_opt;
};

struct isd {
  matrix_t A;
  matrix_t At;

  size_t *perm;
  /* Seeds for pseudo random number generator. */
  uint64_t S0;
  uint64_t S1;

  size_t size_list1;
  LIST_TYPE *list1;
  LIST_TYPE *list1_aux;
  size_t *list1_idx;
  size_t *list1_aux2;
  size_t *list1_lut;

  size_t size_columns1_low;
  LIST_TYPE *columns1_low;

  size_t size_columns1_full;
  uint64_t *columns1_full;

  size_t size_columns2_full;
  uint64_t *columns2_full;

#if !(DUMER_LW)
  uint64_t *s_full;
#endif

  size_t w_target;
  size_t w_solution;
  uint8_t *solution;

  /* Avoid mallocing and freeing all the time. */

  /* Arrays used when building list1. */
  uint8_t *scratch;

  /* Arrays used during collision search. */
  uint64_t *test_syndrome;
#if DUMER_DOOM || DUMER_LW
  uint64_t *current_nosyndrome;
#endif
  uint64_t *current_syndrome;
  uint64_t *xor_pairs;

  uint64_t *xor_rows;
};

typedef struct isd *isd_t;
typedef struct shared *shr_t;

shr_t alloc_shr(size_t n1, size_t n2);
void free_shr(shr_t shr);
void init_shr(shr_t shr, size_t n, size_t k, size_t n1, size_t n2);
isd_t alloc_isd(size_t n, size_t k, size_t r, size_t n1, size_t n2,
                uint64_t nb_combinations1, size_t k_opt);
void free_isd(isd_t isd, size_t r, size_t n);
void init_isd(isd_t isd, enum type current_type, size_t n, size_t k, size_t w,
              uint8_t *mat_h, uint8_t *mat_s);

int dumer(size_t n, size_t k, size_t r, size_t n1, size_t n2, shr_t shr,
          isd_t isd);
void print_solution(size_t n, isd_t isd);
#endif /* DUMER_H */
