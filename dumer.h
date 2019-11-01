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
#include <m4ri/config.h>
#include <m4ri/m4ri.h>
#include <omp.h>
#include <stdint.h>

#ifndef DUMER_L
#define DUMER_L 16
#endif
#ifndef DUMER_P
#define DUMER_P 4
#endif
#define DUMER_P1 (DUMER_P / 2)
#define DUMER_P2 (DUMER_P - DUMER_P1)
#ifndef DUMER_EPS
#define DUMER_EPS 40
#endif
#ifndef DUMER_DOOM
#define DUMER_DOOM 0
#endif
#ifndef DUMER_LW
#define DUMER_LW 0
#endif
#ifndef DUMER_LUT
#define DUMER_LUT 11
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

#if DUMER_L <= 16
#define LIST_WIDTH 16
#define LIST_TYPE uint16_t
#define SORT_TYPE uint16_t
#elif DUMER_L <= 32
#define LIST_WIDTH 32
#define LIST_TYPE uint32_t
#define SORT_TYPE uint32_t
#elif DUMER_L <= 64
#define LIST_WIDTH 64
#define LIST_TYPE uint64_t
#define SORT_TYPE uint64_t
#endif

enum type { QC, SD, LW, GO };

struct shared {
  uint16_t *list1_pos;

  uint16_t *combinations2;
  uint16_t *combinations2_diff;

  uint64_t nb_combinations1;
  uint64_t nb_combinations2;

#if DUMER_LW
  omp_lock_t w_best_lock;
  int w_best;
#endif
};

struct isd {
  mzd_t *A;

  int *perm;
  /* Seeds for pseudo random number generator. */
  uint64_t S0;
  uint64_t S1;

  size_t size_list1;
  LIST_TYPE *list1;
  size_t *list1_idx;
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

  int w_target;
  int w_solution;
  uint8_t *solution;

  /* Avoid mallocing and freeing all the time. */

  /* Arrays used when building list1. */
  LIST_TYPE *scratch;
  LIST_TYPE *stack_syndrome;
  uint16_t *stack_nb_flips;
  uint16_t *stack_maxval;

  /* Arrays used during collision search. */
  uint64_t *test_syndrome;
#if DUMER_DOOM || DUMER_LW
  uint64_t *current_nosyndrome;
#endif
  uint64_t *current_syndrome;
  uint64_t *xor_pairs;
};

typedef struct isd *isd_t;
typedef struct shared *shr_t;

shr_t alloc_shr(int n1, int n2);
void free_shr(shr_t shr);
void init_shr(shr_t shr, int n1, int n2);
isd_t alloc_isd(int n, int k, int r, int n1, int n2, uint64_t nb_combinations1);
void free_isd(isd_t isd);
void init_isd(isd_t isd, enum type current_type, int n, int k, int w,
              int *mat_h, int *mat_s);

int dumer(int n, int k, int r, int n1, int n2, shr_t shr, isd_t isd);
void print_solution(int n, isd_t isd);
#endif /* DUMER_H */
