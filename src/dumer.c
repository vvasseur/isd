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
#include "dumer.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bits.h"
#include "matrix.h"
#include "sort.h"
#include "transpose.h"
#include "xoroshiro128plus.h"

/* Binomial coefficient. */
uint64_t bincoef(size_t n, size_t k) {
  uint64_t res = 1;
  for (size_t i = 0; i < k; ++i) {
    res *= (n - i);
    res /= (i + 1);
  }
  return res;
}

/*
 * List all combinations of 't' elements of a set of 'n' elements.
 *
 * Generate a Chase's sequence: the binary representation of a combination and
 * its successor only differ by two bits that are either consecutive of
 * separated by only one position.
 *
 * See exercise 45 of Knuth's The art of computer programming volume 4A.
 */
static void chase(size_t n, size_t t, uint16_t *combinations, uint16_t *diff) {
  size_t N = 0;
  uint16_t diff_pos = 0;
  uint16_t diff_len = 0;
  int32_t x;
  uint16_t *c = malloc((t + 2) * sizeof(uint16_t));
  uint16_t *z = malloc((t + 2) * sizeof(uint16_t));
  for (size_t j = 1; j <= t + 1; ++j) {
    z[j] = 0;
  }
  for (size_t j = 1; j <= t + 1; ++j) {
    c[j] = n - t - 1 + j;
  }
  /* r is the least subscript with c[r] >= r. */
  size_t r = 1;
  size_t j;

  while (1) {
    for (size_t i = 1; i <= t; ++i) {
      combinations[i - 1 + N * t] = c[i];
    }
    diff[N] = diff_pos + (diff_len - 1) * (n - 1);
    ++N;
    j = r;

  novisit:
    if (z[j]) {
      x = c[j] + 2;
      if (x < z[j]) {
        diff_pos = c[j];
        diff_len = 2;
        c[j] = x;
      } else if (x == z[j] && z[j + 1]) {
        diff_pos = c[j];
        diff_len = 2 - (c[j + 1] % 2);
        c[j] = x - (c[j + 1] % 2);
      } else {
        z[j] = 0;
        ++j;
        if (j <= t)
          goto novisit;
        else
          goto exit;
      }
      if (c[1] > 0) {
        r = 1;
      } else {
        r = j - 1;
      }
    } else {
      x = c[j] + (c[j] % 2) - 2;
      if (x >= (int32_t)j) {
        diff_pos = x;
        diff_len = 2 - (c[j] % 2);
        c[j] = x;
        r = 1;
      } else if (c[j] == j) {
        diff_pos = j - 1;
        diff_len = 1;
        c[j] = j - 1;
        z[j] = c[j + 1] - ((c[j + 1] + 1) % 2);
        r = j;
      } else if (c[j] < j) {
        diff_pos = c[j];
        diff_len = j - c[j];
        c[j] = j;
        z[j] = c[j + 1] - ((c[j + 1] + 1) % 2);
        r = (j > 2) ? j - 1 : 1;
      } else {
        diff_pos = x;
        diff_len = 2 - (c[j] % 2);
        c[j] = x;
        r = j;
      }
    }
  }

exit:
  free(c);
  free(z);
}

/* Apply the same permutation to a matrix and an array. */
static void shuffle_matrix(matrix_t At, size_t *perm, size_t n, size_t n_stop,
                           uint64_t *S0, uint64_t *S1) {
  uint8_t *state = calloc(n, sizeof(uint8_t));

  size_t weight = 0;
  while (weight < n_stop) {
    uint32_t rand = random_lim(n - 1, S0, S1);
    weight += !state[rand];
    state[rand] = 1;
  }

  /* Swap the columns so that all those marked with a one are on the left-hand
   * side of the matrix. */
  size_t i = 0;
  size_t j = n_stop;
  while (i < n_stop && j < n) {
    if (!state[i] && state[j]) {
      matrix_swap_rows(At, i, j);
      size_t swp = perm[i];
      perm[i] = perm[j];
      perm[j] = swp;

      ++i;
      ++j;
    } else {
      if (state[i]) ++i;
      if (!state[j]) ++j;
    }
  }

  free(state);
}

/* Randomly choose an information set and perform a Gaussian elimination. */
static void choose_is(matrix_t A, matrix_t At, size_t *perm, size_t n, size_t k,
                      size_t l, size_t k_opt, int **rev, int **diff,
                      uint64_t *xor_rows, uint64_t *S0, uint64_t *S1) {
  /* Pick a permutation and perform Gaussian elimination.  */
  size_t r = 0;
  while (r < n - k - l) {
    shuffle_matrix(At, perm, n, n - k - l, S0, S1);
#if DUMER_LW
    matrix_transpose_rev_rows(A, At, n, n - k);
#elif !(DUMER_DOOM)  // && !(DUMER_LW)
    matrix_transpose_rev_rows(A, At, n + 1, n - k);
#else                // DUMER_DOOM && !(DUMER_LW)
    matrix_transpose_rev_rows(A, At, n + k, n - k);
#endif

    r = matrix_echelonize_partial(A, n - k, n, k_opt, n - k - l, xor_rows, rev,
                                  diff);
  }
}

/* Extract columns from *A and keep data 32-byte aligned (fitting AVX
 * registers). */
static void get_columns_H_prime_avx(matrix_t A, uint64_t *columns, size_t n,
                                    size_t r) {
  size_t stride = AVX_PADDING(r) / 256;

  for (size_t j = 0; j < n; ++j) {
    copy_avx((uint8_t *)columns, (uint8_t *)(A[j]), stride);
    columns += stride * 4;
  }
}

/* Extract columns from *A and keep data LIST_WIDTH-byte aligned. */
static void get_columns_H_prime(matrix_t A, LIST_TYPE *columns, size_t n) {
  for (size_t j = 0; j < n; ++j) {
    columns[j] = A[j][0] & DUMER_L_MASK;
  }
}

/*
 * See Paul Khuong's
 * https://www.pvk.ca/Blog/2012/07/03/binary-search-star-eliminates-star-branch-mispredictions/
 */
static size_t bin_search(const LIST_TYPE *list, size_t len_list,
                         LIST_TYPE value) {
  if (len_list <= 1) return 0;
  unsigned log = clb(len_list) - 1;
  size_t first_mid = len_list - (1UL << log);
  const LIST_TYPE *low = (list[first_mid] < value) ? list + first_mid : list;
  len_list = 1UL << log;

  for (unsigned i = log; i != 0; i--) {
    len_list /= 2;
    LIST_TYPE mid = low[len_list];
    if (mid < value) low += len_list;
  }

  return (*low == value) ? low - list : low - list + 1;
}

/*
 * Build a list containing the XORs of all possible combinations of 'p'
 * columns.
 */
static void build_list(unsigned n, isd_t isd, const LIST_TYPE *columns,
                       LIST_TYPE *list) {
#if DUMER_P1 > 1
  LIST_TYPE *scratch0 = (LIST_TYPE *)(isd->scratch);
#if DUMER_P1 > 2
  LIST_TYPE *scratch1 =
      (LIST_TYPE *)((uint8_t *)scratch0 + AVX_PADDING(n * LIST_WIDTH) / 8);
#if DUMER_P1 > 3
  LIST_TYPE *scratch2 =
      (LIST_TYPE *)((uint8_t *)scratch1 + AVX_PADDING(n * LIST_WIDTH) / 8);
#if DUMER_P1 > 4
  LIST_TYPE *scratch3 =
      (LIST_TYPE *)((uint8_t *)scratch2 + AVX_PADDING(n * LIST_WIDTH) / 8);
#endif
#endif
#endif
#endif

  const LIST_TYPE *prev_scratch = columns;
#if DUMER_P1 == 1
  for (size_t i0 = 0; i0 < n; ++i0) {
#else /* DUMER_P1 > 1 */
  for (size_t i0 = n; i0-- > DUMER_P1 - 1;) {
#endif
    LIST_TYPE val = prev_scratch[i0];
#if DUMER_P1 >= 2
    xor_bcast(val, (uint8_t *)columns, (uint8_t *)scratch0,
              AVX_PADDING(i0 * LIST_WIDTH) / 256);

    LIST_TYPE *prev_scratch = scratch0;
#if DUMER_P1 == 2
    for (size_t i1 = 0; i1 < i0; ++i1) {
#else /* DUMER_P1 > 2 */
    for (size_t i1 = i0; i1-- > DUMER_P1 - 2;) {
#endif
      LIST_TYPE val = prev_scratch[i1];
#if DUMER_P1 >= 3
      xor_bcast(val, (uint8_t *)columns, (uint8_t *)scratch1,
                AVX_PADDING(i1 * LIST_WIDTH) / 256);

      LIST_TYPE *prev_scratch = scratch1;
#if DUMER_P1 == 3
      for (size_t i2 = 0; i2 < i1; ++i2) {
#else /* DUMER_P1 > 3 */
      for (size_t i2 = i1; i2-- > DUMER_P1 - 3;) {
#endif
        LIST_TYPE val = prev_scratch[i2];
#if DUMER_P1 >= 4
        xor_bcast(val, (uint8_t *)columns, (uint8_t *)scratch2,
                  AVX_PADDING(i2 * LIST_WIDTH) / 256);

        LIST_TYPE *prev_scratch = scratch2;
#if DUMER_P1 == 4
        for (size_t i3 = 0; i3 < i2; ++i3) {
#else /* DUMER_P1 > 4 */
        for (size_t i3 = i2; i3-- > DUMER_P1 - 4;) {
#endif
          LIST_TYPE val = prev_scratch[i3];
          xor_bcast(val, (uint8_t *)columns, (uint8_t *)scratch3,
                    AVX_PADDING(i3 * LIST_WIDTH) / 256);
#endif
#endif
#endif
          *(list++) = val;
#if DUMER_P1 >= 1
        }
#if DUMER_P1 >= 2
      }
#if DUMER_P1 >= 3
    }
#if DUMER_P1 >= 4
  }
#endif
#endif
#endif
#endif
}

/*
 * Build the list of the 'p' positions corresponding to the list built by the
 * preceding function.
 *
 * This list never changes so it is computed only once.
 */
static void build_list_pos(size_t n, uint16_t *pos) {
  /*
   * Backward iterating is better when we need to XOR because it avoids
   * computing an offset. But forward iterating is better when just copying
   * arrays.
   */
#if DUMER_P1 == 1
  for (size_t i0 = 0; i0 < n; ++i0) {
#else /* DUMER_P1 > 1 */
    for (size_t i0 = n; i0-- > DUMER_P1 - 1;) {
#endif
#if DUMER_P1 >= 2
#if DUMER_P1 == 2
    for (size_t i1 = 0; i1 < i0; ++i1) {
#else /* DUMER_P1 > 2 */
      for (size_t i1 = i0; i1-- > DUMER_P1 - 2;) {
#endif
#if DUMER_P1 >= 3
#if DUMER_P1 == 3
      for (size_t i2 = 0; i2 < i1; ++i2) {
#else /* DUMER_P1 > 3 */
        for (size_t i2 = i1; i2-- > DUMER_P1 - 3;) {
#endif
#if DUMER_P1 >= 4
#if DUMER_P1 == 4
        for (size_t i3 = 0; i3 < i2; ++i3) {
#else /* DUMER_P1 > 4 */
          for (size_t i3 = i2; i3-- > DUMER_P1 - 4;) {
#endif
          *(pos++) = i3;
#endif
          *(pos++) = i2;
#endif
          *(pos++) = i1;
#endif
          *(pos++) = i0;
        }
#if DUMER_P1 >= 2
      }
#if DUMER_P1 >= 3
    }
#if DUMER_P1 >= 4
  }
#endif
#endif
#endif
}

/*
 * To accelerate searching a value in a possibly huge list we build a lookup
 * table from the DUMER_LUT most significant bits of the elements of the list.
 *
 * Using this LUT, the binary search is done on a smaller range.
 */
static void build_lut(LIST_TYPE *list, size_t len_list, size_t *lut) {
  lut[0] = 0;
  lut[1 << DUMER_LUT] = len_list;
  size_t step = 1 << DUMER_LUT;
  size_t offset = 1 << (DUMER_LUT - 1);
  size_t nb = 1;
  for (size_t i = 0; i <= DUMER_LUT; ++i) {
    size_t idx = offset;
    for (size_t j = 0; j < nb; ++j) {
      lut[idx] =
          lut[idx - offset] + bin_search(list + lut[idx - offset],
                                         lut[idx + offset] - lut[idx - offset],
                                         idx << DUMER_LUT_SHIFT);
      idx += step;
    }

    step >>= 1;
    offset >>= 1;
    nb <<= 1;
  }
}

static void build_solution(size_t n, size_t r, size_t n1, shr_t shr, isd_t isd,
                           size_t pc, size_t idx1, size_t idx2, size_t shift) {
  size_t left = r - DUMER_L;
  for (size_t i = 0; i < n; ++i) {
    isd->solution[i] = 0;
  }
  for (size_t a = 0; a < DUMER_P1; ++a) {
    size_t column = shr->list1_pos[a + idx1 * DUMER_P1];
    size_t column_permuted = isd->perm[left + column];
    size_t column_shifted =
        column_permuted / r * r + (column_permuted + r - shift) % r;
    isd->solution[column_shifted] ^= 1;
  }
  for (size_t a = 0; a < DUMER_P2; ++a) {
    size_t column = shr->combinations2[a + idx2 * DUMER_P2] + n1 - DUMER_EPS;
    size_t column_permuted = isd->perm[left + column];
    size_t column_shifted =
        column_permuted / r * r + (column_permuted + r - shift) % r;
    isd->solution[column_shifted] ^= 1;
  }
  size_t pos_byte = 0;
  size_t pos_bit = 0;
  for (size_t column = 0; column < r; ++column) {
    if ((isd->test_syndrome[pos_byte] >> pos_bit) & 1) {
      size_t column_permuted = isd->perm[r - 1 - column];
      size_t column_shifted =
          column_permuted / r * r + (column_permuted + r - shift) % r;
      isd->solution[column_shifted] ^= 1;
    }
    pos_bit++;
    if (pos_bit == 64) {
      pos_bit = 0;
      ++pos_byte;
    }
  }
  isd->w_solution = pc;
}

void print_solution(size_t n, isd_t isd) {
#if DUMER_LW
  printf("%ld: ", isd->w_solution);
#endif
  for (size_t i = 0; i < n; ++i) {
    printf("%d", isd->solution[i]);
  }
  printf("\n");
  fflush(stdout);
}

static void xor_pairs(size_t r, size_t n2, isd_t isd) {
  size_t r_padded_bits = AVX_PADDING(r);
  size_t r_padded_qword = r_padded_bits / 64;
  size_t r_padded_ymm = r_padded_bits / 256;
  size_t xor_pairs_pos = 0;
  /* Compute the XORs of consecutive columns. */
  for (size_t i = 0; i < n2 + DUMER_EPS - 1; ++i) {
    xor_avx1((uint8_t *)&isd->columns2_full[i * r_padded_qword],
             (uint8_t *)&isd->columns2_full[(i + 1) * r_padded_qword],
             (uint8_t *)&isd->xor_pairs[xor_pairs_pos++ * r_padded_qword],
             r_padded_ymm);
  }
  /* Compute the XORs of columns distant by 2 positions. */
  for (size_t i = 0; i < n2 + DUMER_EPS - 2; ++i) {
    xor_avx1((uint8_t *)&isd->columns2_full[i * r_padded_qword],
             (uint8_t *)&isd->columns2_full[(i + 2) * r_padded_qword],
             (uint8_t *)&isd->xor_pairs[xor_pairs_pos++ * r_padded_qword],
             r_padded_ymm);
  }
}

static int find_collisions(size_t n, size_t r, size_t n1, shr_t shr,
                           isd_t isd) {
  int ret = 0;

  size_t r_padded_bits = AVX_PADDING(r);
  size_t r_padded_qword = r_padded_bits / 64;
  size_t r_padded_ymm = r_padded_bits / 256;

#if !(DUMER_DOOM) && !(DUMER_LW)
#if DUMER_P2 == 2
  uint16_t pos1 = shr->combinations2[0];
  uint16_t pos2 = shr->combinations2[1];
  xor_avx2((uint8_t *)isd->s_full,
           (uint8_t *)&isd->columns2_full[pos1 * r_padded_qword],
           (uint8_t *)&isd->columns2_full[pos2 * r_padded_qword],
           (uint8_t *)isd->current_syndrome, r_padded_ymm);
#elif DUMER_P2 == 3
    uint16_t pos1 = shr->combinations2[0];
    uint16_t pos2 = shr->combinations2[1];
    uint16_t pos3 = shr->combinations2[2];
    xor_avx3((uint8_t *)isd->s_full,
             (uint8_t *)&isd->columns2_full[pos1 * r_padded_qword],
             (uint8_t *)&isd->columns2_full[pos2 * r_padded_qword],
             (uint8_t *)&isd->columns2_full[pos3 * r_padded_qword],
             (uint8_t *)isd->current_syndrome, r_padded_ymm);
#elif DUMER_P2 == 4
  uint16_t pos1 = shr->combinations2[0];
  uint16_t pos2 = shr->combinations2[1];
  uint16_t pos3 = shr->combinations2[2];
  uint16_t pos4 = shr->combinations2[3];
  xor_avx4((uint8_t *)isd->s_full,
           (uint8_t *)&isd->columns2_full[pos1 * r_padded_qword],
           (uint8_t *)&isd->columns2_full[pos2 * r_padded_qword],
           (uint8_t *)&isd->columns2_full[pos3 * r_padded_qword],
           (uint8_t *)&isd->columns2_full[pos4 * r_padded_qword],
           (uint8_t *)isd->current_syndrome, r_padded_ymm);
#endif
#else  // DUMER_DOOM == 1 || DUMER_LW == 1
#if DUMER_P2 == 2
    uint16_t pos1 = shr->combinations2[0];
    uint16_t pos2 = shr->combinations2[1];
    xor_avx1((uint8_t *)&isd->columns2_full[pos1 * r_padded_qword],
             (uint8_t *)&isd->columns2_full[pos2 * r_padded_qword],
             (uint8_t *)isd->current_nosyndrome, r_padded_ymm);
#elif DUMER_P2 == 3
    uint16_t pos1 = shr->combinations2[0];
    uint16_t pos2 = shr->combinations2[1];
    uint16_t pos3 = shr->combinations2[2];
    xor_avx2((uint8_t *)&isd->columns2_full[pos1 * r_padded_qword],
             (uint8_t *)&isd->columns2_full[pos2 * r_padded_qword],
             (uint8_t *)&isd->columns2_full[pos3 * r_padded_qword],
             (uint8_t *)isd->current_nosyndrome, r_padded_ymm);
#elif DUMER_P2 == 4
    uint16_t pos1 = shr->combinations2[0];
    uint16_t pos2 = shr->combinations2[1];
    uint16_t pos3 = shr->combinations2[2];
    uint16_t pos4 = shr->combinations2[3];
    xor_avx3((uint8_t *)&isd->columns2_full[pos1 * r_padded_qword],
             (uint8_t *)&isd->columns2_full[pos2 * r_padded_qword],
             (uint8_t *)&isd->columns2_full[pos3 * r_padded_qword],
             (uint8_t *)&isd->columns2_full[pos4 * r_padded_qword],
             (uint8_t *)isd->current_nosyndrome, r_padded_ymm);
#endif
#endif

  uint64_t N = 0;
  goto first;

  for (; N < shr->nb_combinations2; ++N) {
#if DUMER_DOOM
    xor_avx1(
        (uint8_t *)isd->current_nosyndrome,
        (uint8_t *)&isd->xor_pairs[shr->combinations2_diff[N] * r_padded_qword],
        (uint8_t *)isd->current_nosyndrome, r_padded_ymm);
#else
      xor_avx1((uint8_t *)isd->current_syndrome,
               (uint8_t *)&isd
                   ->xor_pairs[shr->combinations2_diff[N] * r_padded_qword],
               (uint8_t *)isd->current_syndrome, r_padded_ymm);
#endif

#if !(DUMER_DOOM) || DUMER_LW
  first : {
    size_t shift = 0;
#else
    first:
      for (size_t shift = 0; shift < r; ++shift) {
        xor_avx1((uint8_t *)isd->current_nosyndrome,
                 (uint8_t *)&isd->s_full[shift * r_padded_qword],
                 (uint8_t *)isd->current_syndrome, r_padded_ymm);
#endif

    LIST_TYPE s_low = ((LIST_TYPE *)isd->current_syndrome)[0] & DUMER_L_MASK;

#if (DUMER_LUT) > 0
    size_t idx_lut = isd->list1_lut[s_low >> DUMER_LUT_SHIFT];
    size_t len_lut = isd->list1_lut[(s_low >> DUMER_LUT_SHIFT) + 1] - idx_lut;

    size_t idx_list =
        idx_lut + bin_search(isd->list1 + idx_lut, len_lut, s_low);
#elif (DUMER_LUT_SHIFT) == 0
        size_t idx_list = isd->list1_lut[s_low];
#else
    size_t idx_list = bin_search(isd->list1, shr->nb_combinations1, s_low);
#endif

    while (idx_list < shr->nb_combinations1 && isd->list1[idx_list] == s_low) {
      uint64_t idx_orig = isd->list1_idx[idx_list];

#if DUMER_P1 == 2
      uint16_t pos1 = shr->list1_pos[idx_orig * DUMER_P1];
      uint16_t pos2 = shr->list1_pos[1 + idx_orig * DUMER_P1];
      xor_avx2((uint8_t *)isd->current_syndrome,
               (uint8_t *)&isd->columns1_full[pos1 * r_padded_qword],
               (uint8_t *)&isd->columns1_full[pos2 * r_padded_qword],
               (uint8_t *)isd->test_syndrome, r_padded_ymm);
#elif DUMER_P1 == 3
          uint16_t pos1 = shr->list1_pos[idx_orig * DUMER_P1];
          uint16_t pos2 = shr->list1_pos[1 + idx_orig * DUMER_P1];
          uint16_t pos3 = shr->list1_pos[2 + idx_orig * DUMER_P1];
          xor_avx3((uint8_t *)isd->current_syndrome,
                   (uint8_t *)&isd->columns1_full[pos1 * r_padded_qword],
                   (uint8_t *)&isd->columns1_full[pos2 * r_padded_qword],
                   (uint8_t *)&isd->columns1_full[pos3 * r_padded_qword],
                   (uint8_t *)isd->test_syndrome, r_padded_ymm);
#elif DUMER_P1 == 4
      uint16_t pos1 = shr->list1_pos[idx_orig * DUMER_P1];
      uint16_t pos2 = shr->list1_pos[1 + idx_orig * DUMER_P1];
      uint16_t pos3 = shr->list1_pos[2 + idx_orig * DUMER_P1];
      uint16_t pos4 = shr->list1_pos[3 + idx_orig * DUMER_P1];
      xor_avx4((uint8_t *)isd->current_syndrome,
               (uint8_t *)&isd->columns1_full[pos1 * r_padded_qword],
               (uint8_t *)&isd->columns1_full[pos2 * r_padded_qword],
               (uint8_t *)&isd->columns1_full[pos3 * r_padded_qword],
               (uint8_t *)&isd->columns1_full[pos4 * r_padded_qword],
               (uint8_t *)isd->test_syndrome, r_padded_ymm);
#endif
      size_t pc = popcount(isd->test_syndrome, r_padded_qword, isd->w_target);
      /* Fusion error patterns from both lists. */
      if (pc <= isd->w_target) {
        size_t a1 = 0;
        size_t a2 = 0;
        size_t column1 = shr->list1_pos[idx_orig * DUMER_P1];
        size_t column2 = shr->combinations2[N * DUMER_P2] + n1 - DUMER_EPS;
        while (a2 < DUMER_P2 && a1 < DUMER_P1) {
          if (column1 < column2) {
            ++pc;
            ++a1;
            column1 = shr->list1_pos[a1 + idx_orig * DUMER_P1];
          } else if (column1 > column2) {
            ++pc;
            ++a2;
            column2 = shr->combinations2[a2 + N * DUMER_P2] + n1 - DUMER_EPS;
          } else {
            ++a1;
            ++a2;
            column1 = shr->list1_pos[a1 + idx_orig * DUMER_P1];
            column2 = shr->combinations2[a2 + N * DUMER_P2] + n1 - DUMER_EPS;
          }
        }
        pc += DUMER_P2 + DUMER_P1 - a1 - a2;
      }

      if (pc > 0 && pc <= isd->w_target) {
#if DUMER_LW
        omp_set_lock(&shr->w_best_lock);
        if (pc >= shr->w_best) {
          isd->w_target = shr->w_best - 1;
          omp_unset_lock(&shr->w_best_lock);
          continue;
        } else {
          shr->w_best = pc;
          omp_unset_lock(&shr->w_best_lock);
          isd->w_target = pc - 1;
        }
#endif
        /* Found it! */
        isd->w_solution = pc;
        build_solution(n, r, n1, shr, isd, pc, idx_orig, N, shift);
        ret = 1;
#if !(DUMER_LW) && !(BENCHMARK)
        return ret;
#endif
      }

      ++idx_list;
    }
  }
  }
  return ret;
}

shr_t alloc_shr(size_t n1, size_t n2) {
  shr_t shr = malloc(sizeof(struct shared));

  shr->nb_combinations1 = bincoef(n1 + DUMER_EPS, DUMER_P1);
  shr->nb_combinations2 = bincoef(n2 + DUMER_EPS, DUMER_P2);

  shr->list1_pos = malloc(DUMER_P1 * shr->nb_combinations1 * sizeof(uint16_t));

  shr->combinations2 =
      malloc(shr->nb_combinations2 * DUMER_P2 * sizeof(uint16_t));
  shr->combinations2_diff = malloc(shr->nb_combinations2 * sizeof(uint16_t));

  if (!shr->list1_pos || !shr->combinations2 || !shr->combinations2_diff)
    return NULL;

#if DUMER_LW
  omp_init_lock(&shr->w_best_lock);
  shr->w_best = INT_MAX;
#endif

  matrix_alloc_gray_code(&shr->gray_rev, &shr->gray_diff);

  return shr;
}

void free_shr(shr_t shr) {
  free(shr->list1_pos);
  free(shr->combinations2);
  free(shr->combinations2_diff);
#if DUMER_LW
  omp_destroy_lock(&shr->w_best_lock);
#endif

  matrix_free_gray_code(shr->gray_rev, shr->gray_diff);
  free(shr);
}

void init_shr(shr_t shr, size_t n, size_t k, size_t n1, size_t n2) {
  /*
   * Precompute Chase's sequence.
   *
   * With this sequence, computing all the combinations of P2 elements among n2
   * elements only requires one XOR per new combination by using (2 * N2 - 3)
   * precomputed XORed pairs of columns.
   */
  chase(n2 + DUMER_EPS, DUMER_P2, shr->combinations2, shr->combinations2_diff);

  build_list_pos(n1 + DUMER_EPS, shr->list1_pos);

  matrix_build_gray_code(shr->gray_rev, shr->gray_diff);
  shr->k_opt = matrix_opt_k(n - k, n);
}

isd_t alloc_isd(size_t n, size_t k, size_t r, size_t n1, size_t n2,
                uint64_t nb_combinations1, size_t k_opt) {
  isd_t isd = malloc(sizeof(struct isd));

  /* We make sure that the 32 rows before isd->A are allocated so that we do
   * not have to deal with edge cases during transposition. */
#if DUMER_LW
  (void)(k);
  isd->A = matrix_alloc(r + 64, n);
  matrix_reset(isd->A, r + 64, n);
#elif !(DUMER_DOOM)  // && !(DUMER_LW)
  (void)(k);
  isd->A = matrix_alloc(r + 64, n + 1);
  matrix_reset(isd->A, r + 64, n + 1);
#else                // DUMER_DOOM && !(DUMER_LW)
  isd->A = matrix_alloc(r + 64, n + k);
  matrix_reset(isd->A, r + 64, n + k);
#endif
  isd->At = matrix_alloc(n + k + 64, r + 64);
  matrix_reset(isd->At, n + k + 64, r + 64);
  if (!isd->A || !isd->At) return NULL;
  isd->A = isd->A + 32;
  isd->At = isd->At + 32;

  isd->perm = malloc(n * sizeof(size_t));
  if (!isd->perm) return NULL;

  isd->size_list1 = LIST_WIDTH * nb_combinations1;
  isd->list1 = malloc(isd->size_list1 / 8);
  isd->list1_aux = malloc(isd->size_list1 / 8);
  isd->list1_idx = malloc(nb_combinations1 * sizeof(size_t));
  isd->list1_aux2 = malloc(nb_combinations1 * sizeof(size_t));
  isd->list1_lut = malloc(((1 << DUMER_LUT) + 1) * sizeof(size_t));
  if (!isd->list1 || !isd->list1_aux || !isd->list1_idx || !isd->list1_aux2 ||
      !isd->list1_lut)
    return NULL;

  isd->size_columns1_low = AVX_PADDING(LIST_WIDTH * (n1 + DUMER_EPS));
  isd->columns1_low = aligned_alloc(32, isd->size_columns1_low / 8);
  if (!isd->size_columns1_low || !isd->columns1_low) return NULL;

  isd->size_columns1_full = AVX_PADDING(r) * (n1 + DUMER_EPS);
  isd->size_columns2_full = AVX_PADDING(r) * (n2 + DUMER_EPS);
  isd->columns1_full = aligned_alloc(32, isd->size_columns1_full / 8);
  isd->columns2_full = aligned_alloc(32, isd->size_columns2_full / 8);
  if (!isd->columns1_full || !isd->columns2_full) return NULL;

#if !(DUMER_LW) && !(DUMER_DOOM)
  isd->s_full = aligned_alloc(32, AVX_PADDING(r) / 8);
  if (!isd->s_full) return NULL;
#elif !(DUMER_LW) && DUMER_DOOM
  isd->s_full = aligned_alloc(32, k * AVX_PADDING(r) / 8);
  if (!isd->s_full) return NULL;
#endif

  size_t r_padded_bits = AVX_PADDING(r);
  size_t r_padded_qword = r_padded_bits / 64;
  isd->test_syndrome = aligned_alloc(32, r_padded_qword * sizeof(uint64_t));
#if DUMER_DOOM || DUMER_LW
  isd->current_nosyndrome =
      aligned_alloc(32, r_padded_qword * sizeof(uint64_t));
  if (!isd->current_nosyndrome) return NULL;
#endif
#if DUMER_LW
  isd->current_syndrome = isd->current_nosyndrome;
#else
  isd->current_syndrome =
      aligned_alloc(32, r_padded_qword * sizeof(uint64_t));
#endif
  isd->xor_pairs = aligned_alloc(
      32, (2 * (n2 + DUMER_EPS) - 3) * r_padded_qword * sizeof(uint64_t));

  isd->xor_rows =
      aligned_alloc(32, (1L << k_opt) * AVX_PADDING(n) / 64 * sizeof(uint64_t));

  isd->scratch = aligned_alloc(
      32, DUMER_P1 * AVX_PADDING((n1 + DUMER_EPS) * LIST_WIDTH) / 8);
  if (!isd->test_syndrome || !isd->current_syndrome || !isd->xor_pairs ||
      !isd->xor_rows || !isd->scratch)
    return NULL;

  return isd;
}

void free_isd(isd_t isd, size_t r, size_t n) {
  matrix_free(isd->A - 32, r);
  matrix_free(isd->At - 32, n);
  free(isd->perm);

  free(isd->list1);
  free(isd->list1_aux);
  free(isd->list1_idx);
  free(isd->list1_aux2);
  free(isd->list1_lut);

  free(isd->columns1_low);

  free(isd->columns1_full);

  free(isd->columns2_full);

#if !(DUMER_LW)
  free(isd->s_full);
#endif

  free(isd->solution);

  free(isd->scratch);

  free(isd->test_syndrome);
#if DUMER_DOOM || DUMER_LW
  free(isd->current_nosyndrome);
#endif
#if !(DUMER_LW)
  free(isd->current_syndrome);
#endif
  free(isd->xor_pairs);
  free(isd->xor_rows);

  free(isd);
}

void init_isd(isd_t isd, enum type current_type, size_t n, size_t k, size_t w,
              uint8_t *mat_h, uint8_t *mat_s) {
#if DUMER_LW
  (void)w;
  (void)mat_s;
#endif
  if (!seed_random(&isd->S0, &isd->S1)) exit(EXIT_FAILURE);

  /* Build the M4RI matrix. */
  for (size_t i = 0; i < n - k; ++i) {
    isd->A[i][i / WORD_SIZE] |= 1UL << (i % WORD_SIZE);
  }
  if (current_type == QC) {
    for (size_t j = 0; j < k; ++j) {
      for (size_t i = 0; i < n - k; ++i) {
        if (mat_h[(i - j + k) % k])
          isd->A[i][(k + j) / WORD_SIZE] |= 1L << ((k + j) % WORD_SIZE);
      }
    }
  } else if (current_type == SD || current_type == LW || current_type == GO) {
    for (size_t j = 0; j < k; ++j) {
      for (size_t i = 0; i < n - k; ++i) {
        if (mat_h[i + (n - k) * j])
          isd->A[i][(n - k + j) / WORD_SIZE] |= 1UL
                                                << ((n - k + j) % WORD_SIZE);
      }
    }
  }
  /* Matrix A is extended with the syndrome(s). */
#if !(DUMER_LW) && !(DUMER_DOOM)
  for (size_t i = 0; i < n - k; ++i) {
    if (mat_s[i]) isd->A[i][n / WORD_SIZE] |= 1L << (n % WORD_SIZE);
  }
#elif !(DUMER_LW) && DUMER_DOOM
    /*
     * In quasi-cyclic codes, a circular permutation of a syndrome is the
     * syndrome of the blockwise circularly permuted error pattern.
     */
    for (size_t j = 0; j < k; ++j) {
      for (size_t i = 0; i < n - k; ++i) {
        if (mat_s[(i - j + k) % k])
          isd->A[i][(n + j) / WORD_SIZE] |= 1L << ((n + j) % WORD_SIZE);
      }
    }
#endif
#if DUMER_LW
  matrix_transpose_rev_cols(isd->At, isd->A, n - k, n);
#elif !(DUMER_DOOM)  // && !(DUMER_LW)
  matrix_transpose_rev_cols(isd->At, isd->A, n - k, n + 1);
#else                // DUMER_DOOM && !(DUMER_LW)
  matrix_transpose_rev_cols(isd->At, isd->A, n - k, n + k);
#endif

  for (size_t i = 0; i < n; ++i) {
    isd->perm[i] = i;
  }

  isd->solution = malloc(n * sizeof(uint8_t));
#if DUMER_LW
  isd->w_target = n;
#else
    isd->w_target = w;
#endif
}

int dumer(size_t n, size_t k, size_t r, size_t n1, size_t n2, shr_t shr,
          isd_t isd) {
  /* Choose a random information set and do a Gaussian elimination. */
  choose_is(isd->A, isd->At, isd->perm, n, k, DUMER_L, shr->k_opt,
            shr->gray_rev, shr->gray_diff, isd->xor_rows, &isd->S0, &isd->S1);

#if DUMER_LW
  matrix_transpose_rev_cols(isd->At, isd->A, r, n);
#elif !(DUMER_DOOM)  // && !(DUMER_LW)
  matrix_transpose_rev_cols(isd->At, isd->A, r, n + 1);
#else                // DUMER_DOOM && !(DUMER_LW)
  matrix_transpose_rev_cols(isd->At, isd->A, r, n + k);
#endif

  get_columns_H_prime(isd->At + r - DUMER_L, isd->columns1_low, n1 + DUMER_EPS);

  /*
   * For the first list, we only keep the LIST_WIDTH least significant bits.
   *
   * The full column is then fully computed when there is a collision on the
   * LIST_WIDTH least significant bits in list1 and in list2.
   */
  build_list(n1 + DUMER_EPS, isd, isd->columns1_low, isd->list1);

  /* Keep the original index of an element of the list when sorting. */
  for (uint64_t i = 0; i < shr->nb_combinations1; ++i) {
    isd->list1_idx[i] = i;
  }
  sort(isd->list1, isd->list1_idx, isd->list1_aux, isd->list1_aux2,
       shr->nb_combinations1);
#if (DUMER_LUT) > 0
  /* The lookup table speeds up searching in the sorted list. */
  build_lut(isd->list1, shr->nb_combinations1, isd->list1_lut);
#endif

  get_columns_H_prime_avx(isd->At + r - DUMER_L, isd->columns1_full,
                          n1 + DUMER_EPS, r);
  get_columns_H_prime_avx(isd->At + r - DUMER_L + n1 - DUMER_EPS,
                          isd->columns2_full, n2 + DUMER_EPS, r);

#if !(DUMER_LW) && !(DUMER_DOOM)
  get_columns_H_prime_avx(isd->At + n, isd->s_full, 1, r);
#elif !(DUMER_LW) && DUMER_DOOM
  get_columns_H_prime_avx(isd->At + n, isd->s_full, r, r);
#endif
  xor_pairs(r, n2, isd);

  /*
   * As there is always at least one element matching the LIST_WIDTH least
   * significant bits of an element of list2 in list1, the elements of list2
   * are computed fully from the beginning.
   *
   * Using Chase's sequence, list2 is computed doing only one XOR per
   * element.
   */
  return find_collisions(n, r, n1, shr, isd);
}
