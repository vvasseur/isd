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
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "bits.h"
#include "dumer.h"
#include "xoroshiro128plus.h"

#include "sort.h"

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

/* Apply the same permutation to an M4RI matrix and an array. */
static void fisher_yates_m4ri(mzd_t *A, size_t *perm, rci_t n, size_t n_stop,
                              uint64_t *S0, uint64_t *S1) {
  for (size_t i = 0; i < n_stop; ++i) {
    uint32_t rand = i + random_lim(n - 1 - i, S0, S1);
    mzd_col_swap(A, i, rand);
    size_t swp = perm[i];
    perm[i] = perm[rand];
    perm[rand] = swp;
  }
}

/* Randomly choose an information set and perform a Gaussian elimination. */
static void choose_is(mzd_t *A, size_t *perm, size_t n, size_t k, size_t l,
                      uint64_t *S0, uint64_t *S1) {
  /*
   * Pick a permutation and perform Gaussian elimination.
   *
   * M4RI tries to do a full Gaussian elimination. But as long as the rank is
   * above or equal to n - k - l the permutation is usable.
   */
  size_t r = 0;
  while (r < n - k - l) {
    fisher_yates_m4ri(A, perm, n, n - k - l, S0, S1);
    r = mzd_echelonize_m4ri(A, 1, 0);
  }
}

/* Extract columns from *A and keep data 32-byte aligned (fitting AVX
 * registers). */
static void get_columns_H_prime_avx(mzd_t *A, uint64_t *columns, size_t n,
                                    size_t r, size_t off) {
  size_t pos_byte = 0;
  size_t pos_bit = 0;
  columns[pos_byte] = 0;
  for (size_t j = 0; j < n; ++j) {
    for (size_t i = 0; i < r; ++i) {
      columns[pos_byte] |=
          (mzd_read_bit(A, r - 1 - i, j + off) ? (1L << pos_bit) : 0);
      if (++pos_bit == 64) {
        pos_bit = 0;
        columns[++pos_byte] = 0;
      }
    }
    if (pos_bit != 0) {
      pos_bit = 0;
      columns[++pos_byte] = 0;
    }
    size_t pad = (4 - pos_byte % 4) % 4;
    while (pad--) {
      columns[++pos_byte] = 0;
    }
  }
}

/* Extract columns from *A and keep data LIST_WIDTH-byte aligned. */
static void get_columns_H_prime(mzd_t *A, LIST_TYPE *columns, size_t n,
                                size_t r, size_t l, size_t off) {
  size_t pos_word = 0;
  size_t pos_bit = 0;
  columns[pos_word] = 0;
  for (size_t j = 0; j < n; ++j) {
    for (size_t i = 0; i < l; ++i) {
      columns[pos_word] |=
          (mzd_read_bit(A, r - 1 - i, j + off) ? (1L << pos_bit) : 0);
      if (++pos_bit == LIST_WIDTH) {
        pos_bit = 0;
        columns[++pos_word] = 0;
      }
    }
    if (pos_bit != 0) {
      pos_bit = 0;
      columns[++pos_word] = 0;
    }
  }
}

/* ceil(log2(x)) */
static unsigned lb(unsigned long x) {
  if (x <= 1) return 0;
  return (8 * sizeof(unsigned long)) - __builtin_clzl(x - 1);
}

/*
 * See Paul Khuong's
 * https://www.pvk.ca/Blog/2012/07/03/binary-search-star-eliminates-star-branch-mispredictions/
 */
static size_t bin_search(const LIST_TYPE *list, size_t len_list,
                         LIST_TYPE value) {
  if (len_list <= 1) return 0;
  unsigned log = lb(len_list) - 1;
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
#if DUMER_L <= 8
    xor_bcast_32(((uint32_t)val << 24) | ((uint32_t)val << 16) |
                     ((uint32_t)val << 8) | val,
                 (uint8_t *)columns, (uint8_t *)scratch0,
                 AVX_PADDING(i0 * 8) / 256);
#elif DUMER_L <= 16
    xor_bcast_32(((uint32_t)val << 16) | val, (uint8_t *)columns,
                 (uint8_t *)scratch0, AVX_PADDING(i0 * 16) / 256);
#elif DUMER_L <= 32
    xor_bcast_32(val, (uint8_t *)columns, (uint8_t *)scratch0,
                 AVX_PADDING(i0 * 32) / 256);
#elif DUMER_L <= 64
    xor_bcast_64(val, (uint8_t *)columns, (uint8_t *)scratch0,
                 AVX_PADDING(i0 * 64) / 256);
#endif

    LIST_TYPE *prev_scratch = scratch0;
#if DUMER_P1 == 2
    for (size_t i1 = 0; i1 < i0; ++i1) {
#else /* DUMER_P1 > 2 */
    for (size_t i1 = i0; i1-- > DUMER_P1 - 2;) {
#endif
      LIST_TYPE val = prev_scratch[i1];
#if DUMER_P1 >= 3
#if DUMER_L <= 8
      xor_bcast_32(((uint32_t)val << 24) | ((uint32_t)val << 16) |
                       ((uint32_t)val << 8) | val,
                   (uint8_t *)columns, (uint8_t *)scratch1,
                   AVX_PADDING(i1 * 8) / 256);
#elif DUMER_L <= 16
      xor_bcast_32(((uint32_t)val << 16) | val, (uint8_t *)columns,
                   (uint8_t *)scratch1, AVX_PADDING(i1 * 16) / 256);
#elif DUMER_L <= 32
      xor_bcast_32(val, (uint8_t *)columns, (uint8_t *)scratch1,
                   AVX_PADDING(i1 * 32) / 256);
#elif DUMER_L <= 64
      xor_bcast_64(val, (uint8_t *)columns, (uint8_t *)scratch1,
                   AVX_PADDING(i1 * 64) / 256);
#endif

      LIST_TYPE *prev_scratch = scratch1;
#if DUMER_P1 == 3
      for (size_t i2 = 0; i2 < i1; ++i2) {
#else /* DUMER_P1 > 3 */
      for (size_t i2 = i1; i2-- > DUMER_P1 - 3;) {
#endif
        LIST_TYPE val = prev_scratch[i2];
#if DUMER_P1 >= 4
#if DUMER_L <= 8
        xor_bcast_32(((uint32_t)val << 24) | ((uint32_t)val << 16) |
                         ((uint32_t)val << 8) | val,
                     (uint8_t *)columns, (uint8_t *)scratch2,
                     AVX_PADDING(i2 * 8) / 256);
#elif DUMER_L <= 16
        xor_bcast_32(((uint32_t)val << 16) | val, (uint8_t *)columns,
                     (uint8_t *)scratch2, AVX_PADDING(i2 * 16) / 256);
#elif DUMER_L <= 32
        xor_bcast_32(val, (uint8_t *)columns, (uint8_t *)scratch2,
                     AVX_PADDING(i2 * 32) / 256);
#elif DUMER_L <= 64
        xor_bcast_64(val, (uint8_t *)columns, (uint8_t *)scratch2,
                     AVX_PADDING(i2 * 64) / 256);
#endif

        LIST_TYPE *prev_scratch = scratch2;
#if DUMER_P1 == 4
        for (size_t i3 = 0; i3 < i2; ++i3) {
#else /* DUMER_P1 > 4 */
        for (size_t i3 = i2; i3-- > DUMER_P1 - 4;) {
#endif
          LIST_TYPE val = prev_scratch[i3];
#if DUMER_L <= 8
          xor_bcast_32(((uint32_t)val << 24) | ((uint32_t)val << 16) |
                           ((uint32_t)val << 8) | val,
                       (uint8_t *)columns, (uint8_t *)scratch3,
                       AVX_PADDING(i3 * 8) / 256);
#elif DUMER_L <= 16
          xor_bcast_32(((uint32_t)val << 16) | val, (uint8_t *)columns,
                       (uint8_t *)scratch3, AVX_PADDING(i3 * 16) / 256);
#elif DUMER_L <= 32
        xor_bcast_32(val, (uint8_t *)columns, (uint8_t *)scratch3,
                     AVX_PADDING(i3 * 32) / 256);
#elif DUMER_L <= 64
        xor_bcast_64(val, (uint8_t *)columns, (uint8_t *)scratch3,
                     AVX_PADDING(i3 * 64) / 256);
#endif
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

static int find_collisions(size_t n, size_t r, size_t n1, size_t n2, shr_t shr,
                           isd_t isd) {
  int ret = 0;
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
#if DUMER_L == 64
  uint64_t mask = ~0UL;
#else
    uint64_t mask = (uint64_t)((1UL << DUMER_L) - 1);
#endif

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

    LIST_TYPE s_low = ((LIST_TYPE *)isd->current_syndrome)[0] & mask;

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
#if !(DUMER_LW)
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

  return shr;
}

void free_shr(shr_t shr) {
  free(shr->list1_pos);
  free(shr->combinations2);
  free(shr->combinations2_diff);
#if DUMER_LW
  omp_destroy_lock(&shr->w_best_lock);
#endif
  free(shr);
}

void init_shr(shr_t shr, size_t n1, size_t n2) {
  /*
   * Precompute Chase's sequence.
   *
   * With this sequence, computing all the combinations of P2 elements among n2
   * elements only requires one XOR per new combination by using (2 * N2 - 3)
   * precomputed XORed pairs of columns.
   */
  chase(n2 + DUMER_EPS, DUMER_P2, shr->combinations2, shr->combinations2_diff);

  build_list_pos(n1 + DUMER_EPS, shr->list1_pos);
}

isd_t alloc_isd(size_t n, size_t k, size_t r, size_t n1, size_t n2,
                uint64_t nb_combinations1) {
  isd_t isd = malloc(sizeof(struct isd));

#if DUMER_LW
  isd->A = mzd_init(r, n);
#elif !(DUMER_DOOM)  // && !(DUMER_LW)
    isd->A = mzd_init(r, n + 1);
#else                // DUMER_DOOM && !(DUMER_LW)
  isd->A = mzd_init(r, n + k);
#endif
  if (!isd->A) return NULL;

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

  isd->scratch = aligned_alloc(
      32, DUMER_P1 * AVX_PADDING((n1 + DUMER_EPS) * LIST_WIDTH) / 8);
  if (!isd->test_syndrome || !isd->current_syndrome || !isd->xor_pairs ||
      !isd->scratch)
    return NULL;

  return isd;
}

void free_isd(isd_t isd) {
  mzd_free(isd->A);
  free(isd->perm);

  free(isd->list1);
  free(isd->list1_aux);
  free(isd->list1_idx);
  free(isd->list1_aux2);

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
    mzd_write_bit(isd->A, i, i, 1);
  }
  if (current_type == QC) {
    for (size_t j = 0; j < k; ++j) {
      for (size_t i = 0; i < n - k; ++i) {
        mzd_write_bit(isd->A, i, k + j, mat_h[(i - j + k) % k]);
      }
    }
  } else if (current_type == SD || current_type == LW || current_type == GO) {
    for (size_t j = 0; j < k; ++j) {
      for (size_t i = 0; i < n - k; ++i) {
        mzd_write_bit(isd->A, i, n - k + j, mat_h[i + (n - k) * j]);
      }
    }
  }

  /* Matrix A is extended with the syndrome(s). */
#if !(DUMER_LW) && !(DUMER_DOOM)
  for (size_t i = 0; i < n - k; ++i) {
    mzd_write_bit(isd->A, i, n, mat_s[i]);
  }
#elif !(DUMER_LW) && DUMER_DOOM
    /*
     * In quasi-cyclic codes, a circular permutation of a syndrome is the
     * syndrome of the blockwise circularly permuted error pattern.
     */
    for (size_t j = 0; j < k; ++j) {
      for (size_t i = 0; i < n - k; ++i) {
        mzd_write_bit(isd->A, i, n + j, mat_s[(i - j + k) % k]);
      }
    }
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

size_t dumer(size_t n, size_t k, size_t r, size_t n1, size_t n2, shr_t shr,
             isd_t isd) {
  /* Choose a random information set and do a Gaussian elimination. */
  choose_is(isd->A, isd->perm, n, k, DUMER_L, &isd->S0, &isd->S1);

  get_columns_H_prime(isd->A, isd->columns1_low, n1 + DUMER_EPS, r, DUMER_L,
                      r - DUMER_L);

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

  get_columns_H_prime_avx(isd->A, isd->columns1_full, n1 + DUMER_EPS, r,
                          r - DUMER_L);
  get_columns_H_prime_avx(isd->A, isd->columns2_full, n2 + DUMER_EPS, r,
                          r - DUMER_L + n1 - DUMER_EPS);
#if !(DUMER_LW) && !(DUMER_DOOM)
  get_columns_H_prime_avx(isd->A, isd->s_full, 1, r, n);
#elif !(DUMER_LW) && DUMER_DOOM
    get_columns_H_prime_avx(isd->A, isd->s_full, r, r, n);
#endif

  /*
   * As there is always at least one element matching the LIST_WIDTH least
   * significant bits of an element of list2 in list1, the elements of list2
   * are computed fully from the beginning.
   *
   * Using Chase's sequence, list2 is computed doing only one XOR per
   * element.
   */
  return find_collisions(n, r, n1, n2, shr, isd);
}
