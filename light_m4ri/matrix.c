/******************************************************************************
 *
 *                 M4RI: Linear Algebra over GF(2)
 *
 *    Copyright (C) 2007 Gregory Bard <gregory.bard@ieee.org>
 *    Copyright (C) 2007,2008 Martin Albrecht <malb@informatik.uni-bremen.de>
 *    Copyright (C) 2020 Valentin Vasseur
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
 *  version 2 or higher.
 *
 *    This code is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    General Public License for more details.
 *
 *  The full text of the GPL is available at:
 *
 *                  http://www.gnu.org/licenses/
 ******************************************************************************/
#include "matrix.h"

#include <stdlib.h>

#include "../bits.h"

matrix_t matrix_alloc(size_t rows, size_t cols) {
  matrix_t ret = malloc(rows * sizeof(word_t *));
  size_t cols_padded_avx = AVX_PADDING(cols);
  for (size_t i = 0; i < rows; ++i) {
    ret[i] = aligned_alloc(32, cols_padded_avx * sizeof(word_t) / WORD_SIZE);
  }
  return ret;
}

void matrix_reset(matrix_t M, size_t rows, size_t cols) {
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < (cols + WORD_SIZE - 1) / WORD_SIZE; ++j) {
      M[i][j] = 0;
    }
  }
}

void matrix_free(matrix_t M, size_t rows) {
  for (size_t i = 0; i < rows; ++i) {
    free(M[i]);
  }
  free(M);
}

void matrix_swap_rows(matrix_t M, size_t i, size_t j) {
  word_t *tmp = M[i];
  M[i] = M[j];
  M[j] = tmp;
}

void matrix_swap_cols(matrix_t M, size_t i, size_t j, size_t rows) {
  if (i == j) return;

  size_t iword = i / WORD_SIZE;
  size_t jword = j / WORD_SIZE;

  size_t ibit = i % WORD_SIZE;
  size_t jbit = j % WORD_SIZE;

  if (iword == jword) {
    size_t max_bit;
    size_t min_bit;
    if (ibit < jbit) {
      max_bit = jbit;
      min_bit = ibit;
    } else {
      max_bit = ibit;
      min_bit = jbit;
    }
    size_t offset = max_bit - min_bit;
    word_t mask = (word_t)1 << min_bit;
    for (size_t i = 0; i < rows; ++i) {
      word_t *row = M[i] + iword;
      word_t xor = (*row ^ (*row >> offset)) & mask;
      *row ^= xor;
      *row ^= xor << offset;
    }
  } else {
    size_t max_word;
    size_t min_word;
    size_t max_bit;
    size_t min_bit;
    if (ibit < jbit) {
      max_bit = jbit;
      min_bit = ibit;
      max_word = jword;
      min_word = iword;
    } else {
      max_bit = ibit;
      min_bit = jbit;
      max_word = iword;
      min_word = jword;
    }
    size_t offset = max_bit - min_bit;
    word_t mask = (word_t)1 << min_bit;
    for (size_t i = 0; i < rows; ++i) {
      word_t *min_row = M[i] + min_word;
      word_t *max_row = M[i] + max_word;
      word_t xor = (*min_row ^ (*max_row >> offset)) & mask;
      *min_row ^= xor;
      *max_row ^= xor << offset;
    }
  }
}

static int gray(int i, int k) {
  int lastbit = 0;
  int res = 0;
  for (int j = k; j-- > 0;) {
    int bit = i & (1 << j);
    res |= (lastbit >> 1) ^ bit;
    lastbit = bit;
  }
  return res;
}

void matrix_alloc_gray_code(int ***rev, int ***diff) {
  *rev = malloc((MAX_K + 1) * sizeof(int *));
  *diff = malloc((MAX_K + 1) * sizeof(int *));

  for (size_t k = 0; k <= MAX_K; ++k) {
    (*rev)[k] = malloc((1 << k) * sizeof(int));
    (*diff)[k] = malloc((1 << k) * sizeof(int));
  }
}

void matrix_free_gray_code(int **rev, int **diff) {
  for (size_t k = 0; k <= MAX_K; ++k) {
    free(rev[k]);
    free(diff[k]);
  }
  free(rev);
  free(diff);
}

void matrix_build_gray_code(int **rev, int **diff) {
  for (size_t k = 0; k <= MAX_K; ++k) {
    for (size_t i = 0; i < 1UL << k; ++i) {
      rev[k][gray(i, k)] = i;
    }

    for (size_t i = k + 1; i-- > 0;) {
      for (size_t j = 1; j < (1UL << i) + 1; ++j) {
        diff[k][j * (1 << (k - i)) - 1] = k - i;
      }
    }
  }
}

size_t matrix_gauss_submatrix(matrix_t M, size_t r, size_t c, size_t rows,
                              size_t cols, size_t k) {
  size_t start_row = r;
  size_t j;
  size_t cols_padded_ymm = AVX_PADDING(cols) / 256;
  for (j = c; j < c + k; ++j) {
    int found = 0;
    for (size_t i = start_row; i < rows; ++i) {
      for (size_t l = 0; l < j - c; ++l)
        if ((M[i][(c + l) / WORD_SIZE] >> ((c + l) % WORD_SIZE)) & 1)
          xor_avx1((uint8_t *)M[r + l], (uint8_t *)M[i], (uint8_t *)M[i],
                   cols_padded_ymm);

      if ((M[i][j / WORD_SIZE] >> (j % WORD_SIZE)) & 1) {
        matrix_swap_rows(M, i, start_row);
        for (size_t l = r; l < start_row; ++l) {
          if ((M[l][j / WORD_SIZE] >> (j % WORD_SIZE)) & 1)
            xor_avx1((uint8_t *)M[start_row], (uint8_t *)M[l], (uint8_t *)M[l],
                     cols_padded_ymm);
        }
        ++start_row;
        found = 1;
        break;
      }
    }
    if (found == 0) {
      break;
    }
  }
  return j - c;
}

void matrix_make_table(matrix_t M, size_t r, size_t cols, size_t k, uint64_t *T,
                       int **diff) {
  size_t cols_padded = AVX_PADDING(cols);
  size_t cols_padded_word = cols_padded / 64;
  size_t cols_padded_ymm = cols_padded / 256;

  for (size_t i = 0; i < cols_padded_word; ++i) {
    T[i] = 0L;
  }

  for (size_t i = 0; i + 1 < 1UL << k; ++i) {
    xor_avx1((uint8_t *)M[r + diff[k][i]], (uint8_t *)T,
             (uint8_t *)(T + cols_padded_word), cols_padded_ymm);
    T += cols_padded_word;
  }
}

word_t matrix_read_bits(matrix_t M, size_t x, size_t y, size_t n) {
  int spot = y % WORD_SIZE;
  word_t block = y / WORD_SIZE;
  int spill = spot + n - WORD_SIZE;
  word_t temp = (spill <= 0) ? M[x][block] << -spill
                             : (M[x][block + 1] << (WORD_SIZE - spill)) |
                                   (M[x][block] >> spill);
  return temp >> (WORD_SIZE - n);
}

void matrix_process_rows(matrix_t M, size_t rstart, size_t cstart, size_t rstop,
                         size_t k, size_t cols, uint64_t *T, int **rev) {
  size_t cols_padded = AVX_PADDING(cols);
  size_t cols_padded_ymm = cols_padded / 256;
  size_t cols_padded_word = cols_padded / 64;

  for (size_t r = rstart; r < rstop; ++r) {
    size_t x0 = rev[k][matrix_read_bits(M, r, cstart, k)];
    if (x0)
      xor_avx1((uint8_t *)(T + x0 * cols_padded_word), (uint8_t *)M[r],
               (uint8_t *)M[r], cols_padded_ymm);
  }
}

size_t matrix_echelonize_partial(matrix_t M, size_t rows, size_t cols, size_t k,
                                 size_t rstop, uint64_t *xor_rows, int **rev,
                                 int **diff) {
  size_t kk = k;

  size_t r = 0;
  size_t c = 0;

  while (c < rstop) {
    if (c + kk > rstop) {
      kk = rstop - c;
    }
    size_t kbar;
    kbar = matrix_gauss_submatrix(M, r, c, rows, cols, kk);
    if (kk != kbar) break;

    if (kbar > 0) {
      matrix_make_table(M, r, cols, kbar, xor_rows, diff);
      matrix_process_rows(M, r + kbar, c, rows, kbar, cols, xor_rows, rev);
      matrix_process_rows(M, 0, c, r, kbar, cols, xor_rows, rev);
    }

    r += kbar;
    c += kbar;
  }

  return r;
}

size_t matrix_opt_k(size_t a, size_t b) {
  size_t n = (a < b) ? a : b;
  size_t res = (int)(0.75 * (1 + flb(n)));
  res = (1 > res) ? 1 : res;
  res = (MAX_K < res) ? MAX_K : res;
  return res;
}
