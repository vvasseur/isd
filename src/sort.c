#include "sort.h"

#define RADIX 8
#define BUCKETS (1L << RADIX)

/* LSD radix sort */
void sort(SORT_TYPE *restrict array, size_t *restrict idx,
          SORT_TYPE *restrict aux, size_t *restrict aux_idx, size_t len) {
  const size_t stride = SORT_WIDTH / RADIX;
  size_t count[8][BUCKETS + 1] = {0};

  {
    uint8_t *restrict array8 = ((uint8_t *)array);
    for (size_t i = 0; i < len; ++i) {
      ++count[0][array8[0] + 1];
#if SORT_WIDTH >= 16
      ++count[1][array8[1] + 1];
#endif
#if SORT_WIDTH >= 32
      ++count[2][array8[2] + 1];
      ++count[3][array8[3] + 1];
#endif
#if SORT_WIDTH == 64
      ++count[4][array8[4] + 1];
      ++count[5][array8[5] + 1];
      ++count[6][array8[6] + 1];
      ++count[7][array8[7] + 1];
#endif
      array8 += stride;
    }
  }

  for (size_t j = 1; j < BUCKETS - 1; ++j) {
    count[0][j + 1] += count[0][j];
#if SORT_WIDTH >= 16
    count[1][j + 1] += count[1][j];
#endif
#if SORT_WIDTH >= 32
    count[2][j + 1] += count[2][j];
    count[3][j + 1] += count[3][j];
#endif
#if SORT_WIDTH == 64
    count[4][j + 1] += count[4][j];
    count[5][j + 1] += count[5][j];
    count[6][j + 1] += count[6][j];
    count[7][j + 1] += count[7][j];
#endif
  }

  for (size_t w = 0; w < SORT_WIDTH / RADIX; w++) {
    uint8_t *restrict array8 = ((uint8_t *)array) + w;
    for (size_t i = 0; i < len; ++i) {
      uint8_t byte = array8[stride * i];
      size_t *cnt = &count[w][byte];
      aux[*cnt] = array[i];
      aux_idx[*cnt] = idx[i];
      ++(*cnt);
    }

    /* No need to copy the array as we always do an even number of
     * iterations when SORT_WIDTH >= 16. */
    {
      SORT_TYPE *swp = array;
      array = aux;
      aux = swp;
    }
    {
      size_t *swp = idx;
      idx = aux_idx;
      aux_idx = swp;
    }
  }
#if SORT_WIDTH == 8
  for (size_t i = 0; i < len; ++i) {
    aux[i] = array[i];
    aux_idx[i] = idx[i];
  }
#endif
}
