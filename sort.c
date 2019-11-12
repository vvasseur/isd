#include "sort.h"

#define BITS (sizeof(SORT_TYPE) * 8)
#define RADIX 8
#define BUCKETS (1L << RADIX)
#define DIGIT(A, B) (((A) >> (BITS - ((B) + 1) * RADIX)) & (BUCKETS - 1))

/* LSD radix sort */
void sort(SORT_TYPE *array, size_t *idx, SORT_TYPE *aux, size_t *aux2,
          size_t len) {
  for (size_t w = BITS / RADIX; w-- > 0;) {
    size_t count[BUCKETS + 1] = {0};

    for (size_t i = 0; i < len; ++i) ++count[DIGIT(array[i], w) + 1];

    for (size_t j = 1; j < BUCKETS - 1; ++j) count[j + 1] += count[j];

    for (size_t i = 0; i < len; ++i) {
      size_t cnt = count[DIGIT(array[i], w)];
      aux[cnt] = array[i];
      aux2[cnt] = idx[i];
      ++count[DIGIT(array[i], w)];
    }

    for (size_t i = 0; i < len; ++i) {
      array[i] = aux[i];
      idx[i] = aux2[i];
    }
  }
}
