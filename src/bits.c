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
#include "bits.h"

#ifdef __aarch64__
// https://github.com/kunpengcompute/AvxToNeon
#include "../../AvxToNeon/avx2neon.h"
#else
#include <immintrin.h>
#endif

uint64_t popcount(const uint64_t *buf, unsigned len, unsigned max) {
  uint64_t cnt = 0;
  for (unsigned i = 0; i < len && cnt <= max; i += 4) {
    cnt += __builtin_popcountll(buf[i]);
    cnt += __builtin_popcountll(buf[i + 1]);
    cnt += __builtin_popcountll(buf[i + 2]);
    cnt += __builtin_popcountll(buf[i + 3]);
  }
  return cnt;
}

/* ceil(log2(x)) */
unsigned clb(unsigned long x) {
  if (x <= 1) return 1;
  return (8 * sizeof(unsigned long)) - __builtin_clzl(x - 1);
}

unsigned flb(unsigned long x) {
  if (x < 1) return 0;
  return (8 * sizeof(unsigned long)) - __builtin_clzl(x) - 1;
}

void xor_bcast_16(uint16_t x, uint8_t *y, uint8_t *z, unsigned n) {
  __m256i vec_x = _mm256_set1_epi16(x);

  for (unsigned i = 0; i < n; i += 1) {
    __m256i vec_y = (((__m256i *)y)[i]);
    __m256i vec_z = _mm256_xor_si256(vec_x, vec_y );
    (((__m256i *)z)[i]) = vec_z;
  }
}

void xor_bcast_32(uint32_t x, uint8_t *y, uint8_t *z, unsigned n) {
  __m256i vec_x = _mm256_set1_epi32(x);

  for (unsigned i = 0; i < n; i += 1) {
    __m256i vec_y = (((__m256i *)y)[i]);
    __m256i vec_z = _mm256_xor_si256(vec_x, vec_y );
    (((__m256i *)z)[i]) = vec_z;
  }
}

void xor_bcast_64(uint64_t x, uint8_t *y, uint8_t *z, unsigned n) {
  __m256i vec_x = _mm256_set1_epi64x(x);

  for (unsigned i = 0; i < n; i += 1) {
    __m256i vec_y = (((__m256i *)y)[i]);
    __m256i vec_z = _mm256_xor_si256(vec_x, vec_y );
     (((__m256i *)z)[i]) = vec_z;
  }
}

void xor_avx1(uint8_t *x, uint8_t *y, uint8_t *z, unsigned n) {
  for (unsigned i = 0; i < n; i += 1) {
    __m256i vec_x = (((__m256i *)x)[i]);
    __m256i vec_y = (((__m256i *)y)[i]);
    __m256i vec_z = _mm256_xor_si256(vec_x, vec_y );
    (((__m256i *)z)[i]) = vec_z;
  }
}

void xor_avx2(uint8_t *x, uint8_t *y1, uint8_t *y2, uint8_t *z, unsigned n) {
  for (unsigned i = 0; i < n; i += 1) {
    __m256i vec_x = (((__m256i *)x)[i]);
    __m256i vec_y1 = (((__m256i *)y1)[i]);
    vec_x = _mm256_xor_si256(vec_y1, vec_x );

    __m256i vec_y2 = (((__m256i *)y2)[i]);
    __m256i vec_z = _mm256_xor_si256(vec_y2, vec_x );
    (((__m256i *)z)[i]) = vec_z;
  }
}

void xor_avx3(uint8_t *x, uint8_t *y1, uint8_t *y2, uint8_t *y3, uint8_t *z,
              unsigned n) {
  for (unsigned i = 0; i < n; i += 1) {
    __m256i vec_x = (((__m256i *)x)[i]);
    __m256i vec_y1 = (((__m256i *)y1)[i]);
    vec_x = _mm256_xor_si256(vec_y1, vec_x );

    __m256i vec_y2 = (((__m256i *)y2)[i]);
    vec_x = _mm256_xor_si256(vec_y2, vec_x );

    __m256i vec_y3 = (((__m256i *)y3)[i]);
    __m256i vec_z = _mm256_xor_si256(vec_y3, vec_x );
    (((__m256i *)z)[i]) = vec_z;
  }
}

void xor_avx4(uint8_t *x, uint8_t *y1, uint8_t *y2, uint8_t *y3, uint8_t *y4,
              uint8_t *z, unsigned n) {
  for (unsigned i = 0; i < n; i += 1) {
    __m256i vec_x = (((__m256i *)x)[i]);
    __m256i vec_y1 = (((__m256i *)y1)[i]);
    vec_x = _mm256_xor_si256(vec_y1, vec_x );

    __m256i vec_y2 = (((__m256i *)y2)[i]);
    vec_x = _mm256_xor_si256(vec_y2, vec_x );

    __m256i vec_y3 = (((__m256i *)y3)[i]);
    vec_x = _mm256_xor_si256(vec_y3, vec_x );

    __m256i vec_y4 = (((__m256i *)y4)[i]);
    __m256i vec_z = _mm256_xor_si256(vec_y4, vec_x );
    (((__m256i *)z)[i]) = vec_z;
  }
}

void copy_avx(uint8_t *dst, const uint8_t *src, unsigned n) {
    for (unsigned i = 0; i < n; ++i) {
        __m256i vec_src = _mm256_load_si256((__m256i *)&src[i * 32]);
        _mm256_store_si256((__m256i *)&dst[i * 32], vec_src);
    }
}
