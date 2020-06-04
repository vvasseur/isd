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

#include <immintrin.h>

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

void xor_bcast_32(uint32_t x, uint8_t *y, uint8_t *z, unsigned n) {
  __m256 vec_x;
  asm volatile("vbroadcastss   %[x], %[vec_x]\n\t"
               : [ vec_x ] "=x"(vec_x)
               : [ x ] "xm"(x)
               :);
  for (unsigned i = 0; i < n; i += 1) {
    asm volatile("vxorps   %[vec_yi], %[vec_x], %[vec_z]\n\t"
                 : [ vec_z ] "=xm"(((__m256 *)z)[i])
                 : [ vec_x ] "%x"(vec_x), [ vec_yi ] "xm"(((__m256 *)y)[i])
                 :);
  }
}

void xor_bcast_64(uint64_t x, uint8_t *y, uint8_t *z, unsigned n) {
  __m256 vec_x;
  asm volatile("vbroadcastsd   %[x], %[vec_x]\n\t"
               : [ vec_x ] "=x"(vec_x)
               : [ x ] "xm"(x)
               :);
  for (unsigned i = 0; i < n; i += 1) {
    asm volatile("vxorps   %[vec_yi], %[vec_x], %[vec_z]\n\t"
                 : [ vec_z ] "=xm"(((__m256 *)z)[i])
                 : [ vec_x ] "%x"(vec_x), [ vec_yi ] "xm"(((__m256 *)y)[i])
                 :);
  }
}

void xor_avx1(uint8_t *x, uint8_t *y, uint8_t *z, unsigned n) {
  for (unsigned i = 0; i < n; i += 1) {
    __m256 vec_x;
    asm volatile("vmovaps   %[x], %[vec_x]\n\t"
                 : [ vec_x ] "=x"(vec_x)
                 : [ x ] "xm"(((__m256 *)x)[i])
                 :);
    asm volatile("vxorps   %[vec_yi], %[vec_x], %[vec_z]\n\t"
                 : [ vec_z ] "=xm"(((__m256 *)z)[i])
                 : [ vec_x ] "%x"(vec_x), [ vec_yi ] "xm"(((__m256 *)y)[i])
                 :);
  }
}

void xor_avx2(uint8_t *x, uint8_t *y1, uint8_t *y2, uint8_t *z, unsigned n) {
  for (unsigned i = 0; i < n; i += 1) {
    __m256 vec_x;
    asm volatile("vmovaps   %[x], %[vec_x]\n\t"
                 : [ vec_x ] "=x"(vec_x)
                 : [ x ] "xm"(((__m256 *)x)[i])
                 :);
    asm volatile("vxorps   %[vec_yi], %[vec_x], %[vec_x]\n\t"
                 : [ vec_x ] "+x"(vec_x)
                 : [ vec_yi ] "xm"(((__m256 *)y1)[i])
                 :);
    asm volatile("vxorps   %[vec_yi], %[vec_x], %[vec_z]\n\t"
                 : [ vec_z ] "=xm"(((__m256 *)z)[i])
                 : [ vec_x ] "%x"(vec_x), [ vec_yi ] "xm"(((__m256 *)y2)[i])
                 :);
  }
}

void xor_avx3(uint8_t *x, uint8_t *y1, uint8_t *y2, uint8_t *y3, uint8_t *z,
              unsigned n) {
  for (unsigned i = 0; i < n; i += 1) {
    __m256 vec_x;
    asm volatile("vmovaps   %[x], %[vec_x]\n\t"
                 : [ vec_x ] "=x"(vec_x)
                 : [ x ] "xm"(((__m256 *)x)[i])
                 :);
    asm volatile("vxorps   %[vec_yi], %[vec_x], %[vec_x]\n\t"
                 : [ vec_x ] "+x"(vec_x)
                 : [ vec_yi ] "xm"(((__m256 *)y1)[i])
                 :);
    asm volatile("vxorps   %[vec_yi], %[vec_x], %[vec_x]\n\t"
                 : [ vec_x ] "+x"(vec_x)
                 : [ vec_yi ] "xm"(((__m256 *)y2)[i])
                 :);
    asm volatile("vxorps   %[vec_yi], %[vec_x], %[vec_z]\n\t"
                 : [ vec_z ] "=xm"(((__m256 *)z)[i])
                 : [ vec_x ] "%x"(vec_x), [ vec_yi ] "xm"(((__m256 *)y3)[i])
                 :);
  }
}

void xor_avx4(uint8_t *x, uint8_t *y1, uint8_t *y2, uint8_t *y3, uint8_t *y4,
              uint8_t *z, unsigned n) {
  for (unsigned i = 0; i < n; i += 1) {
    __m256 vec_x;
    asm volatile("vmovaps   %[x], %[vec_x]\n\t"
                 : [ vec_x ] "=x"(vec_x)
                 : [ x ] "xm"(((__m256 *)x)[i])
                 :);
    asm volatile("vxorps   %[vec_yi], %[vec_x], %[vec_x]\n\t"
                 : [ vec_x ] "+x"(vec_x)
                 : [ vec_yi ] "xm"(((__m256 *)y1)[i])
                 :);
    asm volatile("vxorps   %[vec_yi], %[vec_x], %[vec_x]\n\t"
                 : [ vec_x ] "+x"(vec_x)
                 : [ vec_yi ] "xm"(((__m256 *)y2)[i])
                 :);
    asm volatile("vxorps   %[vec_yi], %[vec_x], %[vec_x]\n\t"
                 : [ vec_x ] "+x"(vec_x)
                 : [ vec_yi ] "xm"(((__m256 *)y3)[i])
                 :);
    asm volatile("vxorps   %[vec_yi], %[vec_x], %[vec_z]\n\t"
                 : [ vec_z ] "=xm"(((__m256 *)z)[i])
                 : [ vec_x ] "%x"(vec_x), [ vec_yi ] "xm"(((__m256 *)y4)[i])
                 :);
  }
}
