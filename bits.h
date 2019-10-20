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
#ifndef BITS_H
#define BITS_H
#include <stdint.h>

/*
 * Round a number of bits up to the next multiple of 256 so that data fits in
 * an AVX register.
 */
#define AVX_PADDING(len) (((len + 255) / 256) * 256)

uint64_t popcount(const uint64_t *buf, unsigned len, unsigned max);

void xor_avx1(uint8_t *x, uint8_t *y, uint8_t *z, unsigned n);
void xor_avx2(uint8_t *x, uint8_t *y1, uint8_t *y2, uint8_t *z, unsigned n);
void xor_avx3(uint8_t *x, uint8_t *y1, uint8_t *y2, uint8_t *y3, uint8_t *z,
              unsigned n);
void xor_avx4(uint8_t *x, uint8_t *y1, uint8_t *y2, uint8_t *y3, uint8_t *y4,
              uint8_t *z, unsigned n);

void xor_bcast_32(uint32_t x, uint8_t *y, uint8_t *z, unsigned n);
void xor_bcast_64(uint64_t x, uint8_t *y, uint8_t *z, unsigned n);
#endif /* BITS_H */
