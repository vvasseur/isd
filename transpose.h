/*
   Copyright (c) 2021 Valentin Vasseur

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
#ifndef TRANSPOSE_H
#define TRANSPOSE_H
#include "light_m4ri/matrix.h"

void matrix_transpose(matrix_t At, const matrix_t A, const size_t nrows,
                      const size_t ncols);
void matrix_reverse_rows(matrix_t A, const size_t nrows);
void matrix_transpose_rev_cols(matrix_t At, matrix_t A, const size_t nrows,
                               const size_t ncols);
void matrix_transpose_rev_rows(matrix_t At, matrix_t A, const size_t nrows,
                               const size_t ncols);
#endif /* TRANSPOSE_H */
