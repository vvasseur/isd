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
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "dumer.h"
#ifndef BENCHMARK
#define BENCHMARK 0
#endif

struct timespec timer_start() {
  struct timespec start_time;
  clock_gettime(CLOCK_MONOTONIC, &start_time);
  return start_time;
}

long timer_end(struct timespec start_time) {
  struct timespec end_time;
  clock_gettime(CLOCK_MONOTONIC, &end_time);
  long diff = (end_time.tv_sec - start_time.tv_sec) * (long)1e9 +
              (end_time.tv_nsec - start_time.tv_nsec);
  return diff;
}

static void skip_comment(FILE *file, int *c) {
  *c = getc(file);
  if (*c == '#')
    while ((*c = getc(file)) != EOF && *c != '\n')
      ;
  if (*c == '\n') *c = getc(file);
}

static int read_int(FILE *file, int *c, size_t *n) {
  if (*c == EOF) return 0;
  while (*c != EOF && *c != '\n') {
    if (*c >= '0' && *c <= '9') {
      *n *= 10;
      *n += *c - '0';
    } else
      return 0;
    *c = getc(file);
  }
  return 1;
}

static int read_bin_vector(FILE *file, int *c, uint8_t *m, size_t *len) {
  if (*c == EOF) return 0;
  *len = 0;
  while (*c != EOF && *c != '\n') {
    if (*c == '1' || *c == '0') {
      m[(*len)++] = *c - '0';
    } else
      return 0;
    *c = getc(file);
  }
  return 1;
}

static int read_bin_matrix(FILE *file, int rows, int *c, uint8_t *m,
                           size_t *len) {
  if (*c == EOF) return 0;
  *len = 0;
  while (*c != EOF) {
    if (*c == '1' || *c == '0') {
      m[(*len)++] = *c - '0';
    } else if (*c == '\n')
      --rows;
    else {
      return 0;
    }
    if (!rows) break;
    *c = getc(file);
  }
  return 1;
}

static int parse_input_sd(char *filename, size_t *n, size_t *k, size_t *w,
                          uint8_t **mat_h, size_t *len_h, uint8_t **mat_s,
                          size_t *len_s) {
  int ret = 0;
  *n = 0;
  *k = 0;
  *w = 0;
  size_t seed = 0;
  int c;
  FILE *file;
  file = fopen(filename, "r");
  if (file) {
    skip_comment(file, &c);

    /* Read n. */
    if (!read_int(file, &c, n)) goto end;
    *k = *n / 2;

    skip_comment(file, &c);

    /* Read seed. */
    if (!read_int(file, &c, &seed)) goto end;

    skip_comment(file, &c);

    /* Read w. */
    if (!read_int(file, &c, w)) goto end;

    *mat_h = malloc(*k * *k * sizeof(uint8_t));

    if (!mat_h) goto end;

    skip_comment(file, &c);

    /* Read h. */
    if (!read_bin_matrix(file, *k, &c, *mat_h, len_h)) goto end;

    *mat_s = malloc(*k * sizeof(uint8_t));

    if (!mat_s) goto end;

    skip_comment(file, &c);

    /* Read s. */
    if (!read_bin_vector(file, &c, *mat_s, len_s)) goto end;
  } else {
    return 0;
  }
  ret = 1;

end:
  fclose(file);
  return ret;
}

static int parse_input_go(char *filename, size_t *n, size_t *k, size_t *w,
                          uint8_t **mat_h, size_t *len_h, uint8_t **mat_s,
                          size_t *len_s) {
  int ret = 0;
  *n = 0;
  *k = 0;
  *w = 0;
  int c;
  FILE *file;
  file = fopen(filename, "r");
  if (file) {
    skip_comment(file, &c);

    /* Read n. */
    if (!read_int(file, &c, n)) goto end;

    skip_comment(file, &c);

    /* Read k. */
    if (!read_int(file, &c, k)) goto end;

    skip_comment(file, &c);

    /* Read w. */
    if (!read_int(file, &c, w)) goto end;

    *mat_h = malloc(*n * *k * sizeof(uint8_t));

    if (!mat_h) goto end;

    skip_comment(file, &c);

    /* Read h. */
    if (!read_bin_matrix(file, *k, &c, *mat_h, len_h)) goto end;

    *mat_s = malloc((*n - *k) * sizeof(uint8_t));

    if (!mat_s) goto end;

    skip_comment(file, &c);

    /* Read s. */
    if (!read_bin_vector(file, &c, *mat_s, len_s)) goto end;
  } else {
    return 0;
  }
  ret = 1;

end:
  fclose(file);
  return ret;
}

static int parse_input_qc(char *filename, size_t *n, size_t *k, size_t *w,
                          uint8_t **mat_h, size_t *len_h, uint8_t **mat_s,
                          size_t *len_s) {
  int ret = 0;
  *n = 0;
  *k = 0;
  *w = 0;
  int c;
  FILE *file;
  file = fopen(filename, "r");
  if (file) {
    skip_comment(file, &c);

    /* Read n. */
    if (!read_int(file, &c, n)) goto end;
    *k = *n / 2;

    skip_comment(file, &c);

    /* Read w. */
    if (!read_int(file, &c, w)) goto end;

    *mat_h = malloc(*k * sizeof(uint8_t));

    if (!mat_h) goto end;

    skip_comment(file, &c);

    /* Read h. */
    if (!read_bin_vector(file, &c, *mat_h, len_h)) goto end;

    *mat_s = malloc(*k * sizeof(uint8_t));

    if (!mat_s) goto end;

    skip_comment(file, &c);

    /* Read s. */
    if (!read_bin_vector(file, &c, *mat_s, len_s)) goto end;
  } else {
    return 0;
  }
  ret = 1;

end:
  fclose(file);
  return ret;
}

static int parse_input_lw(char *filename, size_t *n, size_t *k, size_t *w,
                          uint8_t **mat_h, size_t *len_h) {
  int ret = 0;
  *n = 0;
  *k = 0;
  *w = 0;
  size_t seed = 0;
  int c;
  FILE *file;
  file = fopen(filename, "r");
  if (file) {
    skip_comment(file, &c);

    /* Read n. */
    if (!read_int(file, &c, n)) goto end;
    *k = *n / 2;

    skip_comment(file, &c);

    /* Read seed. */
    if (!read_int(file, &c, &seed)) goto end;

    *mat_h = malloc(*k * *k * sizeof(uint8_t));

    if (!mat_h) goto end;

    skip_comment(file, &c);

    /* Read h. */
    if (!read_bin_matrix(file, *k, &c, *mat_h, len_h)) goto end;
  } else {
    return 0;
  }
  ret = 1;

end:
  fclose(file);
  return ret;
}

int main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(stderr,
            "Usage: %s [N_THREADS] [TYPE] [FILE]\n"
            "\n"
            "where TYPE is:\n"
            "         SD for syndrome decoding\n"
            "         QC for quasi-cyclic syndrome decoding\n"
            "         GO for Goppa codes syndrome decoding\n"
            "         LW for low-weight codeword finding\n",
            argv[0]);
    exit(EXIT_FAILURE);
  }

  enum type current_type;

  if (!strcmp(argv[2], "QC"))
    current_type = QC;
  else if (!strcmp(argv[2], "SD"))
    current_type = SD;
  else if (!strcmp(argv[2], "GO"))
    current_type = GO;
  else if (!strcmp(argv[2], "LW"))
    current_type = LW;
  else {
    fprintf(stderr, "Check your arguments!\n");
    exit(EXIT_FAILURE);
  }

#if !(DUMER_LW)
  if (current_type == LW) {
    fprintf(stderr, "No syndrome to decode.\n");
    exit(EXIT_FAILURE);
  }
#endif
#if DUMER_DOOM
  if (current_type != QC) {
    fprintf(stderr,
            "Using DOOM in a non quasi-cyclic setting will "
            "most likely not give any meaningful result!\n");
  }
#endif

  int n_threads = atoi(argv[1]);
  if (n_threads < 0) {
    fprintf(stderr, "N_THREADS should be greater than 0.\n");
    exit(EXIT_FAILURE);
  }

  size_t n, k, w;
  uint8_t *mat_h = NULL;
  uint8_t *mat_s = NULL;
  size_t len_h, len_s;

  int parsed;
  if (current_type == QC)
    parsed =
        parse_input_qc(argv[3], &n, &k, &w, &mat_h, &len_h, &mat_s, &len_s);
  else if (current_type == SD)
    parsed =
        parse_input_sd(argv[3], &n, &k, &w, &mat_h, &len_h, &mat_s, &len_s);
  else if (current_type == GO)
    parsed =
        parse_input_go(argv[3], &n, &k, &w, &mat_h, &len_h, &mat_s, &len_s);
  else if (current_type == LW)
    parsed = parse_input_lw(argv[3], &n, &k, &w, &mat_h, &len_h);

  if (!parsed) {
    fprintf(stderr, "Error parsing file.\n");
    exit(EXIT_FAILURE);
  }

  size_t r = n - k;

  printf("n=%ld ", n);
  printf("k=%ld ", k);
  printf("w=%ld\n", w);
  printf("l=%ld ", DUMER_L);
  printf("p=%ld ", DUMER_P);
  printf("epsilon=%ld ", DUMER_EPS);
  printf("doom=%d\n", DUMER_DOOM);

  /* Birthday decoding */
  size_t n1 = (k + DUMER_L) / 2;
  size_t n2 = k + DUMER_L - n1;
  if (DUMER_EPS > n2 || DUMER_EPS > n1) {
    fprintf(stderr, "Please lower DUMER_EPS.\n");
    exit(EXIT_FAILURE);
  }

  /* Data shared by all threads and computed only once */
  shr_t shr = alloc_shr(n1, n2);
  if (!shr) {
    fprintf(stderr, "Allocation error.\n");
    exit(EXIT_FAILURE);
  }
  init_shr(shr, n, k, n1, n2);

#if (BENCHMARK) <= 0
#pragma omp parallel num_threads(n_threads)
  {
    isd_t isd = alloc_isd(n, k, r, n1, n2, shr->nb_combinations1, shr->k_opt);
    if (!isd) {
      fprintf(stderr, "Allocation error.\n");
      exit(EXIT_FAILURE);
    }
    init_isd(isd, current_type, n, k, w, mat_h, mat_s);

    while (1) {
      int found = dumer(n, k, r, n1, n2, shr, isd);
      if (found) {
        print_solution(n, isd);
#if !(DUMER_LW)
        exit(EXIT_SUCCESS);
#endif
      }
    }
    free_isd(isd, r, n);
  }
#else
  isd_t *isd = malloc(n_threads * sizeof(isd_t));
  if (!isd) {
    fprintf(stderr, "Allocation error.\n");
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < n_threads; i++) {
    isd[i] = alloc_isd(n, k, r, n1, n2, shr->nb_combinations1, shr->k_opt);
    if (!isd[i]) {
      fprintf(stderr, "Allocation error.\n");
      exit(EXIT_FAILURE);
    }
    init_isd(isd[i], current_type, n, k, w, mat_h, mat_s);
  }
  struct timespec vartime = timer_start();  // begin a timer called 'vartime'
#pragma omp parallel num_threads(n_threads)
  {
    int i = omp_get_thread_num();
    for (size_t N = 0; N < (BENCHMARK + i) / n_threads; ++N) {
      dumer(n, k, r, n1, n2, shr, isd[i]);
    }
  }
  long time_elapsed_nanos = timer_end(vartime);
  printf("%ld\n", time_elapsed_nanos);
  for (int i = 0; i < n_threads; i++) {
    free_isd(isd[i], r, n);
  }
#endif

  if (mat_h) free(mat_h);
  if (mat_h) free(mat_s);
  free_shr(shr);
  exit(EXIT_SUCCESS);
}
