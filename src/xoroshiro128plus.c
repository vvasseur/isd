/*  Written in 2016-2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright and
related and neighboring rights to this software to the public domain worldwide.
This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

#include <stdio.h>

#include "xoroshiro128plus.h"

/* This is xoroshiro128+ 1.0, our best and fastest small-state generator
   for floating-point numbers. We suggest to use its upper bits for
   floating-point generation, as it is slightly faster than xoroshiro128**. It
   passes all tests we are aware of except for the four lower bits, which might
   fail linearity tests (and just those), so if low linear complexity is not
   considered an issue (as it is usually the case) it can be used to generate
   64-bit outputs, too; moreover, this generator has a very mild Hamming-weight
   dependency making our test (http://prng.di.unimi.it/hwd.php) fail after 5 TB
   of output; we believe this slight bias cannot affect any application. If you
   are concerned, use xoroshiro128++, xoroshiro128** or xoshiro256+.

   We suggest to use a sign test to extract a random Boolean value, and right
   shifts to extract subsets of bits.

   The state must be seeded so that it is not everywhere zero. If you have a
   64-bit seed, we suggest to seed a splitmix64 generator and use its output to
   fill s.

   NOTE: the parameters (a=24, b=16, b=37) of this version give slightly better
   results in our test than the 2016 version (a=55, b=14, c=36).
*/

static inline uint64_t rotl(const uint64_t x, int k) {
  return (x << k) | (x >> (64 - k));
}

uint64_t random_uint64_t(uint64_t *S0, uint64_t *S1) {
  const uint64_t s0 = *S0;
  uint64_t s1 = *S1;
  const uint64_t result = s0 + s1;

  s1 ^= s0;
  *S0 = rotl(s0, 24) ^ s1 ^ (s1 << 16);  // a, b
  *S1 = rotl(s1, 37);                    // c

  return result;
}

/* This is the jump function for the generator. It is equivalent
   to 2^64 calls to next(); it can be used to generate 2^64 non-overlapping
   subsequences for parallel computations. */
void jump(uint64_t *S0, uint64_t *S1) {
  static const uint64_t JUMP[] = {0xdf900294d8f554a5, 0x170865df4b3201fc};

  uint64_t s0 = 0;
  uint64_t s1 = 0;
  for (unsigned i = 0; i < sizeof(JUMP) / sizeof(*JUMP); i++)
    for (int b = 0; b < 64; b++) {
      if (JUMP[i] & UINT64_C(1) << b) {
        s0 ^= *S0;
        s1 ^= *S1;
      }
      random_uint64_t(S0, S1);
    }

  *S0 = s0;
  *S1 = s1;
}

/* This is the long-jump function for the generator. It is equivalent to
   2^96 calls to next(); it can be used to generate 2^32 starting points, from
   each of which jump() will generate 2^32 non-overlapping subsequences for
   parallel distributed computations. */
void long_jump(uint64_t *S0, uint64_t *S1) {
  static const uint64_t LONG_JUMP[] = {0xd2a98b26625eee7b, 0xdddf9b1090aa7ac1};

  uint64_t s0 = 0;
  uint64_t s1 = 0;
  for (unsigned i = 0; i < sizeof(LONG_JUMP) / sizeof(*LONG_JUMP); i++)
    for (int b = 0; b < 64; b++) {
      if (LONG_JUMP[i] & UINT64_C(1) << b) {
        s0 ^= *S0;
        s1 ^= *S1;
      }
      random_uint64_t(S0, S1);
    }

  *S0 = s0;
  *S1 = s1;
}

int seed_random(uint64_t *S0, uint64_t *S1) {
  uint64_t new_s[2];
  FILE *urandom_fp;

  urandom_fp = fopen("/dev/urandom", "r");
  if (urandom_fp == NULL) return 0;
  if (fread(&new_s, 8, 2, urandom_fp) != 2) return 0;
  fclose(urandom_fp);

  *S0 = new_s[0];
  *S1 = new_s[1];

  return 1;
}

uint64_t random_lim(uint64_t limit, uint64_t *S0, uint64_t *S1) {
  uint64_t divisor = 0xffffffffffffffffUL / (limit + 1);
  uint64_t retval;

  do {
    retval = random_uint64_t(S0, S1) / divisor;
  } while (retval > limit);

  return retval;
}
