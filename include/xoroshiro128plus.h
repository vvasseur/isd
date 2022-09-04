#ifndef XOROSHIRO128PLUS_H
#define XOROSHIRO128PLUS_H
#include <stdint.h>
int seed_random(uint64_t *S0, uint64_t *S1);
uint64_t random_lim(uint64_t limit, uint64_t *S0, uint64_t *S1);
#endif /* XOROSHIRO128PLUS_H */
