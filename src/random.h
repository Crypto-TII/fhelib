#ifndef TIIFHE_RANDOM_H
#define TIIFHE_RANDOM_H

#include <openssl/rand.h>
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct tiifhe_seed tiifhe_Seed;
struct tiifhe_seed {
	uint64_t state[25];
	unsigned char value[16];
	size_t pos;
	int init;
};

#include "config.h"

/* Return len uniform random bytes to the buffer out using seed. */
void tiifhe_random_bytes(unsigned char *out, size_t len, tiifhe_Seed *seed);

/* Return len bytes uniform from {-1, 0, 0, 1} given a seed to the buffer out
 * of signed 8-bit integers (i8) or singed 64-bit integers (i64). */
void tiifhe_random_cbd1_i8(int8_t *out, size_t len, tiifhe_Seed *seed);
void tiifhe_random_cbd1_i64(int64_t *out, size_t len, tiifhe_Seed *seed);

/* Return len random bytes from a centered binomial distribution with variance
 * sigma^2 = 21 / 2 = 10.5 given a seed to the buffer out of signed 8-bit
 * integers (i8) or signed 64-bit integers (i64). */
void tiifhe_random_cbd21_i8(int8_t *out, size_t len, tiifhe_Seed *seed);
void tiifhe_random_cbd21_i64(int64_t *out, size_t len, tiifhe_Seed *seed);

/* Return a uniform random seed using the OpenSSL RAND_bytes function. */
void tiifhe_random_seed(tiifhe_Seed *seed);

/* Return len uniform random elements in Z_mod to the buffer out of
 * signed (i64) or unsigned (u64) 64-bit integers given a seed. Uses
 * rejection sampling. */
void tiifhe_random_uniform_i64(int64_t *out, size_t len, uint64_t mod, tiifhe_Seed *seed);
void tiifhe_random_uniform_u64(uint64_t *out, size_t len, uint64_t mod, tiifhe_Seed *seed);

#endif /* TIIFHE_RANDOM_H */
