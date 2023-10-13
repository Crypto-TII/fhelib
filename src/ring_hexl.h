#ifndef TIIFHE_RING_H
#define TIIFHE_RING_H

#include <assert.h>
#include <gmp.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"

typedef uint64_t           tiifhe_Digit;
typedef struct tiifhe_mod  tiifhe_Mod;
typedef struct tiifhe_poly tiifhe_Poly;

struct tiifhe_mod {
	uint64_t value;
	uint64_t t;
	uint64_t P;
	void *ntt;
	size_t idx;
};

struct tiifhe_poly {
	uint64_t value[TIIFHE_N];
};

extern tiifhe_Mod tiifhe_q[TIIFHE_QPLEN];
extern const tiifhe_Digit tiifhe_rns[TIIFHE_QPLEN];

tiifhe_Digit tiifhe_mod_const(size_t idx, const tiifhe_Digit op);
void         tiifhe_mod_deinit(size_t idx);
void         tiifhe_mod_init_mpz(size_t idx, mpz_t t);
void         tiifhe_mod_init_u64(size_t idx, uint64_t t);
tiifhe_Digit tiifhe_mod_inv(size_t idx, const tiifhe_Digit op);
tiifhe_Digit tiifhe_mod_mul(size_t idx, const tiifhe_Digit op1, const tiifhe_Digit op2);
tiifhe_Digit tiifhe_mod_neg(size_t idx, const tiifhe_Digit op);

void tiifhe_poly_add(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Poly *op2);
void tiifhe_poly_addmul(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Poly *op2, const tiifhe_Poly *op3, const tiifhe_Poly *op4);
void tiifhe_poly_cmod(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *p);
void tiifhe_poly_deinit(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *p);
void tiifhe_poly_ext(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly **p, const size_t *base, const size_t *drops, size_t len);
void tiifhe_poly_init(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *p);
void tiifhe_poly_init_error(size_t idx, tiifhe_Poly *p, const tiifhe_Poly *e);
void tiifhe_poly_intt(size_t idx, tiifhe_Poly *p);
void tiifhe_poly_mod_mpz(size_t idx, tiifhe_Poly *rop, const mpz_t *p);
void tiifhe_poly_mod_u64(size_t idx, tiifhe_Poly *rop, const uint64_t *p);
void tiifhe_poly_mul(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Poly *op2);
void tiifhe_poly_muladd(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Poly *op2, const tiifhe_Poly *op3);
void tiifhe_poly_mulc(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Digit op2);
void tiifhe_poly_mulcadd(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Digit op2, const tiifhe_Poly *op3);
void tiifhe_poly_neg(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *p);
void tiifhe_poly_ntt(size_t idx, tiifhe_Poly *p);
void tiifhe_poly_rot_intt(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *p, size_t steps);
void tiifhe_poly_rot_ntt(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *p, size_t steps);
void tiifhe_poly_sample_error(tiifhe_Poly *e, tiifhe_Seed *seed);
void tiifhe_poly_sample_secret(tiifhe_Poly *s, tiifhe_Seed *seed);
void tiifhe_poly_sample_uniform(size_t idx, tiifhe_Poly *p, tiifhe_Seed *seed);
void tiifhe_poly_sub(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Poly *op2);

#endif /* TIIFHE_RING_H */
