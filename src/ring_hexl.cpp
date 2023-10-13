#include <gmp.h>
extern "C" {
	#include "ring_hexl.h"
}

#include "hexl/hexl.hpp"

tiifhe_Mod tiifhe_q[TIIFHE_QPLEN];
const tiifhe_Digit tiifhe_rns[TIIFHE_QPLEN] = {
	TIIFHE_Q,
	TIIFHE_P
};

tiifhe_Digit
tiifhe_mod_const(size_t idx, const tiifhe_Digit op)
{

	(void)idx;

	return op;
}

void
tiifhe_mod_deinit(size_t idx)
{

	delete (intel::hexl::NTT *)tiifhe_q[idx].ntt;
}

void
tiifhe_mod_init_mpz(size_t idx, mpz_t t)
{
	mpz_t tmp;

	mpz_init(tmp);

	mpz_mod_ui(tmp, t, tiifhe_rns[idx]);
	tiifhe_mod_init_u64(idx, mpz_get_ui(tmp));

	mpz_clear(tmp);
}

void
tiifhe_mod_init_u64(size_t idx, uint64_t t)
{
	int64_t mod, P;
	size_t i;

	mod = tiifhe_rns[idx];
	tiifhe_q[idx].value = mod;
	tiifhe_q[idx].t = t;

	P = 1;
	for (i = TIIFHE_QLEN; i < TIIFHE_QPLEN; ++i)
		P = intel::hexl::MultiplyMod(P, tiifhe_rns[i], mod);
	tiifhe_q[idx].P = P;

	tiifhe_q[idx].ntt = (void *)(new intel::hexl::NTT(TIIFHE_N, mod));
	tiifhe_q[idx].idx = idx;
}

tiifhe_Digit
tiifhe_mod_inv(size_t idx, const tiifhe_Digit op)
{

	return intel::hexl::InverseMod(op, tiifhe_rns[idx]);
}

tiifhe_Digit
tiifhe_mod_mul(size_t idx, const tiifhe_Digit op1, const tiifhe_Digit op2)
{

	return intel::hexl::MultiplyMod(op1, op2, tiifhe_rns[idx]);
}

tiifhe_Digit
tiifhe_mod_neg(size_t idx, const tiifhe_Digit op)
{

	return tiifhe_rns[idx] - op;
}

void
tiifhe_poly_add(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Poly *op2)
{

	intel::hexl::EltwiseAddMod(rop->value, op1->value, op2->value, TIIFHE_N, tiifhe_rns[idx]);
}

void
tiifhe_poly_addmul(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Poly *op2, const tiifhe_Poly *op3, const tiifhe_Poly *op4)
{
	tiifhe_Poly tmp;

	intel::hexl::EltwiseAddMod(tmp.value, op1->value, op2->value, TIIFHE_N, tiifhe_rns[idx]);
	intel::hexl::EltwiseAddMod(rop->value, op3->value, op4->value, TIIFHE_N, tiifhe_rns[idx]);
	intel::hexl::EltwiseMultMod(rop->value, rop->value, tmp.value, TIIFHE_N, tiifhe_rns[idx], 1);
}

void
tiifhe_poly_cmod(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *p)
{

	(void)idx;
	memcpy(rop, p, sizeof *rop);
}

void
tiifhe_poly_deinit(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *p)
{

	(void)idx;
	memcpy(rop, p, sizeof *rop);
}

void
tiifhe_poly_ext(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly **p, const size_t *base, const size_t *drops, size_t len)
{
	tiifhe_Digit inv[TIIFHE_QPLEN], div[TIIFHE_QPLEN], tmp;
	size_t i;

	for (i = 0; i < len; ++i) {
		size_t ii;

		inv[i] = 1;
		for (ii = 0; ii < len; ++ii) {
			if (drops[ii] != 0)
				continue;

			tmp = intel::hexl::InverseMod(tiifhe_rns[base[ii]], tiifhe_rns[base[i]]);
			inv[i] = intel::hexl::MultiplyMod(inv[i], tmp, tiifhe_rns[base[i]]);
		}

		div[i] = intel::hexl::InverseMod(tiifhe_rns[base[i]], tiifhe_rns[idx]);
		for (ii = 0; ii < len; ++ii) {
			if (idx == base[i] && idx == base[ii])
				continue;

			if (drops[ii] != 0)
				continue;

			div[i] = intel::hexl::MultiplyMod(div[i], tiifhe_rns[base[ii]], tiifhe_rns[idx]);
		}
	}

	memset(rop, 0, sizeof *rop);
	for (i = 0; i < len; ++i) {
		tiifhe_Poly tmp;

		if (drops[i] != 0)
			continue;

		tiifhe_poly_mulc(base[i], &tmp, p[i], inv[i]);
		tiifhe_poly_mulcadd(idx, rop, &tmp, div[i], rop);
	}
}

void
tiifhe_poly_init(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *p)
{
	size_t j;

	/* check for arithmetic shift */
	assert(-1 >> 1 == -1);

	for (j = 0; j < TIIFHE_N; ++j) {
		int64_t val = p->value[j];
		rop->value[j] = val + ((val >> 63) & tiifhe_rns[idx]);
	}
}

void
tiifhe_poly_init_error(size_t idx, tiifhe_Poly *p, const tiifhe_Poly *e)
{

	tiifhe_poly_init(idx, p, e);
	tiifhe_poly_mulc(idx, p, p, tiifhe_q[idx].t);
	tiifhe_poly_ntt(idx, p);
}

void
tiifhe_poly_intt(size_t idx, tiifhe_Poly *p)
{
	intel::hexl::NTT *ntt;

	ntt = (intel::hexl::NTT *)tiifhe_q[idx].ntt;
	ntt->ComputeInverse(p->value, p->value, 1, 1);
}

void
tiifhe_poly_mod_mpz(size_t idx, tiifhe_Poly *rop, const mpz_t *p)
{
	mpz_t tmp;
	size_t j;

	mpz_init(tmp);

	for (j = 0; j < TIIFHE_N; ++j) {
		mpz_mod_ui(tmp, p[j], tiifhe_rns[idx]);
		rop->value[j] = mpz_get_ui(tmp);
	}

	mpz_clear(tmp);
}

void
tiifhe_poly_mod_u64(size_t idx, tiifhe_Poly *rop, const uint64_t *p)
{

	intel::hexl::EltwiseReduceMod(rop->value, p, TIIFHE_N, tiifhe_rns[idx], 1, 1);
}

void
tiifhe_poly_mul(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Poly *op2)
{

	intel::hexl::EltwiseMultMod(rop->value, op1->value, op2->value, TIIFHE_N, tiifhe_rns[idx], 1);
}

void
tiifhe_poly_muladd(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Poly *op2, const tiifhe_Poly *op3)
{
	tiifhe_Poly tmp;

	intel::hexl::EltwiseMultMod(tmp.value, op1->value, op2->value, TIIFHE_N, tiifhe_rns[idx], 1);
	intel::hexl::EltwiseAddMod(rop->value, tmp.value, op3->value, TIIFHE_N, tiifhe_rns[idx]);
}

void
tiifhe_poly_mulc(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Digit op2)
{

	intel::hexl::EltwiseFMAMod(rop->value, op1->value, op2, 0, TIIFHE_N, tiifhe_rns[idx], 1);
}

void
tiifhe_poly_mulcadd(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Digit op2, const tiifhe_Poly *op3)
{

	intel::hexl::EltwiseFMAMod(rop->value, op1->value, op2, op3->value, TIIFHE_N, tiifhe_rns[idx], 1);
}

void
tiifhe_poly_neg(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *p)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j)
		rop->value[j] = tiifhe_rns[idx] - p->value[j];
}

void
tiifhe_poly_ntt(size_t idx, tiifhe_Poly *p)
{
	intel::hexl::NTT *ntt;

	ntt = (intel::hexl::NTT *)tiifhe_q[idx].ntt;
	ntt->ComputeForward(p->value, p->value, 1, 1);
}

void
tiifhe_poly_rot_intt(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *p, size_t steps)
{
	size_t step, i, j;

	assert(rop != p);

	steps %= TIIFHE_N / 2;
	if (steps == 0) {
		memcpy(rop, p, sizeof *rop);
		return;
	}

	step = 3;
	for (i = 1; i < steps; ++i) {
		step *= 3;
		step %= 2 * TIIFHE_N;
	}

	i = 0;
	for (j = 0; j < TIIFHE_N; ++j) {
		tiifhe_Digit tmp = p->value[j];

		if (i / TIIFHE_N & 1)
			tmp = tiifhe_rns[idx] - tmp;
		rop->value[i % TIIFHE_N] = tmp;

		i += step;
		i %= 2 * TIIFHE_N;
	}
}

void
tiifhe_poly_rot_ntt(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *p, size_t steps)
{
	size_t dlog, start, step, i, j;

	assert(rop != p);
	(void)idx;

	dlog = tiifhe_util_msb64(TIIFHE_N);

	steps %= TIIFHE_N / 2;
	if (steps == 0) {
		memcpy(rop, p, sizeof *rop);
		return;
	}

	start = 1;
	step = 3;
	for (i = 1; i < steps; ++i) {
		start = (start * 3 + 1) % TIIFHE_N;
		step = (step * 3) % TIIFHE_N;
	}

	i = start;
	for (j = 0; j < TIIFHE_N; ++j) {
		size_t idx1 = tiifhe_util_bitrev(i, dlog);
		size_t idx2 = tiifhe_util_bitrev(j, dlog);
		rop->value[idx2] = p->value[idx1];
		i = (i + step) % TIIFHE_N;
	}
}

void
tiifhe_poly_sample_error(tiifhe_Poly *e, tiifhe_Seed *seed)
{

	tiifhe_random_cbd21_i64((int64_t *)e->value, TIIFHE_N, seed);
}

void
tiifhe_poly_sample_secret(tiifhe_Poly *s, tiifhe_Seed *seed)
{

	tiifhe_random_cbd1_i64((int64_t *)s->value, TIIFHE_N, seed);
}

void
tiifhe_poly_sample_uniform(size_t idx, tiifhe_Poly *p, tiifhe_Seed *seed)
{

	tiifhe_random_uniform_u64(p->value, TIIFHE_N, tiifhe_rns[idx], seed);
}

void
tiifhe_poly_sub(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Poly *op2)
{

	intel::hexl::EltwiseSubMod(rop->value, op1->value, op2->value, TIIFHE_N, tiifhe_rns[idx]);
}
