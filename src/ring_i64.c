#include "ring_i64.h"

tiifhe_Mod   tiifhe_q[TIIFHE_QPLEN];
tiifhe_Digit tiifhe_rns_mods[TIIFHE_QPLEN][TIIFHE_QPLEN];
tiifhe_Digit tiifhe_rns_invs[TIIFHE_QPLEN][TIIFHE_QPLEN];
const tiifhe_Digit tiifhe_rns[TIIFHE_QPLEN] = {
	TIIFHE_Q,
	TIIFHE_P
};

__extension__ typedef          __int128  int128_t;
__extension__ typedef unsigned __int128 uint128_t;

#if defined(__x86_64__)
	#include <x86intrin.h>
	#define BSR64(x) __bsrq(x)
#else
	#define BSR64(x) tiifhe_util_msb64(x)
#endif

/* https://eprint.iacr.org/2018/039 */
static int64_t barrett_precomp(int64_t mod);
static int64_t barrett_reduce(int64_t a, int64_t mod, int64_t precomp);

/* idx-based functions */
static tiifhe_Digit digit_add(size_t idx, tiifhe_Digit a, tiifhe_Digit b);
static tiifhe_Digit digit_deinit(size_t idx, tiifhe_Digit a);
static tiifhe_Digit digit_init(size_t idx, tiifhe_Digit a);
static tiifhe_Digit digit_inv(size_t idx, tiifhe_Digit a);
static tiifhe_Digit digit_mul(size_t idx, tiifhe_Digit a, tiifhe_Digit b);
static tiifhe_Digit digit_neg(size_t idx, tiifhe_Digit a);
static tiifhe_Digit digit_sub(size_t idx, tiifhe_Digit a, tiifhe_Digit b);

/* https://eprint.iacr.org/2018/039 */
static int64_t montgomery_deinit(int64_t a, int64_t mod, int64_t inv);
static int64_t montgomery_init(int64_t a, int64_t mod, int64_t inv, int64_t pow);
static int64_t montgomery_mulred(int64_t a, int64_t b, int64_t mod, int64_t inv);
static int64_t montgomery_powred(int64_t a, size_t exp, int64_t mod, int64_t inv, int64_t pow);
static int64_t montgomery_precomp_inv(int64_t mod);
static int64_t montgomery_precomp_pow(int64_t mod);
static int64_t montgomery_reduce(int128_t a, int64_t mod, int64_t inv);

/* https://eprint.iacr.org/2018/039 */
static void    seiler_intt(int64_t *poly, size_t degree, int64_t mod, int64_t root, int64_t barrett, int64_t inv, int64_t pow);
static void    seiler_ntt(int64_t *poly, size_t degree, int64_t mod, int64_t root, int64_t barrett, int64_t inv, int64_t pow);
static int64_t seiler_precomp(size_t degree, int64_t mod, int64_t barrett, int64_t inv, int64_t pow);

static int64_t
barrett_precomp(int64_t mod)
{
	uint128_t precomp;

	precomp = (uint128_t)1 << (63 + BSR64(mod));
	precomp = (precomp + (mod >> 1)) / mod;

	return -(int64_t)precomp;
}

static int64_t
barrett_reduce(int64_t a, int64_t mod, int64_t precomp)
{
	int128_t tmp128;
	int64_t tmp64;

	/* check for arithmetic shift */
	assert(-1 >> 1 == -1);

	tmp128 = a;
	tmp64 = (tmp128 * precomp) >> 64;
	tmp128 = tmp64 >> (BSR64(mod) - 1);
	tmp64 = (int64_t)tmp128 * mod + a;

	return tmp64 + ((tmp64 >> 63) & mod);
}

tiifhe_Digit
digit_add(size_t idx, tiifhe_Digit a, tiifhe_Digit b)
{

	return barrett_reduce(a + b, tiifhe_q[idx].value, tiifhe_q[idx].barrett);
}

tiifhe_Digit
digit_deinit(size_t idx, tiifhe_Digit a)
{

	return montgomery_deinit(a, tiifhe_q[idx].value, tiifhe_q[idx].inv);
}

tiifhe_Digit
digit_init(size_t idx, tiifhe_Digit a)
{

	return montgomery_init(a, tiifhe_q[idx].value, tiifhe_q[idx].inv, tiifhe_q[idx].pow);
}

tiifhe_Digit
digit_inv(size_t idx, tiifhe_Digit a)
{
	int64_t inv;

	inv = tiifhe_util_invmod(digit_deinit(idx, a), tiifhe_q[idx].value);
	return digit_init(idx, inv);
}

tiifhe_Digit
digit_mul(size_t idx, tiifhe_Digit a, tiifhe_Digit b)
{

	return montgomery_mulred(a, b, tiifhe_q[idx].value, tiifhe_q[idx].inv);
}

tiifhe_Digit
digit_neg(size_t idx, tiifhe_Digit a)
{

	(void)idx;

	return -a;
}

tiifhe_Digit
digit_sub(size_t idx, tiifhe_Digit a, tiifhe_Digit b)
{

	return barrett_reduce(a - b, tiifhe_q[idx].value, tiifhe_q[idx].barrett);
}

static int64_t
montgomery_deinit(int64_t a, int64_t mod, int64_t inv)
{
	int64_t tmp;

	/* check for arithmetic shift */
	assert(-1 >> 1 == -1);

	tmp = montgomery_mulred(a, 1, mod, inv);

	return tmp + ((tmp >> 63) & mod);
}

static int64_t
montgomery_init(int64_t a, int64_t mod, int64_t inv, int64_t pow)
{

	return montgomery_mulred(a, pow, mod, inv);
}

static int64_t
montgomery_mulred(int64_t a, int64_t b, int64_t mod, int64_t inv)
{

	return montgomery_reduce((int128_t)a * b, mod, inv);
}

static int64_t
montgomery_powred(int64_t a, size_t exp, int64_t mod, int64_t inv, int64_t pow)
{
	int64_t ret;

	ret = montgomery_init(1, mod, inv, pow);
	while (exp > 0) {
		if (exp & 1)
			ret = montgomery_mulred(a, ret, mod, inv);
		a = montgomery_mulred(a, a, mod, inv);

		exp >>= 1;
	}

	return ret;
}

static int64_t
montgomery_precomp_inv(int64_t mod)
{
	uint64_t inv;
	int i;

	/* https://crypto.stackexchange.com/questions/47493 */
	inv = mod & 1;
	for (i = 0; i < 6; ++i)
		inv *= 2 - mod * inv;

	return inv;
}

static int64_t
montgomery_precomp_pow(int64_t mod)
{
	int64_t pow;
	int i;

	pow = 1;
	for (i = 0; i < 128; ++i)
		pow <<= 1, pow %= mod;

	return pow;
}

static int64_t
montgomery_reduce(int128_t a, int64_t mod, int64_t inv)
{
	int64_t tmp;

	/* check for arithmetic shift */
	assert(-1 >> 1 == -1);

	tmp = (int128_t)(int64_t)a * inv;
	tmp = (a - (int128_t)tmp * mod) >> 64;

	return tmp;
}

static void
seiler_intt(int64_t *poly, size_t degree, int64_t mod, int64_t root, int64_t barrett, int64_t inv, int64_t pow)
{
	size_t dlog, dinv, j, k, l, s;

	dlog = BSR64(degree);
	dinv = montgomery_init(tiifhe_util_invmod(degree, mod), mod, inv, pow);

	k = 0;
	for (l = 1; l < degree; l <<= 1) {
		for (s = 0; s < degree; s = j + l) {
			int64_t r;
			size_t exp;

			exp = 1 + tiifhe_util_bitrev(k++, dlog);
			r = montgomery_powred(root, exp, mod, inv, pow);
			r = montgomery_deinit(r, mod, inv);
			r = montgomery_init(tiifhe_util_invmod(r, mod), mod, inv, pow);

			for (j = s; j < s + l; ++j) {
				int64_t tmp;

				tmp = poly[j];
				poly[j] = barrett_reduce(tmp + poly[j + l], mod, barrett);
				poly[j + l] = barrett_reduce(tmp - poly[j + l], mod, barrett);
				poly[j + l] = montgomery_mulred(r, poly[j + l], mod, inv);
			}
		}
	}

	for (j = 0; j < degree; ++j)
		poly[j] = montgomery_mulred(poly[j], dinv, mod, inv);
}

static void
seiler_ntt(int64_t *poly, size_t degree, int64_t mod, int64_t root, int64_t barrett, int64_t inv, int64_t pow)
{
	size_t dlog, j, k, l, s;

	dlog = BSR64(degree);

	k = 1;
	for (l = degree >> 1; l > 0; l >>= 1) {
		for (s = 0; s < degree; s = j + l) {
			int64_t r;
			size_t exp;

			exp = tiifhe_util_bitrev(k++, dlog);
			r = montgomery_powred(root, exp, mod, inv, pow);

			for (j = s; j < s + l; ++j) {
				int64_t tmp;

				tmp = montgomery_mulred(r, poly[j + l], mod, inv);
				poly[j + l] = barrett_reduce(poly[j] - tmp, mod, barrett);
				poly[j] = barrett_reduce(poly[j] + tmp, mod, barrett);
			}
		}
	}
}

static int64_t
seiler_precomp(size_t degree, int64_t mod, int64_t barrett, int64_t inv, int64_t pow)
{
	int64_t root, cnt, one;
	size_t ord, exp;

	ord = mod - 1;
	exp = ord >> (BSR64(degree) + 1);
	one = montgomery_init(1, mod, inv, pow);

	cnt = one;
	for (;;) {
		cnt = barrett_reduce(cnt + one, mod, barrett);
		if (montgomery_powred(cnt, ord, mod, inv, pow) != one)
			continue;

		root = montgomery_powred(cnt, exp, mod, inv, pow);
		if (montgomery_powred(root, degree, mod, inv, pow) != one)
			break;
	}

	return root;
}

tiifhe_Digit
tiifhe_mod_const(size_t idx, const tiifhe_Digit op)
{

	return digit_init(idx, op);
}

void
tiifhe_mod_deinit(size_t idx)
{

	(void)idx;
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
	int64_t barrett, mod, inv, P, pow;
	size_t i;

	mod = tiifhe_rns[idx];
	barrett = barrett_precomp(mod);
	inv = montgomery_precomp_inv(mod);
	pow = montgomery_precomp_pow(mod);

	tiifhe_q[idx].value = mod;
	tiifhe_q[idx].barrett = barrett;
	tiifhe_q[idx].inv = inv;
	tiifhe_q[idx].pow = pow;
	tiifhe_q[idx].root = seiler_precomp(TIIFHE_N, mod, barrett, inv, pow);
	tiifhe_q[idx].t = montgomery_init(t, mod, inv, pow);
	tiifhe_q[idx].idx = idx;

	for (i = 0; i < TIIFHE_QPLEN; ++i) {
		tiifhe_Digit tmp;

		tmp = digit_init(idx, tiifhe_rns[i]);

		tiifhe_rns_mods[idx][i] = tmp;
		if (idx == i)
			tiifhe_rns_invs[idx][i] = digit_init(idx, 1);
		else
			tiifhe_rns_invs[idx][i] = digit_inv(idx, tmp);
	}

	P = digit_init(idx, 1);
	for (i = TIIFHE_QLEN; i < TIIFHE_QPLEN; ++i)
		P = digit_mul(idx, P, tiifhe_rns_mods[idx][i]);
	tiifhe_q[idx].P = P;
}

tiifhe_Digit
tiifhe_mod_inv(size_t idx, const tiifhe_Digit op)
{

	return digit_inv(idx, op);
}

tiifhe_Digit
tiifhe_mod_mul(size_t idx, const tiifhe_Digit op1, const tiifhe_Digit op2)
{

	return digit_mul(idx, op1, op2);
}

tiifhe_Digit
tiifhe_mod_neg(size_t idx, const tiifhe_Digit op)
{

	return digit_neg(idx, op);
}

void
tiifhe_poly_add(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Poly *op2)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j)
		rop->value[j] = digit_add(idx, op1->value[j], op2->value[j]);
}

void
tiifhe_poly_addmul(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Poly *op2, const tiifhe_Poly *op3, const tiifhe_Poly *op4)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j) {
		int64_t add12, add34;

		add12 = digit_add(idx, op1->value[j], op2->value[j]);
		add34 = digit_add(idx, op3->value[j], op4->value[j]);
		rop->value[j] = digit_mul(idx, add12, add34);
	}
}

void
tiifhe_poly_cmod(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *p)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j) {
		rop->value[j] = p->value[j];
		if (p->value[j] > tiifhe_q[idx].value / 2)
			rop->value[j] -= tiifhe_q[idx].value;
		if (p->value[j] < -tiifhe_q[idx].value / 2)
			rop->value[j] += tiifhe_q[idx].value;
	}
}

void
tiifhe_poly_deinit(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *p)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j)
		rop->value[j] = digit_deinit(idx, p->value[j]);
}

void
tiifhe_poly_ext(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly **p, const size_t *base, const size_t *drops, size_t len)
{
	tiifhe_Digit inv[TIIFHE_QPLEN], div[TIIFHE_QPLEN];
	size_t i, j;

	for (i = 0; i < len; ++i) {
		size_t ii;

		inv[i] = digit_init(base[i], 1);
		for (ii = 0; ii < len; ++ii) {
			if (drops[ii] != 0)
				continue;
			inv[i] = digit_mul(base[i], inv[i], tiifhe_rns_invs[base[i]][base[ii]]);
		}

		div[i] = tiifhe_rns_invs[idx][base[i]];
		for (ii = 0; ii < len; ++ii) {
			if (idx == base[i] && idx == base[ii])
				continue;
			if (drops[ii] != 0)
				continue;
			div[i] = digit_mul(idx, div[i], tiifhe_rns_mods[idx][base[ii]]);
		}
	}

	for (j = 0; j < TIIFHE_N; ++j) {
		tiifhe_Digit sum;

		sum = 0;
		for (i = 0; i < len; ++i) {
			tiifhe_Digit tmp;

			if (drops[i] != 0)
				continue;

			tmp = digit_mul(base[i], p[i]->value[j], inv[i]);
			tmp = digit_deinit(base[i], tmp);
			tmp = digit_init(idx, tmp);
			tmp = digit_mul(idx, tmp, div[i]);
			sum = digit_add(idx, sum, tmp);
		}
		rop->value[j] = sum;
	}
}

void
tiifhe_poly_init(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *p)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j)
		rop->value[j] = digit_init(idx, p->value[j]);
}

void
tiifhe_poly_init_error(size_t idx, tiifhe_Poly *p, const tiifhe_Poly *e)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j) {
		tiifhe_Digit init = digit_init(idx, e->value[j]);
		p->value[j] = digit_mul(idx, init, tiifhe_q[idx].t);
	}
	tiifhe_poly_ntt(idx, p);
}

void
tiifhe_poly_intt(size_t idx, tiifhe_Poly *p)
{

	seiler_intt(p->value, TIIFHE_N, tiifhe_q[idx].value, tiifhe_q[idx].root, tiifhe_q[idx].barrett, tiifhe_q[idx].inv, tiifhe_q[idx].pow);
}

void
tiifhe_poly_mod_mpz(size_t idx, tiifhe_Poly *rop, const mpz_t *p)
{
	mpz_t tmp;
	size_t j;

	mpz_init(tmp);

	for (j = 0; j < TIIFHE_N; ++j) {
		mpz_mod_ui(tmp, p[j], tiifhe_q[idx].value);
		rop->value[j] = digit_init(idx, mpz_get_ui(tmp));
	}

	mpz_clear(tmp);
}

void
tiifhe_poly_mod_u64(size_t idx, tiifhe_Poly *rop, const uint64_t *p)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j)
		rop->value[j] = digit_init(idx, p[j] % tiifhe_q[idx].value);
}

void
tiifhe_poly_mul(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Poly *op2)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j)
		rop->value[j] = digit_mul(idx, op1->value[j], op2->value[j]);
}

void
tiifhe_poly_muladd(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Poly *op2, const tiifhe_Poly *op3)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j)
		rop->value[j] = digit_add(idx, digit_mul(idx, op1->value[j], op2->value[j]), op3->value[j]);
}

void
tiifhe_poly_mulc(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Digit op2)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j)
		rop->value[j] = digit_mul(idx, op1->value[j], op2);
}

void
tiifhe_poly_mulcadd(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Digit op2, const tiifhe_Poly *op3)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j)
		rop->value[j] = digit_add(idx, digit_mul(idx, op1->value[j], op2), op3->value[j]);
}

void
tiifhe_poly_neg(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *p)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j)
		rop->value[j] = digit_neg(idx, p->value[j]);
}

void
tiifhe_poly_ntt(size_t idx, tiifhe_Poly *p)
{

	seiler_ntt(p->value, TIIFHE_N, tiifhe_q[idx].value, tiifhe_q[idx].root, tiifhe_q[idx].barrett, tiifhe_q[idx].inv, tiifhe_q[idx].pow);
}

void
tiifhe_poly_rot_intt(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *p, size_t steps)
{
	size_t step, i, j;

	assert(rop != p);
	(void)idx;

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
			tmp = -tmp;
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

	tiifhe_random_cbd21_i64(e->value, TIIFHE_N, seed);
}

void
tiifhe_poly_sample_secret(tiifhe_Poly *s, tiifhe_Seed *seed)
{

	tiifhe_random_cbd1_i64(s->value, TIIFHE_N, seed);
}

void
tiifhe_poly_sample_uniform(size_t idx, tiifhe_Poly *p, tiifhe_Seed *seed)
{

	tiifhe_random_uniform_i64(p->value, TIIFHE_N, tiifhe_q[idx].value, seed);
}

void
tiifhe_poly_sub(size_t idx, tiifhe_Poly *rop, const tiifhe_Poly *op1, const tiifhe_Poly *op2)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j)
		rop->value[j] = digit_sub(idx, op1->value[j], op2->value[j]);
}
