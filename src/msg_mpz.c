#include "msg_mpz.h"

tiifhe_MessageMod tiifhe_t;

#if defined(__x86_64__)
	#include <x86intrin.h>
	#define BSR64(x) __bsrq(x)
#else
	#define BSR64(x) tiifhe_util_msb64(x)
#endif

/* https://en.wikipedia.org/wiki/Chinese_remainder_theorem */
static void mpz_cmod(mpz_t rop, const mpz_t op, const mpz_t mod);
static void mpz_crt(mpz_t *poly, size_t degree, const uint64_t **rns, const uint64_t *mods, size_t len);
static void mpz_precomp(mpz_t root, size_t degree, const mpz_t mod);
static void mpz_ntt(mpz_t *poly, size_t degree, const mpz_t mod, const mpz_t root);
static void mpz_intt(mpz_t *poly, size_t degree, const mpz_t mod, const mpz_t root);

static void
mpz_cmod(mpz_t rop, const mpz_t op, const mpz_t mod)
{
	mpz_t bound;

	mpz_init_set(bound, mod);
	mpz_tdiv_q_2exp(bound, mod, 1);

	mpz_mod(rop, op, mod);
	if (mpz_cmp_ui(op, 0) < 0)
		mpz_add(rop, rop, mod);
	if (mpz_cmp(rop, bound) >= 0)
		mpz_sub(rop, rop, mod);

	mpz_clear(bound);
}

static void
mpz_crt(mpz_t *poly, size_t degree, const uint64_t **rns, const uint64_t *mods, size_t len)
{
	mpz_t mod, inv, tmp;
	size_t i, j;

	mpz_inits(mod, inv, tmp, 0);

	mpz_set_ui(mod, mods[0]);
	for (j = 0; j < degree; ++j)
		mpz_set_ui(poly[j], rns[0][j]);

	for (i = 1; i < len; ++i) {
		mpz_set_ui(inv, mods[i]);
		mpz_invert(inv, inv, mod);
		mpz_mul_ui(mod, mod, mods[i]);

		for (j = 0; j < degree; ++j) {
			mpz_sub_ui(poly[j], poly[j], rns[i][j]);
			mpz_mul_ui(poly[j], poly[j], mods[i]);
			mpz_mod(poly[j], poly[j], mod);

			mpz_mul(poly[j], poly[j], inv);
			mpz_add_ui(poly[j], poly[j], rns[i][j]);
			mpz_mod(poly[j], poly[j], mod);
		}
	}

	mpz_clears(mod, inv, tmp, 0);
}

static void
mpz_precomp(mpz_t root, size_t degree, const mpz_t mod)
{
	mpz_t ord, exp, pow, tmp;

	mpz_inits(ord, exp, pow, tmp, 0);

	mpz_sub_ui(ord, mod, 1);
	mpz_tdiv_q_2exp(exp, ord, BSR64(degree) + 1);
	mpz_set_ui(tmp, 1);

	for (;;) {
		mpz_add_ui(tmp, tmp, 1);
		mpz_mod(tmp, tmp, mod);

		mpz_powm(pow, tmp, ord, mod);
		if (mpz_cmp_ui(pow, 1) != 0)
			continue;

		mpz_powm(root, tmp, exp, mod);
		mpz_powm_ui(pow, root, degree, mod);
		if (mpz_cmp_ui(pow, 1) != 0)
			break;
	}

	mpz_clears(ord, exp, pow, tmp, 0);
}

static void
mpz_ntt(mpz_t *poly, size_t degree, const mpz_t mod, const mpz_t root)
{
	mpz_t r, tmp;
	size_t dlog, exp, j, k, l, s;

	mpz_inits(r, tmp, 0);
	dlog = BSR64(degree);

	k = 1;
	for (l = degree >> 1; l > 0; l >>= 1) {
		for (s = 0; s < degree; s = j + l) {
			exp = tiifhe_util_bitrev(k++, dlog);
			mpz_powm_ui(r, root, exp, mod);

			for (j = s; j < s + l; ++j) {
				mpz_mul(tmp, r, poly[j + l]);
				mpz_mod(tmp, tmp, mod);

				mpz_sub(poly[j + l], poly[j], tmp);
				mpz_mod(poly[j + l], poly[j + l], mod);

				mpz_add(poly[j], poly[j], tmp);
				mpz_mod(poly[j], poly[j], mod);
			}
		}
	}

	mpz_clears(r, tmp, 0);
}

static void
mpz_intt(mpz_t *poly, size_t degree, const mpz_t mod, const mpz_t root)
{
	mpz_t dinv, r, tmp;
	size_t dlog, exp, j, k, l, s;

	mpz_inits(dinv, r, tmp, 0);
	mpz_set_ui(dinv, degree);
	mpz_invert(dinv, dinv, mod);
	dlog = BSR64(degree);

	k = 0;
	for (l = 1; l < degree; l <<= 1) {
		for (s = 0; s < degree; s = j + l) {
			exp = 1 + tiifhe_util_bitrev(k++, dlog);
			mpz_powm_ui(r, root, exp, mod);
			mpz_invert(r, r, mod);

			for (j = s; j < s + l; ++j) {
				mpz_set(tmp, poly[j]);

				mpz_add(poly[j], poly[j + l], tmp);
				mpz_mod(poly[j], poly[j], mod);

				mpz_sub(poly[j + l], tmp, poly[j + l]);
				mpz_mod(poly[j + l], poly[j + l], mod);

				mpz_mul(poly[j + l], poly[j + l], r);
				mpz_mod(poly[j + l], poly[j + l], mod);
			}
		}
	}

	for (j = 0; j < degree; ++j) {
		mpz_mul(poly[j], poly[j], dinv);
		mpz_mod(poly[j], poly[j], mod);
	}

	mpz_clears(dinv, r, tmp, 0);
}


tiifhe_Message *
tiifhe_msg_alloc(size_t len)
{
	tiifhe_Message *m;
	size_t i, j;

	m = tiifhe_util_alloc(len, sizeof *m);
	for (i = 0; i < len; ++i)
		for (j = 0; j < TIIFHE_N; ++j)
			mpz_init_set_ui(m[i].value[j], 0);

	return m;
}

void
tiifhe_msg_dealloc(tiifhe_Message *m, size_t len)
{
	size_t i, j;

	for (i = 0; i < len; ++i)
		for (j = 0; j < TIIFHE_N; ++j)
			mpz_clear(m[i].value[j]);
	free(m);
}

void
tiifhe_msg_addc(tiifhe_Message *rop, const tiifhe_Message *m, uint64_t c)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j) {
		mpz_add_ui(rop->value[j], m->value[j], c);
		mpz_mod(rop->value[j], rop->value[j], tiifhe_t.value);
	}
}

void
tiifhe_msg_cmod(tiifhe_Message *rop, const tiifhe_Message *m, const mpz_t mod)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j)
		mpz_cmod(rop->value[j], m->value[j], mod);
}

void
tiifhe_msg_copy(tiifhe_Message *rop, const tiifhe_Message *m)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j)
		mpz_set(rop->value[j], m->value[j]);
}

void
tiifhe_msg_crt_i64(tiifhe_Message *m, const int64_t **rns, const int64_t *mods, size_t len)
{
	const uint64_t **rns2;
	uint64_t *mods2;
	size_t i;

	mods2 = tiifhe_util_alloc(len, sizeof *mods2);
	rns2 = tiifhe_util_alloc(len, sizeof *rns2);

	for (i = 0; i < len; ++i) {
		mods2[i] = mods[i];
		rns2[i] = (uint64_t *)rns[i];
	}
	tiifhe_msg_crt_u64(m, rns2, mods2, len);

	tiifhe_util_dealloc(mods2);
	tiifhe_util_dealloc(rns2);
}

void
tiifhe_msg_crt_u64(tiifhe_Message *m, const uint64_t **rns, const uint64_t *mods, size_t len)
{

	mpz_crt(m->value, TIIFHE_N, rns, mods, len);
}

void
tiifhe_msg_fprintdiv(FILE *f, const tiifhe_Message *m, const mpz_t mod)
{
	mpz_t tmp;
	size_t j;

	mpz_init(tmp);

	for (j = 0; j < TIIFHE_N - 1; ++j) {
		mpz_tdiv_q(tmp, m->value[j], mod);
		gmp_fprintf(f, "%Zd, ", tmp);
	}
	mpz_tdiv_q(tmp, m->value[j], mod);
	gmp_fprintf(f, "%Zd\n", tmp);

	mpz_clear(tmp);
}

void
tiifhe_msg_mulc(tiifhe_Message *rop, const tiifhe_Message *m, uint64_t c)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j) {
		mpz_mul_ui(rop->value[j], m->value[j], c);
		mpz_mod(rop->value[j], rop->value[j], tiifhe_t.value);
	}
}

void
tiifhe_msg_mulc_inv(tiifhe_Message *rop, const tiifhe_Message *m, uint64_t c)
{
	mpz_t tmp;
	size_t j;

	mpz_init_set_ui(tmp, c);
	mpz_invert(tmp, tmp, tiifhe_t.value);

	for (j = 0; j < TIIFHE_N; ++j) {
		mpz_mul(rop->value[j], m->value[j], tmp);
		mpz_mod(rop->value[j], rop->value[j], tiifhe_t.value);
	}

	mpz_clear(tmp);
}

void
tiifhe_msg_pack(tiifhe_Message *rop, const tiifhe_Message *m)
{
	tiifhe_Message *cpy;
	size_t dlog, idx, j;

	cpy = tiifhe_msg_alloc(1);

	/* reorder slots */
	idx = 1;
	dlog = BSR64(TIIFHE_N);
	for (j = 0; j < TIIFHE_N / 2; ++j) {
		size_t idx1, idx2;

		idx1 = tiifhe_util_bitrev((idx - 1) / 2, dlog);
		idx2 = tiifhe_util_bitrev((2 * TIIFHE_N - idx - 1) / 2, dlog);

		mpz_set(cpy->value[idx1], m->value[j]);
		mpz_set(cpy->value[idx2], m->value[j | TIIFHE_N / 2]);

		idx *= 3;
		idx %= 2 * TIIFHE_N;
	}

	/* actual packing */
	mpz_intt(cpy->value, TIIFHE_N, tiifhe_t.value, tiifhe_t.root);
	for (j = 0; j < TIIFHE_N; ++j)
		mpz_set(rop->value[j], cpy->value[j]);

	tiifhe_msg_dealloc(cpy, 1);
}

void
tiifhe_msg_pmod(tiifhe_Message *rop, const tiifhe_Message *m, const mpz_t mod)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j)
		mpz_mod(rop->value[j], m->value[j], mod);
}

void
tiifhe_msg_unpack(tiifhe_Message *rop, const tiifhe_Message *m)
{
	tiifhe_Message *cpy;
	size_t dlog, idx, j;

	cpy = tiifhe_msg_alloc(1);
	tiifhe_msg_copy(cpy, m);

	/* actual unpacking */
	mpz_ntt(cpy->value, TIIFHE_N, tiifhe_t.value, tiifhe_t.root);

	/* reorder slots */
	idx = 1;
	dlog = BSR64(TIIFHE_N);
	for (j = 0; j < TIIFHE_N / 2; ++j) {
		size_t idx1, idx2;

		idx1 = tiifhe_util_bitrev((idx - 1) / 2, dlog);
		idx2 = tiifhe_util_bitrev((2 * TIIFHE_N - idx - 1) / 2, dlog);

		mpz_set(rop->value[j], cpy->value[idx1]);
		mpz_set(rop->value[j | TIIFHE_N / 2], cpy->value[idx2]);

		idx *= 3;
		idx %= 2 * TIIFHE_N;
	}

	tiifhe_msg_dealloc(cpy, 1);
}

void
tiifhe_msgmod_deinit(void)
{

	mpz_clears(tiifhe_t.value, tiifhe_t.root, 0);
}

void
tiifhe_msgmod_init(void)
{

	mpz_inits(tiifhe_t.value, tiifhe_t.root, 0);
	mpz_set_str(tiifhe_t.value, TIIFHE_T, 0);
	mpz_precomp(tiifhe_t.root, TIIFHE_N, tiifhe_t.value);
}
