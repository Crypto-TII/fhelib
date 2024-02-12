#include "bgv.h"

size_t tiifhe_decomp_idx[TIIFHE_OMEGA][TIIFHE_KMAX];
size_t tiifhe_decomp_len[TIIFHE_OMEGA];

#define MAX(a, b) ((a) > (b) ? (a) : (b))

static void ciphertext_intt(size_t idx, tiifhe_Ciphertext *ct);
static void ciphertext_ntt(size_t idx, tiifhe_Ciphertext *ct);
static void ciphertext_switch_ntt(size_t idx, tiifhe_CiphertextSwitch *csw);
static void delta_ext(size_t idx, tiifhe_Delta *rop, tiifhe_Delta *delta, const size_t *base, size_t len);

static void
ciphertext_intt(size_t idx, tiifhe_Ciphertext *ct)
{

	if (ct->drop != 0)
		return;
	if (ct->ntt == 0)
		return;
	ct->ntt = 0;

	tiifhe_poly_intt(idx, &ct->poly[0]);
	tiifhe_poly_intt(idx, &ct->poly[1]);

	if (ct->degree == 2)
		return;

	tiifhe_poly_intt(idx, &ct->poly[2]);
}

static void
ciphertext_ntt(size_t idx, tiifhe_Ciphertext *ct)
{

	if (ct->drop != 0)
		return;
	if (ct->ntt == 1)
		return;
	ct->ntt = 1;

	tiifhe_poly_ntt(idx, &ct->poly[0]);
	tiifhe_poly_ntt(idx, &ct->poly[1]);

	if (ct->degree == 2)
		return;

	tiifhe_poly_ntt(idx, &ct->poly[2]);
}

static void
ciphertext_switch_ntt(size_t idx, tiifhe_CiphertextSwitch *csw)
{
	size_t k;

	if (csw->ct.drop != 0)
		return;
	if (csw->ct.ntt == 1)
		return;

	ciphertext_ntt(idx, &csw->ct);
	for (k = 0; k < TIIFHE_OMEGA; ++k)
		tiifhe_poly_ntt(idx, &csw->ext[k]);
}

static void
delta_ext(size_t idx, tiifhe_Delta *rop, tiifhe_Delta *delta, const size_t *base, size_t len)
{
	const tiifhe_Poly *p[TIIFHE_QPLEN];
	size_t drops[TIIFHE_QPLEN], i;
	tiifhe_Digit tmp;

	ciphertext_intt(idx, &delta->ct);
	for (i = 0; i < len; ++i)
		drops[i] = delta[i].ct.drop;

	rop->ct.degree = delta->ct.degree;
	rop->ct.drop = 0;
	rop->ct.ntt = 0;
	rop->inv = tiifhe_mod_const(idx, 1);
	for (i = 0; i < len; ++i) {
		if (drops[i] != 0)
			continue;

		tmp = tiifhe_mod_inv(idx, tiifhe_mod_const(idx, tiifhe_q[base[i]].value));
		rop->inv = tiifhe_mod_mul(idx, rop->inv, tmp);
	}

	for (i = 0; i < len; ++i)
		p[i] = &delta[i].ct.poly[0];
	tiifhe_poly_ext(idx, &rop->ct.poly[0], p, base, drops, len);

	for (i = 0; i < len; ++i)
		p[i] = &delta[i].ct.poly[1];
	tiifhe_poly_ext(idx, &rop->ct.poly[1], p, base, drops, len);

	if (rop->ct.degree == 2)
		return;

	for (i = 0; i < len; ++i)
		p[i] = &delta[i].ct.poly[2];
	tiifhe_poly_ext(idx, &rop->ct.poly[2], p, base, drops, len);
}

tiifhe_Ciphertext *
tiifhe_bgv_alloc_ct(void)
{

	return tiifhe_util_alloc(TIIFHE_QPLEN, sizeof(tiifhe_Ciphertext));
}

tiifhe_KeyPublic *
tiifhe_bgv_alloc_pk(void)
{

	return tiifhe_util_alloc(TIIFHE_QLEN, sizeof(tiifhe_KeyPublic));
}

tiifhe_KeySecret *
tiifhe_bgv_alloc_sk(void)
{

	return tiifhe_util_alloc(1, sizeof(tiifhe_KeySecret));
}

tiifhe_KeySwitch *
tiifhe_bgv_alloc_ksw(void)
{

	return tiifhe_util_alloc(TIIFHE_QPLEN, sizeof(tiifhe_KeySwitch));
}

void
tiifhe_bgv_dealloc(void *ptr)
{

	tiifhe_util_dealloc(ptr);
}

void
tiifhe_bgv_dealloc_ct(tiifhe_Ciphertext *ct)
{

	tiifhe_util_dealloc(ct);
}

void
tiifhe_bgv_dealloc_pk(tiifhe_KeyPublic *pk)
{

	tiifhe_util_dealloc(pk);
}

void
tiifhe_bgv_dealloc_sk(tiifhe_KeySecret *sk)
{

	tiifhe_util_dealloc(sk);
}

void
tiifhe_bgv_dealloc_ksw(tiifhe_KeySwitch *ksw)
{

	tiifhe_util_dealloc(ksw);
}

void
tiifhe_bgv_deinit(void)
{
	size_t i;

	tiifhe_msgmod_deinit();
	for (i = 0; i < TIIFHE_QPLEN; ++i)
		tiifhe_mod_deinit(i);
}

void
tiifhe_bgv_init(void)
{
	mpz_t tmp;
	size_t div, rest, len, i, k;

	mpz_init(tmp);

	tiifhe_msgmod_init();

	for (i = 0; i < TIIFHE_QPLEN; ++i)
		tiifhe_mod_init(i, tiifhe_t.value);

	/* initialize chunks */
	div = TIIFHE_QLEN / TIIFHE_OMEGA;
	rest = TIIFHE_QLEN % TIIFHE_OMEGA;

	for (k = 0; k < TIIFHE_OMEGA; ++k) {
		len = div;
		if (rest > 0)
			++len, --rest;
		tiifhe_decomp_len[k] = len;
	}

	for (i = 0; i < TIIFHE_QLEN; ++i)
		tiifhe_decomp_idx[i % TIIFHE_OMEGA][i / TIIFHE_OMEGA] = i;

	mpz_clear(tmp);
}

void
tiifhe_bgv_add(tiifhe_Ciphertext *rop, tiifhe_Ciphertext *op1, tiifhe_Ciphertext *op2)
{
	size_t i;

	for (i = 0; i < TIIFHE_QLEN; ++i)
		tiifhe_bgv_idx_add(i, &rop[i], &op1[i], &op2[i]);
}

void
tiifhe_bgv_addc(tiifhe_Ciphertext *rop, tiifhe_Ciphertext *ct, const tiifhe_Message *m)
{
	tiifhe_Message *cpy;
	tiifhe_Encoding p;
	size_t i, it;

	cpy = tiifhe_msg_alloc(1);
	tiifhe_msg_copy(cpy, m);

	for (i = 0; i < TIIFHE_QLEN; ++i) {
		for (it = 0; it < ct[i].drop; ++it)
			tiifhe_msg_mulc_inv(cpy, cpy, tiifhe_q[i].value);
	}

	for (i = 0; i < TIIFHE_QLEN; ++i) {
		tiifhe_bgv_idx_encode(i, &p, cpy);
		tiifhe_bgv_idx_addc(i, &rop[i], &ct[i], &p);
	}

	tiifhe_msg_dealloc(cpy, 1);
}

void
tiifhe_bgv_copy(tiifhe_Ciphertext *rop, const tiifhe_Ciphertext *ct)
{

	memcpy(rop, ct, TIIFHE_QPLEN * sizeof *rop);
}

void
tiifhe_bgv_decrypt(tiifhe_Message *m, const tiifhe_KeySecret *sk, tiifhe_Ciphertext *ct)
{
	tiifhe_Encoding *p;
	size_t i;

	p = tiifhe_util_alloc(TIIFHE_QLEN, sizeof *p);

	for (i = 0; i < TIIFHE_QLEN; ++i)
		tiifhe_bgv_idx_decrypt(i, &p[i], sk, &ct[i]);

	tiifhe_bgv_decode(m, p);
	tiifhe_util_dealloc(p);
}

void
tiifhe_bgv_drop(tiifhe_Ciphertext *ct, const tiifhe_Ciphertext *cmp)
{
	tiifhe_Message *c;
	size_t drops[TIIFHE_QLEN];
	size_t drop, i, it;

	c = tiifhe_msg_alloc(1);
	tiifhe_msg_addc(c, c, 1);

	drop = 0;
	for (i = 0; i < TIIFHE_QLEN; ++i) {
		if (ct[i].drop > cmp[i].drop) {
			fprintf(stdout, "fhelib: drop target has a higher level at idx %lu\n", i);
			exit(1);
		}

		drops[i] = cmp[i].drop - ct[i].drop;
		drop += drops[i];

		for (it = 1; it < drops[i]; ++it)
			tiifhe_msg_mulc_inv(c, c, tiifhe_q[i].value);
	}

	if (drop == 0)
		return;

	tiifhe_msg_pack(c, c);
	tiifhe_bgv_mulc(ct, ct, c);
	for (i = 0; i < TIIFHE_QLEN; ++i)
		if (drops[i])
			tiifhe_bgv_modswitch(ct);

	tiifhe_msg_dealloc(c, 1);
}

void
tiifhe_bgv_encrypt(tiifhe_Ciphertext *ct, tiifhe_KeyPublic *pk, const tiifhe_Message *m)
{
	tiifhe_Seed seed;
	tiifhe_Encoding p;
	size_t i;

	tiifhe_random_seed(&seed);
	for (i = 0; i < TIIFHE_QLEN; ++i) {
		tiifhe_bgv_idx_encode(i, &p, m);
		tiifhe_bgv_idx_encrypt(i, &ct[i], &pk[i], &p, &seed);
	}
}

void
tiifhe_bgv_keygen_public(tiifhe_KeyPublic *pk, const tiifhe_KeySecret *sk)
{
	tiifhe_Seed seed;
	size_t i;

	tiifhe_random_seed(&seed);
	for (i = 0; i < TIIFHE_QLEN; ++i)
		tiifhe_bgv_idx_keygen_public(i, &pk[i], sk, &seed);
}

void
tiifhe_bgv_keygen_secret(tiifhe_KeySecret *sk)
{
	tiifhe_Seed seed;
	size_t i;
	(void)i;

	tiifhe_random_seed(&seed);
	#if TIIFHE_EXPAND_SEED
	for (i = 0; i < TIIFHE_QPLEN; ++i) {
		tiifhe_poly_sample_secret(&sk->poly[i], &seed);
		tiifhe_poly_init(i, &sk->poly[i], &sk->poly[i]);
		tiifhe_poly_ntt(i, &sk->poly[i]);
		seed.init = 1;
	}
	#else
	tiifhe_poly_sample_secret(&sk->poly, &seed);
	#endif /* TIIFHE_EXPAND_SEED */
}

void
tiifhe_bgv_keygen_switch_mul(tiifhe_KeySwitch *ksw, const tiifhe_KeySecret *sk)
{
	tiifhe_Seed seed[TIIFHE_OMEGA];
	size_t i, k;

	for (k = 0; k < TIIFHE_OMEGA; ++k)
		tiifhe_random_seed(&seed[k]);
	for (i = 0; i < TIIFHE_QPLEN; ++i)
		tiifhe_bgv_idx_keygen_switch_mul(i, &ksw[i], sk, seed);
}

void
tiifhe_bgv_keygen_switch_rot(tiifhe_KeySwitch *ksw, const tiifhe_KeySecret *sk, size_t steps)
{
	tiifhe_Seed seed[TIIFHE_OMEGA];
	size_t i, k;

	for (k = 0; k < TIIFHE_OMEGA; ++k)
		tiifhe_random_seed(&seed[k]);
	for (i = 0; i < TIIFHE_QPLEN; ++i)
		tiifhe_bgv_idx_keygen_switch_rot(i, &ksw[i], sk, steps, seed);
}

void
tiifhe_bgv_keyswitch_mul(tiifhe_Ciphertext *ct, tiifhe_KeySwitch *ksw)
{

	tiifhe_bgv_keyswitch(ct, ksw, 2);
}

void
tiifhe_bgv_keyswitch_rot(tiifhe_Ciphertext *ct, tiifhe_KeySwitch *ksw)
{

	tiifhe_bgv_keyswitch(ct, ksw, 1);
}

void
tiifhe_bgv_modswitch(tiifhe_Ciphertext *ct)
{
	tiifhe_Delta delta;
	size_t idx, i;

	for (i = 0; i < TIIFHE_QLEN; ++i)
		if (ct[i].drop == 0)
			break;

	idx = i;
	if (idx == TIIFHE_QLEN) {
		fprintf(stdout, "fhelib: no modulus left to switch\n");
		exit(1);
	}

	tiifhe_bgv_idx_modswitch_delta(idx, &delta, &ct[idx]);
	for (i = 0; i < TIIFHE_QLEN; ++i)
		tiifhe_bgv_idx_modswitch_ext(i, &ct[i], &delta, 1);
}

void
tiifhe_bgv_mul(tiifhe_Ciphertext *rop, tiifhe_Ciphertext *op1, tiifhe_Ciphertext *op2)
{
	size_t i;

	for (i = 0; i < TIIFHE_QLEN; ++i)
		tiifhe_bgv_idx_mul(i, &rop[i], &op1[i], &op2[i]);
}

void
tiifhe_bgv_mulc(tiifhe_Ciphertext *rop, tiifhe_Ciphertext *ct, const tiifhe_Message *m)
{
	tiifhe_Encoding p;
	size_t i;

	for (i = 0; i < TIIFHE_QLEN; ++i) {
		tiifhe_bgv_idx_encode(i, &p, m);
		tiifhe_bgv_idx_mulc(i, &rop[i], &ct[i], &p);
	}
}

void
tiifhe_bgv_rot(tiifhe_Ciphertext *rop, tiifhe_Ciphertext *ct, size_t steps)
{
	tiifhe_Ciphertext *cpy;
	size_t i;

	cpy = tiifhe_util_alloc(TIIFHE_QLEN, sizeof *cpy);
	memcpy(cpy, ct, TIIFHE_QLEN * sizeof *cpy);

	for (i = 0; i < TIIFHE_QLEN; ++i)
		tiifhe_bgv_idx_rot(i, &rop[i], &cpy[i], steps);

	tiifhe_util_dealloc(cpy);
}

void
tiifhe_bgv_keyswitch(tiifhe_Ciphertext *ct, tiifhe_KeySwitch *ksw, size_t dsw)
{
	tiifhe_CiphertextSwitch *csw;
	tiifhe_Delta *delta;
	size_t i;

	csw = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *csw);
	delta = tiifhe_util_alloc(TIIFHE_PLEN, sizeof *delta);

	tiifhe_bgv_keyswitch_ext(csw, ct, dsw);
	for (i = 0; i < TIIFHE_QPLEN; ++i)
		tiifhe_bgv_idx_keyswitch_dot(i, &ct[i], &csw[i], &ksw[i]);
	for (i = TIIFHE_QLEN; i < TIIFHE_QPLEN; ++i)
		tiifhe_bgv_idx_keyswitch_delta(i, &delta[i - TIIFHE_QLEN], &ct[i]);
	for (i = 0; i < TIIFHE_QLEN; ++i)
		tiifhe_bgv_idx_keyswitch_switch(i, &ct[i], delta);

	tiifhe_util_dealloc(delta);
	tiifhe_util_dealloc(csw);
}

void
tiifhe_bgv_decode(tiifhe_Message *m, const tiifhe_Encoding *p)
{
	tiifhe_Poly *cpy;
	const tiifhe_Digit *rns[TIIFHE_QLEN];
	tiifhe_Digit mods[TIIFHE_QLEN];
	mpz_t mod;

	size_t drops[TIIFHE_QLEN];
	size_t len, i, ii;

	cpy = tiifhe_util_alloc(TIIFHE_QLEN, sizeof *cpy);
	mpz_init_set_ui(mod, 1);

	len = 0;
	for (i = 0; i < TIIFHE_QLEN; ++i) {
		drops[i] = p[i].drop;
		if (drops[i] != 0)
			continue;

		memcpy(&cpy[i], &p[i], sizeof cpy[i]);
		tiifhe_poly_intt(i, &cpy[i]);
		tiifhe_poly_deinit(i, &cpy[i], &cpy[i]);

		mpz_mul_ui(mod, mod, tiifhe_q[i].value);
		mods[len] = tiifhe_q[i].value;
		rns[len] = cpy[i].value;
		++len;
	}
	tiifhe_msg_crt(m, rns, mods, len);
	tiifhe_msg_cmod(m, m, mod);

	#if TIIFHE_LOG_ERROR
	do {
		static FILE *f = 0;

		if (f == 0) {
			f = fopen("error.csv", "w+");
			if (f == 0) {
				perror("fhelib: failed to open error.csv");
				exit(1);
			}
		}

		tiifhe_msg_fprintdiv(f, m, tiifhe_t.value);
	} while (0);
	#endif /* TIIFHE_LOG_ERROR */

	for (i = 0; i < TIIFHE_QLEN; ++i)
		for (ii = 0; ii < drops[i]; ++ii)
			tiifhe_msg_mulc(m, m, tiifhe_q[i].value);

	tiifhe_msg_pmod(m, m, tiifhe_t.value);

	tiifhe_util_dealloc(cpy);
	mpz_clear(mod);
}

void
tiifhe_bgv_keyswitch_ext(tiifhe_CiphertextSwitch *csw, tiifhe_Ciphertext *ct, size_t poly)
{
	const tiifhe_Poly *p[TIIFHE_KMAX];
	size_t drops[TIIFHE_QLEN];
	size_t idx, len, i, k;

	for (idx = 0; idx < TIIFHE_QLEN; ++idx) {
		ciphertext_intt(idx, &ct[idx]);
		memcpy(&csw[idx].ct, &ct[idx], sizeof csw[idx].ct);
	}

	for (idx = 0; idx < TIIFHE_QPLEN; ++idx) {
		if (csw[idx].ct.drop != 0)
			continue;

		for (k = 0; k < TIIFHE_OMEGA; ++k) {
			len = tiifhe_decomp_len[k];
			for (i = 0; i < len; ++i) {
				p[i] = &ct[tiifhe_decomp_idx[k][i]].poly[poly];
				drops[i] = ct[tiifhe_decomp_idx[k][i]].drop;
			}
			tiifhe_poly_ext(idx, &csw[idx].ext[k], p, tiifhe_decomp_idx[k], drops, len);
		}
		csw[idx].poly = poly;
	}
}

void
tiifhe_bgv_idx_add(size_t idx, tiifhe_Ciphertext *rop, tiifhe_Ciphertext *op1, tiifhe_Ciphertext *op2)
{

	if (op1->drop != op2->drop) {
		fprintf(stdout, "fhelib: trying to add ciphertexts with different scaling factors at idx %lu (%lu, %lu)\n", idx, op1->drop, op2->drop);
		exit(1);
	}

	if (op1->ntt != op2->ntt) {
		ciphertext_ntt(idx, op1);
		ciphertext_ntt(idx, op2);
	}

	rop->degree = MAX(op1->degree, op2->degree);
	rop->drop = op1->drop;
	rop->ntt = op1->ntt;

	if (rop->drop != 0)
		return;

	tiifhe_poly_add(idx, &rop->poly[0], &op1->poly[0], &op2->poly[0]);
	tiifhe_poly_add(idx, &rop->poly[1], &op1->poly[1], &op2->poly[1]);

	if (rop->degree == 2)
		return;

	tiifhe_poly_add(idx, &rop->poly[2], &op1->poly[2], &op2->poly[2]);
}

void
tiifhe_bgv_idx_addc(size_t idx, tiifhe_Ciphertext *rop, tiifhe_Ciphertext *ct, const tiifhe_Encoding *p)
{

	ciphertext_ntt(idx, ct);

	rop->degree = ct->degree;
	rop->drop = ct->drop;
	rop->ntt = ct->ntt;

	if (rop->drop != 0)
		return;

	tiifhe_poly_add(idx, &rop->poly[0], &ct->poly[0], &p->value);
}

void
tiifhe_bgv_idx_decrypt(size_t idx, tiifhe_Encoding *p, const tiifhe_KeySecret *sk, tiifhe_Ciphertext *ct)
{
	tiifhe_Poly *s;

	ciphertext_ntt(idx, ct);
	p->drop = ct->drop;

	if (p->drop != 0)
		return;

	s = tiifhe_util_alloc(1, sizeof *s);

	#if TIIFHE_EXPAND_SEED
	memcpy(s, &sk->poly[idx], sizeof *s);
	#else
	tiifhe_poly_init(idx, s, &sk->poly);
	tiifhe_poly_ntt(idx, s);
	#endif /* TIIFHE_EXPAND_SEED */
	tiifhe_poly_muladd(idx, &p->value, &ct->poly[1], s, &ct->poly[0]);

	if (ct->degree == 3) {
		tiifhe_poly_mul(idx, s, s, s);
		tiifhe_poly_muladd(idx, &p->value, &ct->poly[2], s, &p->value);
	}

	tiifhe_util_dealloc(s);
}

void
tiifhe_bgv_idx_encode(size_t idx, tiifhe_Encoding *p, const tiifhe_Message *m)
{
	
	tiifhe_poly_mod(idx, &p->value, m->value);
	tiifhe_poly_ntt(idx, &p->value);
	p->drop = 0;
}

void
tiifhe_bgv_idx_encrypt(size_t idx, tiifhe_Ciphertext *ct, tiifhe_KeyPublic *pk, const tiifhe_Encoding *p, tiifhe_Seed *seed)
{
	tiifhe_Poly *coins;

	ct->degree = 2;
	ct->drop = p->drop;
	ct->ntt = 1;

	if (ct->drop != 0)
		return;

	/* sample and init coins: a, u, e0, e1 */
	coins = tiifhe_util_alloc(4, sizeof *coins);
	#if TIIFHE_EXPAND_SEED
	memcpy(&coins[0], &pk->rand, sizeof coins[0]);
	#else
	tiifhe_poly_sample_uniform(idx, &coins[0], &pk->seed);
	pk->seed.init = 1;
	#endif /* TIIFHE_EXPAND_SEED */

	tiifhe_poly_sample_secret(&coins[1], seed);
	tiifhe_poly_init(idx, &coins[1], &coins[1]);
	tiifhe_poly_ntt(idx, &coins[1]);

	tiifhe_poly_sample_error(&coins[2], seed);
	tiifhe_poly_init_error(idx, &coins[2], &coins[2]);

	tiifhe_poly_sample_error(&coins[3], seed);
	tiifhe_poly_init_error(idx, &coins[3], &coins[3]);
	seed->init = 1;

	/* compute ciphertext */
	tiifhe_poly_muladd(idx, &ct->poly[0], &pk->poly, &coins[1], &coins[2]);
	tiifhe_poly_add(idx, &ct->poly[0], &ct->poly[0], &p->value);
	tiifhe_poly_muladd(idx, &ct->poly[1], &coins[0], &coins[1], &coins[3]);
	memset(&ct->poly[2], 0, sizeof ct->poly[2]);

	tiifhe_util_dealloc(coins);
}

void
tiifhe_bgv_idx_keygen_public(size_t idx, tiifhe_KeyPublic *pk, const tiifhe_KeySecret *sk, tiifhe_Seed *seed)
{
	tiifhe_Poly *p;
	tiifhe_Seed sd;

	p = tiifhe_util_alloc(3, sizeof *p);

	/* sample and init polynomials */
	tiifhe_random_seed(&sd);
	tiifhe_poly_sample_uniform(idx, &p[0], &sd);
	sd.init = 1;

	#if TIIFHE_EXPAND_SEED
	memcpy(&pk->rand, &p[0], sizeof pk->rand);
	memcpy(&p[1], &sk->poly[idx], sizeof p[1]);
	#else
	memcpy(&pk->seed, &sd, sizeof pk->seed);
	tiifhe_poly_init(idx, &p[1], &sk->poly);
	tiifhe_poly_ntt(idx, &p[1]);
	#endif /* TIIFHE_EXPAND_SEED */

	tiifhe_poly_sample_error(&p[2], seed);
	tiifhe_poly_init_error(idx, &p[2], &p[2]);
	seed->init = 1;

	/* compute public key */
	tiifhe_poly_neg(idx, &p[0], &p[0]);
	tiifhe_poly_muladd(idx, &pk->poly, &p[0], &p[1], &p[2]);

	tiifhe_util_dealloc(p);
}

void
tiifhe_bgv_idx_keygen_switch(size_t idx, tiifhe_KeySwitch *ksw, const tiifhe_KeySecret *sk, const tiifhe_KeySecret *sk2, tiifhe_Seed *seed)
{
	tiifhe_Poly *a, *e, *s;
	tiifhe_Seed sd;
	size_t k;

	a = tiifhe_util_alloc(TIIFHE_OMEGA, sizeof *a);
	e = tiifhe_util_alloc(TIIFHE_OMEGA, sizeof *e);
	s = tiifhe_util_alloc(2, sizeof *s);

	tiifhe_random_seed(&sd);
	for (k = 0; k < TIIFHE_OMEGA; ++k) {
		tiifhe_poly_sample_uniform(idx, &a[k], &sd);
	}
	sd.init = 1;

	#if TIIFHE_EXPAND_SEED
	memcpy(&ksw->rand, a, TIIFHE_OMEGA * sizeof(tiifhe_Poly));
	#else
	memcpy(&ksw->seed, &sd, sizeof ksw->seed);
	#endif /* TIIFHE_EXPAND_SEED */

	for (k = 0; k < TIIFHE_OMEGA; ++k) {
		tiifhe_poly_sample_error(&e[k], &seed[k]);
		tiifhe_poly_init_error(idx, &e[k], &e[k]);
		seed[k].init = 1;
	}

	#if TIIFHE_EXPAND_SEED
	memcpy(&s[0], &sk->poly[idx], sizeof s[0]);
	#else
	tiifhe_poly_init(idx, &s[0], &sk->poly);
	tiifhe_poly_ntt(idx, &s[0]);
	#endif /* TIIFHE_EXPAND_SEED */

	#if TIIFHE_EXPAND_SEED
	memcpy(&s[1], &sk2->poly[idx], sizeof s[1]);
	#else
	tiifhe_poly_init(idx, &s[1], &sk2->poly);
	tiifhe_poly_ntt(idx, &s[1]);
	#endif /* TIIFHE_EXPAND_SEED */

	for (k = 0; k < TIIFHE_OMEGA; ++k) {
		tiifhe_poly_neg(idx, &a[k], &a[k]);
		tiifhe_poly_muladd(idx, &ksw->poly[k], &a[k], &s[0], &e[k]);

		if (idx % TIIFHE_OMEGA == k)
			tiifhe_poly_mulcadd(idx, &ksw->poly[k], &s[1], tiifhe_q[idx].P, &ksw->poly[k]);
	}

	tiifhe_util_dealloc(a);
	tiifhe_util_dealloc(e);
	tiifhe_util_dealloc(s);
}

void
tiifhe_bgv_idx_keygen_switch_mul(size_t idx, tiifhe_KeySwitch *ksw, const tiifhe_KeySecret *sk, tiifhe_Seed *seed)
{
	tiifhe_KeySecret *s;
	size_t i;
	(void)i;

	#if TIIFHE_EXPAND_SEED
	s = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *s);
	for (i = 0; i < TIIFHE_QPLEN; ++i) {
		tiifhe_poly_mul(idx, &s->poly[idx], &sk->poly[idx], &sk->poly[idx]);
	}
	#else
	s = tiifhe_util_alloc(1, sizeof *s);
	tiifhe_poly_init(idx, &s->poly, &sk->poly);
	tiifhe_poly_ntt(idx, &s->poly);
	tiifhe_poly_mul(idx, &s->poly, &s->poly, &s->poly);
	tiifhe_poly_intt(idx, &s->poly);
	tiifhe_poly_deinit(idx, &s->poly, &s->poly);
	tiifhe_poly_cmod(idx, &s->poly, &s->poly);
	#endif /* TIIFHE_EXPAND_SEED */

	tiifhe_bgv_idx_keygen_switch(idx, ksw, sk, s, seed);

	tiifhe_util_dealloc(s);
}

void
tiifhe_bgv_idx_keygen_switch_rot(size_t idx, tiifhe_KeySwitch *ksw, const tiifhe_KeySecret *sk, size_t steps, tiifhe_Seed *seed)
{
	tiifhe_KeySecret *s;
	size_t i;
	(void)i;

	#if TIIFHE_EXPAND_SEED
	s = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *s);
	for (i = 0; i < TIIFHE_QPLEN; ++i) {
		tiifhe_poly_rot_ntt(idx, &s->poly[idx], &sk->poly[idx], steps);
	}
	#else
	s = tiifhe_util_alloc(1, sizeof *s);
	tiifhe_poly_rot_intt(idx, &s->poly, &sk->poly, steps);
	#endif /* TIIFHE_EXPAND_SEED */

	tiifhe_bgv_idx_keygen_switch(idx, ksw, sk, s, seed);

	tiifhe_util_dealloc(s);
}

void
tiifhe_bgv_idx_keyswitch_delta(size_t idx, tiifhe_Delta *delta, tiifhe_Ciphertext *ct)
{

	tiifhe_bgv_idx_modswitch_delta(idx, delta, ct);
}

void
tiifhe_bgv_idx_keyswitch_dot(size_t idx, tiifhe_Ciphertext *ct, tiifhe_CiphertextSwitch *csw, tiifhe_KeySwitch *ksw)
{
	tiifhe_Poly *a;
	size_t k;

	ciphertext_switch_ntt(idx, csw);
	memset(ct, 0, sizeof *ct);

	ct->degree = 2;
	ct->drop = csw->ct.drop;
	ct->ntt = 1;

	if (ct->drop != 0)
		return;

	a = tiifhe_util_alloc(1, sizeof *a);

	for (k = 0; k < TIIFHE_OMEGA; ++k)
		tiifhe_poly_muladd(idx, &ct->poly[0], &csw->ext[k], &ksw->poly[k], &ct->poly[0]);
	tiifhe_poly_mulcadd(idx, &ct->poly[0], &csw->ct.poly[0], tiifhe_q[idx].P, &ct->poly[0]);

	for (k = 0; k < TIIFHE_OMEGA; ++k) {
		#if TIIFHE_EXPAND_SEED
		tiifhe_poly_muladd(idx, &ct->poly[1], &csw->ext[k], &ksw->rand[k], &ct->poly[1]);
		#else
		tiifhe_poly_sample_uniform(idx, a, &ksw->seed);
		tiifhe_poly_muladd(idx, &ct->poly[1], &csw->ext[k], a, &ct->poly[1]);
		#endif /* TIIFHE_EXPAND_SEED */
	}
	tiifhe_poly_mulcadd(idx, &ct->poly[1], &csw->ct.poly[1], csw->poly == 2 ? tiifhe_q[idx].P : 0, &ct->poly[1]);
	#if !TIIFHE_EXPAND_SEED
	ksw->seed.init = 1;
	#endif /* !TIIFHE_EXPAND_SEED */

	tiifhe_util_dealloc(a);
}

void
tiifhe_bgv_idx_keyswitch_switch(size_t idx, tiifhe_Ciphertext *ct, tiifhe_Delta *delta)
{
	tiifhe_Delta *tmp;
	size_t P[TIIFHE_PLEN], i;

	if (ct->drop != 0)
		return;

	tmp = tiifhe_util_alloc(1, sizeof *tmp);

	for (i = 0; i < TIIFHE_PLEN; ++i)
		P[i] = TIIFHE_QLEN + i;

	delta_ext(idx, tmp, delta, P, TIIFHE_PLEN);
	tiifhe_bgv_idx_modswitch(idx, ct, tmp);

	tiifhe_util_dealloc(tmp);
}

void
tiifhe_bgv_idx_modswitch(size_t idx, tiifhe_Ciphertext *ct, tiifhe_Delta *delta)
{

	ciphertext_ntt(idx, ct);
	ciphertext_ntt(idx, &delta->ct);

	if (ct->drop != 0)
		return;

	assert(ct->degree == delta->ct.degree);

	tiifhe_poly_mulcadd(idx, &ct->poly[0], &delta->ct.poly[0], tiifhe_q[idx].t, &ct->poly[0]);
	tiifhe_poly_mulc(idx, &ct->poly[0], &ct->poly[0], delta->inv);

	tiifhe_poly_mulcadd(idx, &ct->poly[1], &delta->ct.poly[1], tiifhe_q[idx].t, &ct->poly[1]);
	tiifhe_poly_mulc(idx, &ct->poly[1], &ct->poly[1], delta->inv);

	if (ct->degree == 2)
		return;

	tiifhe_poly_mulcadd(idx, &ct->poly[2], &delta->ct.poly[2], tiifhe_q[idx].t, &ct->poly[2]);
	tiifhe_poly_mulc(idx, &ct->poly[2], &ct->poly[2], delta->inv);
}


void
tiifhe_bgv_idx_modswitch_delta(size_t idx, tiifhe_Delta *delta, tiifhe_Ciphertext *ct)
{
	tiifhe_Digit inv;

	ciphertext_ntt(idx, ct);

	delta->ct.degree = ct->degree;
	delta->ct.drop = ct->drop;
	delta->ct.ntt = 0;
	delta->inv = 0;
	delta->idx = idx;

	if (delta->ct.drop != 0)
		return;

	inv = tiifhe_mod_neg(idx, tiifhe_mod_inv(idx, tiifhe_q[idx].t));
	ct->drop = 1;

	tiifhe_poly_mulc(idx, &delta->ct.poly[0], &ct->poly[0], inv);
	tiifhe_poly_intt(idx, &delta->ct.poly[0]);

	tiifhe_poly_mulc(idx, &delta->ct.poly[1], &ct->poly[1], inv);
	tiifhe_poly_intt(idx, &delta->ct.poly[1]);

	if (delta->ct.degree == 2)
		return;

	tiifhe_poly_mulc(idx, &delta->ct.poly[2], &ct->poly[2], inv);
	tiifhe_poly_intt(idx, &delta->ct.poly[2]);
}

void
tiifhe_bgv_idx_modswitch_ext(size_t idx, tiifhe_Ciphertext *ct, tiifhe_Delta *delta, size_t len)
{
	tiifhe_Delta *cpy;
	size_t base[TIIFHE_QPLEN], i;

	if (ct->drop != 0)
		return;

	cpy = tiifhe_util_alloc(1, sizeof *cpy);
	for (i = 0; i < len; ++i)
		base[i] = delta[i].idx;

	delta_ext(idx, cpy, delta, base, len);
	tiifhe_bgv_idx_modswitch(idx, ct, cpy);

	tiifhe_util_dealloc(cpy);
}

void
tiifhe_bgv_idx_mul(size_t idx, tiifhe_Ciphertext *rop, tiifhe_Ciphertext *op1, tiifhe_Ciphertext *op2)
{

	if (op1->degree != 2 || op2->degree != 2) {
		fprintf(stdout, "fhelib: trying to multiply ciphertexts with degree greater two at idx %lu (%lu, %lu)\n", idx, op1->degree, op2->degree);
		exit(1);
	}

	ciphertext_ntt(idx, op1);
	ciphertext_ntt(idx, op2);

	rop->degree = 3;
	rop->drop = op1->drop + op2->drop;
	rop->ntt = 1;

	if (rop->drop != 0)
		return;

	tiifhe_poly_mul(idx, &rop->poly[2], &op1->poly[1], &op2->poly[1]);
	tiifhe_poly_addmul(idx, &rop->poly[1], &op1->poly[0], &op1->poly[1], &op2->poly[0], &op2->poly[1]);
	tiifhe_poly_mul(idx, &rop->poly[0], &op1->poly[0], &op2->poly[0]);
	tiifhe_poly_sub(idx, &rop->poly[1], &rop->poly[1], &rop->poly[0]);
	tiifhe_poly_sub(idx, &rop->poly[1], &rop->poly[1], &rop->poly[2]);
}

void
tiifhe_bgv_idx_mulc(size_t idx, tiifhe_Ciphertext *rop, tiifhe_Ciphertext *ct, const tiifhe_Encoding *p)
{

	ciphertext_ntt(idx, ct);

	rop->degree = ct->degree;
	rop->drop = ct->drop;
	rop->ntt = ct->ntt;

	if (rop->drop != 0)
		return;

	tiifhe_poly_mul(idx, &rop->poly[0], &ct->poly[0], &p->value);
	tiifhe_poly_mul(idx, &rop->poly[1], &ct->poly[1], &p->value);

	if (rop->degree == 2)
		return;

	tiifhe_poly_mul(idx, &rop->poly[2], &ct->poly[2], &p->value);
}

void
tiifhe_bgv_idx_rot(size_t idx, tiifhe_Ciphertext *rop, tiifhe_Ciphertext *ct, size_t steps)
{

	if (rop == ct) {
		fprintf(stdout, "fhelib: trying to rotate ciphertext in-place at idx %lu\n", idx);
		exit(1);
	}

	rop->degree = ct->degree;
	rop->drop = ct->drop;
	rop->ntt = ct->ntt;

	if (rop->drop != 0)
		return;

	if (rop->ntt) {
		tiifhe_poly_rot_ntt(idx, &rop->poly[0], &ct->poly[0], steps);
		tiifhe_poly_rot_ntt(idx, &rop->poly[1], &ct->poly[1], steps);

		if (rop->degree == 2)
			return;

		tiifhe_poly_rot_ntt(idx, &rop->poly[2], &ct->poly[2], steps);
	} else {
		tiifhe_poly_rot_intt(idx, &rop->poly[0], &ct->poly[0], steps);
		tiifhe_poly_rot_intt(idx, &rop->poly[1], &ct->poly[1], steps);

		if (rop->degree == 2)
			return;

		tiifhe_poly_rot_intt(idx, &rop->poly[2], &ct->poly[2], steps);
	}
}

void
tiifhe_bgv_idx_rot_csw(size_t idx, tiifhe_CiphertextSwitch *rop, tiifhe_CiphertextSwitch *csw, size_t steps)
{
	size_t k;

	if (rop == csw) {
		fprintf(stdout, "fhelib: trying to rotate extended ciphertext in-place at idx %lu\n", idx);
		exit(1);
	}

	tiifhe_bgv_idx_rot(idx, &rop->ct, &csw->ct, steps);
	if (rop->ct.drop != 0)
		return;

	for (k = 0; k < TIIFHE_OMEGA; ++k) {
		if (rop->ct.ntt)
			tiifhe_poly_rot_ntt(idx, &rop->ext[k], &csw->ext[k], steps);
		else
			tiifhe_poly_rot_intt(idx, &rop->ext[k], &csw->ext[k], steps);
	}
}

void
tiifhe_bgv_idx_rot_csw_inplace(size_t idx, tiifhe_CiphertextSwitch *csw, size_t steps)
{
	tiifhe_CiphertextSwitch *cpy;

	cpy = tiifhe_util_alloc(1, sizeof *cpy);
	memcpy(cpy, csw, sizeof *cpy);

	tiifhe_bgv_idx_rot_csw(idx, csw, cpy, steps);

	tiifhe_util_dealloc(cpy);
}

void
tiifhe_bgv_idx_rot_inplace(size_t idx, tiifhe_Ciphertext *ct, size_t steps)
{
	tiifhe_Ciphertext *cpy;

	cpy = tiifhe_util_alloc(1, sizeof *cpy);
	memcpy(cpy, ct, sizeof *cpy);

	tiifhe_bgv_idx_rot(idx, ct, cpy, steps);

	tiifhe_util_dealloc(cpy);
}
