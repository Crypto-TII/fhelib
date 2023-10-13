#include "tiifhe.h"

#define LEN(a) (sizeof a / sizeof *a)
gmp_randstate_t state;

void
msg_add(tiifhe_Message *rop, const tiifhe_Message *op1, const tiifhe_Message *op2)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j) {
		mpz_add(rop->value[j], op1->value[j], op2->value[j]);
		mpz_mod(rop->value[j], rop->value[j], tiifhe_t.value);
	}
}

tiifhe_Message *
msg_alloc_rand(void)
{
	tiifhe_Message *m;
	size_t j;

	m = tiifhe_msg_alloc(1);
	for (j = 0; j < TIIFHE_N; ++j) {
		mpz_urandomm(m->value[j], state, tiifhe_t.value);
	}

	return m;
}

int
msg_equal(const tiifhe_Message *m1, const tiifhe_Message *m2)
{
	size_t j;
	int cmp;

	cmp = 0;
	for (j = 0; j < TIIFHE_N; ++j) {
		cmp += mpz_cmp(m1->value[j], m2->value[j]);
	}

	return cmp == 0;
}

void
msg_mul(tiifhe_Message *rop, const tiifhe_Message *op1, const tiifhe_Message *op2)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j) {
		mpz_mul(rop->value[j], op1->value[j], op2->value[j]);
		mpz_mod(rop->value[j], rop->value[j], tiifhe_t.value);
	}
}

void
msg_rot(tiifhe_Message *rop, const tiifhe_Message *m, size_t steps)
{
	mpz_t tmp;
	size_t i, j;

	mpz_init(tmp);

	for (i = 0; i < steps; ++i) {
		mpz_set(tmp, m->value[0]);
		for (j = 0; j < TIIFHE_N / 2 - 1; ++j) {
			mpz_set(rop->value[j], m->value[j + 1]);
		}
		mpz_set(rop->value[j], tmp);

		mpz_set(tmp, m->value[TIIFHE_N / 2]);
		for (j = TIIFHE_N / 2; j < TIIFHE_N - 1; ++j) {
			mpz_set(rop->value[j], m->value[j + 1]);
		}
		mpz_set(rop->value[j], tmp);
	}

	mpz_clear(tmp);
}

void
test_pack(void)
{
	tiifhe_Message *cmp, *m;
	tiifhe_Encoding *p;
	size_t i;

	cmp = msg_alloc_rand();
	m = tiifhe_msg_alloc(1);
	p = tiifhe_util_alloc(TIIFHE_QLEN, sizeof *p);

	tiifhe_msg_pack(m, cmp);
	for (i = 0; i < TIIFHE_QLEN; ++i) {
		tiifhe_bgv_idx_encode(i, &p[i], m);
	}
	tiifhe_bgv_decode(m, p);
	tiifhe_msg_unpack(m, m);

	assert(msg_equal(cmp, m));

	tiifhe_util_dealloc(p);
	tiifhe_msg_dealloc(cmp, 1);
	tiifhe_msg_dealloc(m, 1);
}

void
test_encrypt(void)
{
	tiifhe_Seed seed[2];
	tiifhe_KeySecret *sk;
	tiifhe_KeyPublic *pk;

	tiifhe_Message *cmp, *m;
	tiifhe_Ciphertext *ct;
	tiifhe_Encoding *p;
	size_t i;

	sk = tiifhe_util_alloc(1, sizeof *sk);
	pk = tiifhe_util_alloc(1, sizeof *pk);
	cmp = msg_alloc_rand();
	m = tiifhe_msg_alloc(1);
	ct = tiifhe_util_alloc(1, sizeof *ct);
	p = tiifhe_util_alloc(TIIFHE_QLEN, sizeof *p);

	tiifhe_bgv_keygen_secret(sk);
	for (i = 0; i < LEN(seed); ++i) {
		tiifhe_random_seed(&seed[i]);
	}

	tiifhe_msg_pack(m, cmp);
	for (i = 0; i < TIIFHE_QLEN; ++i) {
		tiifhe_bgv_idx_encode(i, &p[i], m);
		tiifhe_bgv_idx_keygen_public(i, pk, sk, &seed[0]);
		tiifhe_bgv_idx_encrypt(i, ct, pk, &p[i], &seed[1]);
		tiifhe_bgv_idx_decrypt(i, &p[i], sk, ct);
	}
	tiifhe_bgv_decode(m, p);
	tiifhe_msg_unpack(m, m);
	assert(msg_equal(cmp, m));

	tiifhe_util_dealloc(sk);
	tiifhe_util_dealloc(pk);
	tiifhe_msg_dealloc(cmp, 1);
	tiifhe_msg_dealloc(m, 1);
	tiifhe_util_dealloc(ct);
	tiifhe_util_dealloc(p);
}

void
test_arithmetic(void)
{
	tiifhe_Seed seed[4];
	tiifhe_KeySecret *sk;
	tiifhe_KeyPublic *pk;

	tiifhe_Message *cmp, *m[4];
	tiifhe_Ciphertext *ct;
	tiifhe_Encoding *p;
	size_t i;

	sk = tiifhe_util_alloc(1, sizeof *sk);
	pk = tiifhe_util_alloc(1, sizeof *pk);
	cmp = tiifhe_msg_alloc(1);
	m[0] = msg_alloc_rand();
	m[1] = msg_alloc_rand();
	m[2] = msg_alloc_rand();
	m[3] = msg_alloc_rand();
	ct = tiifhe_util_alloc(3, sizeof *ct);
	p = tiifhe_util_alloc(TIIFHE_QLEN, sizeof *p);

	tiifhe_bgv_keygen_secret(sk);
	for (i = 0; i < LEN(seed); ++i) {
		tiifhe_random_seed(&seed[i]);
	}

	msg_mul(cmp, m[0], m[3]);
	msg_add(cmp, cmp, m[1]);
	msg_mul(cmp, cmp, m[2]);

	tiifhe_msg_pack(m[0], m[0]);
	tiifhe_msg_pack(m[1], m[1]);
	tiifhe_msg_pack(m[2], m[2]);
	tiifhe_msg_pack(m[3], m[3]);

	for (i = 0; i < TIIFHE_QLEN; ++i) {
		tiifhe_bgv_idx_keygen_public(i, pk, sk, &seed[0]);

		tiifhe_bgv_idx_encode(i, &p[i], m[0]);
		tiifhe_bgv_idx_encrypt(i, &ct[0], pk, &p[i], &seed[1]);
		tiifhe_bgv_idx_encode(i, &p[i], m[1]);
		tiifhe_bgv_idx_encrypt(i, &ct[1], pk, &p[i], &seed[2]);
		tiifhe_bgv_idx_encode(i, &p[i], m[2]);
		tiifhe_bgv_idx_encrypt(i, &ct[2], pk, &p[i], &seed[3]);
		tiifhe_bgv_idx_encode(i, &p[i], m[3]);

		tiifhe_bgv_idx_mulc(i, &ct[0], &ct[0], &p[i]);
		tiifhe_bgv_idx_add(i, &ct[0], &ct[0], &ct[1]);
		tiifhe_bgv_idx_mul(i, &ct[0], &ct[0], &ct[2]);

		tiifhe_bgv_idx_decrypt(i, &p[i], sk, &ct[0]);
	}

	tiifhe_bgv_decode(m[0], p);
	tiifhe_msg_unpack(m[0], m[0]);
	assert(msg_equal(cmp, m[0]));

	tiifhe_util_dealloc(sk);
	tiifhe_util_dealloc(pk);
	tiifhe_msg_dealloc(cmp, 1);
	for (i = 0; i < LEN(m); ++i) {
		tiifhe_msg_dealloc(m[i], 1);
	}
	tiifhe_util_dealloc(ct);
	tiifhe_util_dealloc(p);
}

void
test_modswitch(void)
{
	tiifhe_Seed seed[2];
	tiifhe_KeySecret *sk;
	tiifhe_KeyPublic *pk;

	tiifhe_Message *cmp, *m;
	tiifhe_Ciphertext *ct;
	tiifhe_Delta *delta;
	tiifhe_Encoding *p;
	size_t i;

	sk = tiifhe_util_alloc(1, sizeof *sk);
	pk = tiifhe_util_alloc(1, sizeof *pk);
	cmp = msg_alloc_rand();
	m = tiifhe_msg_alloc(1);
	ct = tiifhe_util_alloc(TIIFHE_QLEN, sizeof *ct);
	delta = tiifhe_util_alloc(2, sizeof *delta);
	p = tiifhe_util_alloc(TIIFHE_QLEN, sizeof *p);

	tiifhe_bgv_keygen_secret(sk);
	for (i = 0; i < LEN(seed); ++i) {
		tiifhe_random_seed(&seed[i]);
	}

	tiifhe_msg_pack(m, cmp);
	for (i = 0; i < TIIFHE_QLEN; ++i) {
		tiifhe_bgv_idx_keygen_public(i, pk, sk, &seed[0]);
		tiifhe_bgv_idx_encode(i, &p[i], m);
		tiifhe_bgv_idx_encrypt(i, &ct[i], pk, &p[i], &seed[1]);

		if (i == 0) {
			tiifhe_bgv_idx_modswitch_delta(i, &delta[i], &ct[i]);
		}
		tiifhe_bgv_idx_modswitch_ext(i, &ct[i], &delta[0], 1);

		if (i == 1) {
			tiifhe_bgv_idx_modswitch_delta(i, &delta[i], &ct[i]);
		}
		tiifhe_bgv_idx_modswitch_ext(i, &ct[i], &delta[1], 1);

		tiifhe_bgv_idx_decrypt(i, &p[i], sk, &ct[i]);
	}
	tiifhe_bgv_decode(m, p);
	tiifhe_msg_unpack(m, m);
	assert(msg_equal(cmp, m));

	tiifhe_msg_pack(m, cmp);
	msg_mul(cmp, cmp, cmp);
	msg_add(cmp, cmp, cmp);
	for (i = 0; i < TIIFHE_QLEN; ++i) {
		tiifhe_bgv_idx_keygen_public(i, pk, sk, &seed[0]);
		tiifhe_bgv_idx_encode(i, &p[i], m);
		tiifhe_bgv_idx_encrypt(i, &ct[i], pk, &p[i], &seed[1]);

		if (i == 0) {
			tiifhe_bgv_idx_modswitch_delta(i, &delta[i], &ct[i]);
		}
		tiifhe_bgv_idx_modswitch_ext(i, &ct[i], &delta[0], 1);
		tiifhe_bgv_idx_mul(i, &ct[i], &ct[i], &ct[i]);

		if (i == 1) {
			tiifhe_bgv_idx_modswitch_delta(i, &delta[i], &ct[i]);
		}
		tiifhe_bgv_idx_modswitch_ext(i, &ct[i], &delta[1], 1);
		tiifhe_bgv_idx_add(i, &ct[i], &ct[i], &ct[i]);

		tiifhe_bgv_idx_decrypt(i, &p[i], sk, &ct[i]);
	}
	tiifhe_bgv_decode(m, p);
	tiifhe_msg_unpack(m, m);
	assert(msg_equal(cmp, m));

	tiifhe_msg_pack(m, cmp);
	for (i = 0; i < TIIFHE_QLEN; ++i) {
		tiifhe_bgv_idx_keygen_public(i, pk, sk, &seed[0]);
		tiifhe_bgv_idx_encode(i, &p[i], m);
		tiifhe_bgv_idx_encrypt(i, &ct[i], pk, &p[i], &seed[1]);
	}
	tiifhe_bgv_idx_modswitch_delta(0, &delta[0], &ct[0]);
	tiifhe_bgv_idx_modswitch_delta(1, &delta[1], &ct[1]);
	for (i = 0; i < TIIFHE_QLEN; ++i) {
		tiifhe_bgv_idx_modswitch_ext(i, &ct[i], delta, 2);
		tiifhe_bgv_idx_decrypt(i, &p[i], sk, &ct[i]);
	}
	tiifhe_bgv_decode(m, p);
	tiifhe_msg_unpack(m, m);
	assert(msg_equal(cmp, m));

	tiifhe_util_dealloc(delta);
	tiifhe_util_dealloc(ct);
	tiifhe_util_dealloc(p);
	tiifhe_msg_dealloc(cmp, 1);
	tiifhe_msg_dealloc(m, 1);
	tiifhe_util_dealloc(pk);
	tiifhe_util_dealloc(sk);
}

void
test_keyswitch(void)
{
	tiifhe_Seed seed[TIIFHE_OMEGA + 3];
	tiifhe_KeySecret *sk;
	tiifhe_KeyPublic *pk;
	tiifhe_KeySwitch *ksw, *ksw2, *kswr;

	tiifhe_Message *cmp, *m;
	tiifhe_Ciphertext *ct;
	tiifhe_CiphertextSwitch *csw, *cswr;
	tiifhe_Delta *delta;
	tiifhe_Encoding *p;
	size_t i;

	sk = tiifhe_util_alloc(1, sizeof *sk);
	pk = tiifhe_util_alloc(1, sizeof *pk);
	ksw = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *ksw);
	ksw2 = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *ksw2);
	kswr = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *kswr);
	cmp = msg_alloc_rand();
	m = tiifhe_msg_alloc(1);
	ct = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *ct);
	csw = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *csw);
	cswr = tiifhe_util_alloc(1, sizeof *cswr);
	delta = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *delta);
	p = tiifhe_util_alloc(TIIFHE_QLEN, sizeof *p);

	tiifhe_bgv_keygen_secret(sk);

	for (i = 0; i < LEN(seed); ++i) {
		tiifhe_random_seed(&seed[i]);
	}
	for (i = 0; i < TIIFHE_QPLEN; ++i) {
		tiifhe_bgv_idx_keygen_switch(i, &ksw[i], sk, sk, &seed[3]);
	}

	for (i = 3; i < LEN(seed); ++i) {
		tiifhe_random_seed(&seed[i]);
	}
	for (i = 0; i < TIIFHE_QPLEN; ++i) {
		tiifhe_bgv_idx_keygen_switch_rot(i, &kswr[i], sk, 1, &seed[3]);
	}

	for (i = 3; i < LEN(seed); ++i) {
		tiifhe_random_seed(&seed[i]);
	}
	for (i = 0; i < TIIFHE_QPLEN; ++i) {
		tiifhe_bgv_idx_keygen_switch_mul(i, &ksw2[i], sk, &seed[3]);
	}

	tiifhe_msg_pack(m, cmp);
	for (i = 0; i < TIIFHE_QLEN; ++i) {
		tiifhe_bgv_idx_keygen_public(i, pk, sk, &seed[0]);
		tiifhe_bgv_idx_encode(i, &p[i], m);
		tiifhe_bgv_idx_encrypt(i, &ct[i], pk, &p[i], &seed[1]);
	}
	tiifhe_bgv_keyswitch(ct, ksw, 1);
	for (i = 0; i < TIIFHE_QLEN; ++i) {
		tiifhe_bgv_idx_decrypt(i, &p[i], sk, &ct[i]);
	}
	tiifhe_bgv_decode(m, p);
	tiifhe_msg_unpack(m, m);
	assert(msg_equal(cmp, m));

	tiifhe_msg_pack(m, cmp);
	msg_rot(cmp, cmp, 1);
	for (i = 0; i < TIIFHE_QLEN; ++i) {
		tiifhe_bgv_idx_keygen_public(i, pk, sk, &seed[0]);
		tiifhe_bgv_idx_encode(i, &p[i], m);
		tiifhe_bgv_idx_encrypt(i, &ct[i], pk, &p[i], &seed[1]);
		tiifhe_bgv_idx_rot_inplace(i, &ct[i], 1);
	}
	tiifhe_bgv_keyswitch(ct, kswr, 1);
	for (i = 0; i < TIIFHE_QLEN; ++i) {
		tiifhe_bgv_idx_decrypt(i, &p[i], sk, &ct[i]);
	}
	tiifhe_bgv_decode(m, p);
	tiifhe_msg_unpack(m, m);
	assert(msg_equal(cmp, m));

	tiifhe_msg_pack(m, cmp);
	for (i = 0; i < TIIFHE_QLEN; ++i) {
		tiifhe_bgv_idx_keygen_public(i, pk, sk, &seed[0]);
		tiifhe_bgv_idx_encode(i, &p[i], m);
		tiifhe_bgv_idx_encrypt(i, &ct[i], pk, &p[i], &seed[1]);
		tiifhe_bgv_idx_mul(i, &ct[i], &ct[i], &ct[i]);
	}
	tiifhe_bgv_keyswitch(ct, ksw2, 2);
	for (i = 0; i < TIIFHE_QLEN; ++i)
		tiifhe_bgv_idx_decrypt(i, &p[i], sk, &ct[i]);
	tiifhe_bgv_decode(m, p);
	tiifhe_msg_unpack(m, m);
	msg_mul(cmp, cmp, cmp);
	assert(msg_equal(cmp, m));

	tiifhe_msg_pack(m, cmp);
	for (i = 0; i < TIIFHE_QLEN; ++i) {
		tiifhe_bgv_idx_keygen_public(i, pk, sk, &seed[0]);
		tiifhe_bgv_idx_encode(i, &p[i], m);
		tiifhe_bgv_idx_encrypt(i, &ct[i], pk, &p[i], &seed[1]);

		if (i == 0) {
			tiifhe_bgv_idx_modswitch_delta(0, &delta[0], &ct[0]);
		}
		tiifhe_bgv_idx_modswitch_ext(i, &ct[i], &delta[0], 1);

		if (i == 1) {
			tiifhe_bgv_idx_modswitch_delta(1, &delta[1], &ct[1]);
		}
		tiifhe_bgv_idx_modswitch_ext(i, &ct[i], &delta[1], 1);
	}
	tiifhe_bgv_keyswitch(ct, ksw, 1);
	for (i = 0; i < TIIFHE_QLEN; ++i) {
		tiifhe_bgv_idx_decrypt(i, &p[i], sk, &ct[i]);
	}
	tiifhe_bgv_decode(m, p);
	tiifhe_msg_unpack(m, m);
	assert(msg_equal(cmp, m));

	tiifhe_msg_pack(m, cmp);
	for (i = 0; i < TIIFHE_QLEN; ++i) {
		tiifhe_bgv_idx_keygen_public(i, pk, sk, &seed[0]);
		tiifhe_bgv_idx_encode(i, &p[i], m);
		tiifhe_bgv_idx_encrypt(i, &ct[i], pk, &p[i], &seed[1]);
	}
	tiifhe_bgv_keyswitch_ext(csw, ct, 1);
	for (i = 0; i < TIIFHE_QPLEN; ++i) {
		tiifhe_bgv_idx_rot_csw(i, cswr, &csw[i], 1);
		tiifhe_bgv_idx_keyswitch_dot(i, &ct[i], cswr, &kswr[i]);
	}
	for (i = TIIFHE_QLEN; i < TIIFHE_QPLEN; ++i) {
		tiifhe_bgv_idx_keyswitch_delta(i, &delta[i - TIIFHE_QLEN], &ct[i]);
	}
	for (i = 0; i < TIIFHE_QLEN; ++i) {
		tiifhe_bgv_idx_keyswitch_switch(i, &ct[i], delta);
		tiifhe_bgv_idx_decrypt(i, &p[i], sk, &ct[i]);
	}
	tiifhe_bgv_decode(m, p);
	tiifhe_msg_unpack(m, m);
	msg_rot(cmp, cmp, 1);
	assert(msg_equal(cmp, m));

	tiifhe_util_dealloc(sk);
	tiifhe_util_dealloc(pk);
	tiifhe_util_dealloc(ksw);
	tiifhe_util_dealloc(ksw2);
	tiifhe_util_dealloc(kswr);

	tiifhe_msg_dealloc(cmp, 1);
	tiifhe_msg_dealloc(m, 1);
	tiifhe_util_dealloc(ct);
	tiifhe_util_dealloc(csw);
	tiifhe_util_dealloc(cswr);
	tiifhe_util_dealloc(delta);
	tiifhe_util_dealloc(p);
}

int
main(void)
{

	tiifhe_bgv_init();
	gmp_randinit_default(state);

	fputs("[+] Testing BGV packing:      ", stderr);
	test_pack();
	fputs("1/1.\n", stderr);

	fputs("[+] Testing BGV encrypt:      ", stderr);
	test_encrypt();
	fputs("1/1.\n", stderr);

	fputs("[+] Testing BGV arithmetic:   ", stderr);
	test_arithmetic();
	fputs("1/1.\n", stderr);

	fputs("[+] Testing BGV modswitch:    ", stderr);
	test_modswitch();
	fputs("3/3.\n", stderr);

	fputs("[+] Testing BGV keyswitch:    ", stderr);
	test_keyswitch();
	fputs("5/5.\n", stderr);

	tiifhe_bgv_deinit();
	gmp_randclear(state);

	return 0;
}
