#include "tiifhe.h"

void
msg_print(const char *prefix, const tiifhe_Message *m)
{
	size_t j;

	printf("%s: ", prefix);
	for (j = 0; j < 9; ++j)
		gmp_printf("%Zd, ", m->value[j]);
	gmp_printf("%Zd, ...\n", m->value[j]);
}

int
main(void)
{
	/* Example: Polynomial Evaluation
	 *
	 * In the following, we will evaluate the polynomial
	 *
	 *   x^2 + 2 * rot(x, 1) + 3
	 *
	 * on a packed message to showcase how to use the high-level API of
	 * the library. If following the example is confusing during any part,
	 * please open an issue so we can improve its documentation. We start
	 * by declaring all variables used in the process.
	 */
	tiifhe_KeySecret *sk;
	tiifhe_KeyPublic *pk;
	tiifhe_KeySwitch *kmul, *krot;
	tiifhe_Ciphertext *ct0, *ct1;
	tiifhe_Message *m, *c;
	size_t j;

	/* Initialization
	 *
	 * First, we need to initialize the library, precomputing required
	 * values for the build time parameters.
	 */
	tiifhe_bgv_init();

	/* Memory Allocation
	 *
	 * Since we are working with pure C, memory management is manual.
	 * However, for ease of use, we provide allocation functions for all
	 * important data structures checking for allocation errors. The
	 * allocation function for the message is provided by the message
	 * layer with a slightly different prefix; we can also use it to
	 * allocate more than one message if necessary.
	 */
	sk = tiifhe_bgv_alloc_sk();
	pk = tiifhe_bgv_alloc_pk();
	ct0 = tiifhe_bgv_alloc_ct();
	ct1 = tiifhe_bgv_alloc_ct();
	kmul = tiifhe_bgv_alloc_ksw();
	krot = tiifhe_bgv_alloc_ksw();
	m = tiifhe_msg_alloc(1);
	c = tiifhe_msg_alloc(1);

	/* Key Generation
	 *
	 * Before we can start with the homomorphic computation, we have to
	 * generate keys. First, we generate a secret key and a corresponding
	 * public key. For multiplications, we require a key switching key
	 * as well as for every rotation we want to perform.
	 */
	tiifhe_bgv_keygen_secret(sk);
	tiifhe_bgv_keygen_public(pk, sk);
	tiifhe_bgv_keygen_switch_mul(kmul, sk);
	tiifhe_bgv_keygen_switch_rot(krot, sk, 1);

	/* Message Initialization
	 *
	 * For the mpz message layer, we can use the existing mpz functions
	 * to set the message to any desired values. Since we want to operate
	 * on multiple values in parallel (on TIIFHE_N values to be precise),
	 * we also pack the message.
	 *
	 * You can think of the packed message as a 2xN/2-array, that is two
	 * row vectors:
	 *   [  0,       1,       2, ..., N/2 - 1]
	 *   [N/2, N/2 + 1, N/2 + 2, ...,   N - 1]
	 *
	 * Arithmetic operations are performed slot-wise while rotations are
	 * applied row-wise, that is each row is rotated as specified.
	 */
	for (j = 0; j < TIIFHE_N; ++j)
		mpz_set_ui(m->value[j], j);
	msg_print("msg", m);
	tiifhe_msg_pack(m, m);

	/* Encryption
	 *
	 * We now create two ciphertexts encrypting the same underlying
	 * message. Using the copy function is cheaper than performing
	 * another encryption.
	 */
	tiifhe_bgv_encrypt(ct0, pk, m);
	tiifhe_bgv_copy(ct1, ct0);

	/* Decryption
	 *
	 * Here, we quickly check that the above worked. We can always use the
	 * decryption function with the secret key to encrypt a ciphertext.
	 * Remember to unpack a message if you packed it before encryption.
	 */
	tiifhe_bgv_decrypt(m, sk, ct0);
	tiifhe_msg_unpack(m, m);
	msg_print("ct0", m);

	tiifhe_bgv_decrypt(m, sk, ct1);
	tiifhe_msg_unpack(m, m);
	msg_print("ct1", m);

	/* Multiplication
	 *
	 * We proceed with the most complex part of the polynomial evaluation,
	 * the multiplication. After each multiplication, we perform a key
	 * switching to reduce the ciphertext from degree three back to a
	 * ciphertext of degree two. Note that our decryption function can
	 * also handle ciphertexts of degree three, so if you do not perform
	 * any further operations with it, key switching is optional.
	 *
	 * We also perform a scaling which is recommended due to the large
	 * error growth during evaluation. Usually, scaling is performed
	 * after every multiplication and parameters have to be chosen
	 * accordingly. However, one can perform scaling at any point in the
	 * homomorphic circuit. We start removing a modulus at the lowest
	 * index that is still available. In the following, this corresponds
	 * to scaling away the modulus at index 0 as specified at build time
	 * since we are using a fresh encryption.
	 *
	 * A different example would be a ciphertext which already dropped the
	 * primes at index 0 and index 1. Then, a call to the scaling function
	 * would remove the prime at index 2.
	 *
	 * We check that the multiplication did indeed happen as expected,
	 * printing out the encrypted vector with squared values.
	 */
	tiifhe_bgv_mul(ct0, ct0, ct0);
	tiifhe_bgv_keyswitch_mul(ct0, kmul);
	tiifhe_bgv_modswitch(ct0);

	tiifhe_bgv_decrypt(m, sk, ct0);
	tiifhe_msg_unpack(m, m);
	msg_print("mul", m);

	/* Rotation & Constant Multiplication
	 *
	 * We now perform the rotation which in itself does not increase the
	 * error. However, in contrast to multiplication, we do need a key
	 * switching even if we decrypt instantly afterward and slightly
	 * increases the error.
	 *
	 * For the constant multiplication, we prepare a message as if we
	 * would encrypt it, in this case setting all of our message slots
	 * to the value two. Then, we pack the message as before.
	 */
	tiifhe_bgv_rot(ct1, ct1, 1);
	tiifhe_bgv_keyswitch_rot(ct1, krot);

	for (j = 0; j < TIIFHE_N; ++j)
		mpz_set_ui(c->value[j], 2);
	tiifhe_msg_pack(c, c);
	tiifhe_bgv_mulc(ct1, ct1, c);

	tiifhe_bgv_decrypt(m, sk, ct1);
	tiifhe_msg_unpack(m, m);
	msg_print("rot", m);

	/* Addition
	 *
	 * We are now almost ready to add the two ciphertext and a constant
	 * to compute the final result. However, there is one issue: our two
	 * ciphertexts are not at the same level due to the scaling which we
	 * performed after multiplication. To fix this issue, we can drop a
	 * ciphertext to the level of another one.
	 *
	 * For the constant, we proceed as with the constant multiplication,
	 * setting all slots of the constant message to three. Finally, we
	 * add up the ciphertexts and add the constant to obtain our final
	 * result.
	 */
	tiifhe_bgv_drop(ct1, ct0);
	tiifhe_bgv_add(ct0, ct0, ct1);

	for (j = 0; j < TIIFHE_N; ++j)
		mpz_set_ui(c->value[j], 3);
	tiifhe_msg_pack(c, c);
	tiifhe_bgv_addc(ct0, ct0, c);

	tiifhe_bgv_decrypt(m, sk, ct0);
	tiifhe_msg_unpack(m, m);
	msg_print("add", m);

	/* Resource Cleanup
	 *
	 * In this example, cleaning up resources is optional since we exit
	 * the program afterward anyway. However, in more complex scenarios,
	 * cleaning up your resources might be a good idea especially with
	 * larger parameters.
	 */
	tiifhe_bgv_dealloc(sk);
	tiifhe_bgv_dealloc(pk);
	tiifhe_bgv_dealloc(ct0);
	tiifhe_bgv_dealloc(ct1);
	tiifhe_bgv_dealloc(kmul);
	tiifhe_bgv_dealloc(krot);
	tiifhe_msg_dealloc(m, 1);
	tiifhe_msg_dealloc(c, 1);

	tiifhe_bgv_deinit();

	return 0;
}
