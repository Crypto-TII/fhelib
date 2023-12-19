#ifndef TIIFHE_MSG_H
#define TIIFHE_MSG_H

#include <assert.h>
#include <gmp.h>
#include <stdio.h>

#include "config.h"

typedef struct tiifhe_message     tiifhe_Message;
typedef struct tiifhe_message_mod tiifhe_MessageMod;

struct tiifhe_message {
	mpz_t value[TIIFHE_N];
};

struct tiifhe_message_mod {
	mpz_t value;
	mpz_t root;
};

extern tiifhe_MessageMod tiifhe_t;

tiifhe_Message *tiifhe_msg_alloc(size_t len);
void tiifhe_msg_dealloc(tiifhe_Message *m, size_t len);

void tiifhe_msg_addc(tiifhe_Message *rop, const tiifhe_Message *m, uint64_t c);
void tiifhe_msg_cmod(tiifhe_Message *rop, const tiifhe_Message *m, const mpz_t mod);
void tiifhe_msg_copy(tiifhe_Message *rop, const tiifhe_Message *m);
void tiifhe_msg_crt_i64(tiifhe_Message *m, const int64_t **rns, const int64_t *mods, size_t len);
void tiifhe_msg_crt_u64(tiifhe_Message *m, const uint64_t **rns, const uint64_t *mods, size_t len);
void tiifhe_msg_fprintdiv(FILE *f, const tiifhe_Message *m, const mpz_t mod);
void tiifhe_msg_mulc(tiifhe_Message *rop, const tiifhe_Message *m, uint64_t c);
void tiifhe_msg_mulc_inv(tiifhe_Message *rop, const tiifhe_Message *m, uint64_t c);
void tiifhe_msg_pack(tiifhe_Message *rop, const tiifhe_Message *m);
void tiifhe_msg_pmod(tiifhe_Message *rop, const tiifhe_Message *m, const mpz_t mod);
void tiifhe_msg_unpack(tiifhe_Message *rop, const tiifhe_Message *m);

void tiifhe_msgmod_deinit(void);
void tiifhe_msgmod_init(void);

#endif /* TIIFHE_MSG_H */
