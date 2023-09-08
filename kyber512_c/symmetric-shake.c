#include <stddef.h>
#include <stdint.h>
#include "params.h"
#include "fips202.h"
#include "symmetric.h"

/*************************************************
* Name:        kyber_shake128_absorb
*
* Description: Absorb step of the SHAKE128 specialized for the Kyber context.
*
* Arguments:   - keccak_state *state: pointer to (uninitialized) output
*                                     Keccak state
*              - const uint8_t *seed: pointer to SYMBYTES input
*                                     to be absorbed into state
*              - uint8_t i            additional byte of input
*              - uint8_t j            additional byte of input

XOF function : Sample matrix 'A' and 'Tr(A)'
**************************************************/
void kyber_shake128_absorb(keccak_state *state,
                           const uint8_t seed[SYMBYTES],
                           uint8_t x,
                           uint8_t y)
{
  unsigned int i;
  uint8_t extseed[SYMBYTES+2];

  for(i=0;i<SYMBYTES;i++)
    extseed[i] = seed[i];
  extseed[i++] = x;
  extseed[i]   = y;

  shake128_absorb(state, extseed, sizeof(extseed));
}

/*************************************************
* Name:        kyber_shake256_prf
*
* Description: Usage of SHAKE256 as a PRF, concatenates secret and public input
*              and then generates outlen bytes of SHAKE256 output
*
* Arguments:   - uint8_t *out:       pointer to output
*              - size_t outlen:      number of requested output bytes
*              - const uint8_t *key: pointer to the key
*                                    (of length SYMBYTES)
*              - uint8_t nonce:      single-byte nonce (public PRF input)

PRF function : Sample vector 's, e' in keygen part and 'r, e1, e2' in enc part
**************************************************/
void kyber_shake256_prf(uint8_t *out,
                        size_t outlen,
                        const uint8_t key[SYMBYTES],
                        uint8_t nonce)
{
  unsigned int i;
  uint8_t extkey[SYMBYTES+1];

  for(i=0;i<SYMBYTES;i++)
    extkey[i] = key[i];
  extkey[i] = nonce;

  shake256(out, outlen, extkey, sizeof(extkey));
}