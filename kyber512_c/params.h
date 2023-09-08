#ifndef PARAMS_H
#define PARAMS_H

#ifndef K
#define K                       2
#endif
#define N                       256
#define Q                       3329

#define SYMBYTES                32

#define POLYBYTES               384 // After encoding
#define POLYVECBYTES            384 * K

#if K == 2
#define ETA1                    3 // Kyber512 ETA1
#else
#define ETA1                    2 // Kyber768/1024 ETA1
#endif
#define ETA2                    2 // Kyber512/768/1024 ETA2

#define POLYCOMPRESSEDBYTES     128 // After compressing, encoding (d_v = 4)
#define POLYVECCOMPRESSEDBYTES  320 * K // After compressing, encoding (d_u = 10)

#define MSGBYTES                32  // 32-bytes message space 
#define PUBLICKEYBYTES          POLYVECBYTES + SYMBYTES // (Compress of polyvec 'b' || 32-bytes array 'publicseed')
#define SECRETKEYBYTES          POLYVECBYTES // (Compress of polyvec 's')
#define CIPHERTEXTBYTES         POLYVECCOMPRESSEDBYTES + POLYCOMPRESSEDBYTES // (ct_u after compress, encode || ct_v after Compress, Encode) 

#endif