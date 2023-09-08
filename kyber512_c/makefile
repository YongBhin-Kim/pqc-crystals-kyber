# Makefile for
# 1. Kyber512 PKE


CC						= gcc
LDFLAGS					= -lcrypto
LDFLAGS					+="-L/usr/local/opt/openssl@3/lib"
CFLAGS					= "-I/usr/local/opt/openssl@3/include"


# Kyber
PKE	 					= pke
PKE_SRCS				= pke.c
PKE_HEADERS				= pke.h

POLY_SRCS 	 			= poly.c
POLY_HEADERS 	 		= poly.h

POLYVEC_SRCS 			= polyvec.c 
POLYVEC_HEADERS 		= polyvec.h

CBD_SRCS				= cbd.c 
CBD_HEADERS				= cbd.h 

AES256_CTR_SRCS			= fips202.c
AES256_CTR_HEADERS		= fips202.h

SHA3_512_SRCS			= rng.c
SHA3_512_HEADERS		= rng.h

REDUCE_SRCS				= reduce.c
REDUCE_HEADERS			= reduce.h

PARAM_HEADERS			= params.h

SHAKE_SRCS				= symmetric-shake.c
SHAKE_HEADERS			= symmetric.h

PRINT_SRCS				= print.c
PRINT_HEADERS			= print.h

SOURCES = $(PKE_SRCS) $(POLY_SRCS) $(POLYVEC) $(POLYVEC_SRCS) $(AES256_CTR_SRCS) $(SHA3_512_SRCS) $(REDUCE_SRCS) $(SHAKE_SRCS) $(CBD_SRCS) 						$(PRINT_SRCS)
HEADERS = $(POLY_HEADERS) $(POLYVEC_HEADERS) $(AES256_CTR_HEADERS) $(SHA3_512_HEADERS) $(REDUCE_HEADERS) $(PARAM_HEADERS) $(SHAKE_HEADERS)	$(CBD_HEADERS)		$(PRINT_HEADERS) $(PKE_HEADERS)

all : $(PKE)

$(PKE) : $(SOURCES) $(HEADERS)
	$(CC) -o $@ $(SOURCES) $(LDFLAGS) $(CFLAGS)

clean:
	rm -rf $(PKE)

# EOF