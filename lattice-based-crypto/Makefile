# Makefile for
# 1. lwe_example // lwe example
# 2. poly_mul // poly_mul example (over ring)
# 3. polyvec_mul // polyvec_mul example (over ring)
# 4. module_rlwe_example // module-RLWE example (over module-ring)


CC						= gcc

# LWE example
LWE			 			= lwe_example
LWE_SRCS	 			= lwe_example.c

# Polynomial multiplications - poly * poly (element in R_q)
# - schoolbook multiplication ( O(n^2) )
POLYMUL 	 			= poly_mul_example
POLYMUL_SRCS 			= poly_mul.c poly_mul_example.c
POLYMUL_HEADERS 		= poly_mul.h

# Polynomial multiplications - polyvec * polyvec and matrix * polyvec (element in R_q)
POLYVECMUL				= polyvec_mul_example
POLYVECMUL_SRCS			= polyvec_mul_example.c poly_mul.c
POLYVECMUL_HEADERS		= polyvec_mul_example.h poly_mul.h

# RLWE example
MODULE_LWE 					= module_lwe_example
MODULE_LWE_SRCS 			= module_lwe_example.c poly_mul.c
MODULE_LWE_HEADERS 			= module_lwe_example.h poly_mul.h

# module-RLWE example
MODULE_RLWE 					= module_rlwe_example
MODULE_RLWE_SRCS 				= module_rlwe_example.c poly_mul.c
MODULE_RLWE_HEADERS 			= module_rlwe_example.h poly_mul.h


all : $(LWE) $(POLYMUL) $(POLYVECMUL) $(MODULE_LWE) $(MODULE_RLWE)

$(LWE): $(LWE_SRCS) # $(UTIL_SRCS) $(UTIL_HEADERS)
	$(CC) -o $@ $(LWE_SRCS) $(LDFLAGS) 

$(POLYMUL): $(POLYMUL_SRCS) $(POLYMUL_HEADERS) 
	$(CC) -o $@ $(POLYMUL_SRCS) $(LDFLAGS) 

$(POLYVECMUL) : $(POLYVECMUL_SRCS) $(POLYVECMUL_HEADERS) $(POLYMUL_SRCS) $(POLYMUL_HEADERS)
	$(CC) -o $@ $(POLYVECMUL_SRCS) $(LDFLAGS)

$(MODULE_LWE) : $(MODULE_LWE_SRCS) $(MODULE_LWE_HEADERS) $(POLYMUL_SRCS) $(POLYMUL_HEADERS)
	$(CC) -o $@ $(MODULE_LWE_SRCS) $(LDFLAGS)

$(MODULE_RLWE) : $(MODULE_RLWE_SRCS) $(MODULE_RLWE_HEADERS) $(POLYMUL_SRCS) $(POLYMUL_HEADERS)
	$(CC) -o $@ $(MODULE_RLWE_SRCS) $(LDFLAGS)


clean:
	rm -rf $(LWE) 
	rm -rf $(POLYMUL)
	rm -rf $(MODULE_LWE)
	rm -rf $(MODULE_RLWE)

# EOF