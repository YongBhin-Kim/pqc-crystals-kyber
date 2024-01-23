# pqc-crystals-kyber
My implementations of NIST PQC(Post-Quantom Cryptography) competition PKE/KEM finalist Crystals-Kyber

**1. lattice-based-crypto** <br>
Crystals-Kyber에 사용되는 격자 기반 문제에 관한 구현 <br>
컴파일 & 빌드 방법 <br>
```
$ make
$ ./poly_mul_example
$ ./polyvec_mul_example
$ ./lwe_example
$ ./module_lwe_example
$ ./module_rlwe_example
```
<br>

**2. kyber512_c** <br>
kyber512(NIST security level 1, 128-bit security)의 C 구현 <br>
컴파일 & 빌드 방법 <br>
```
$ make
$ ./pke
```
<br>

**3. kyber512_py** <br>
kyber512(NIST security level 1, 128-bit security)의 Python 구현 <br>
컴파일 & 빌드 방법 <br>
```
$ python kyber512.py
$ python kyber512_debug.py
```
