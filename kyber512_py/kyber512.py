import random
import hashlib

########################################################################################################################
# Parameter sets
########################################################################################################################
"""
    Kyber parameters
"""
global Q, N, K, ETA1, ETA2
K                      = 2
N                      = 256
Q                      = 3329

global SYMBYTES
SYMBYTES               = 32 

global POLYBYTES, POLYVECBYTES
POLYBYTES              = 384                                     # 512 *12/8 (By function 'Encode_12')
POLYVECBYTES           = K * POLYBYTES                           # Depended on strength of Kyber

global POLYCOMPRESSEDBYTES, POLYVECCOMPRESSEDBYTES # length of ct, length of ct1, length of ct2
if K == 2:
    ETA1                   = 3
    POLYCOMPRESSEDBYTES    = 128
    POLYVECCOMPRESSEDBYTES = K * 320
elif K == 3:
    ETA1                   = 2
    POLYCOMPRESSEDBYTES    = 128
    POLYVECCOMPRESSEDBYTES = K * 320
    
elif K == 4:
    POLYCOMPRESSEDBYTES    = 160
    POLYVECCOMPRESSEDBYTES = K * 352
ETA2                       = 2


global MSGBYTES, PUBLICKEYBYTES, SECRETKEYBYTES
MSGBYTES                   = SYMBYTES
PUBLICKEYBYTES             = POLYVECBYTES + SYMBYTES
SECRETKEYBYTES             = POLYVECBYTES
    
CIPHERTEXTBYTES            = POLYVECCOMPRESSEDBYTES + POLYCOMPRESSEDBYTES


"""
    Debug options
"""
global DEBUGTYPE, HEX, DECIMAL
HEX         = 16
DECIMAL     = 10
DEBUGTYPE   = DECIMAL

""" 
    Define type of message
        - If you want to use an already defined test vector as a message, define this value(MESSAGETYPE) as TESTVECTOR, 
          and if you want a message as standard input, define this value(MESSAGETYPE) as STDIN
"""
global MESSAGETYPE, TESTVECTOR, STDIN
STDIN       = 1
TESTVECTOR  = 2
MESSAGETYPE = TESTVECTOR


########################################################################################################################
# Functions
########################################################################################################################

"""
    For cbd1
    (24-bit integer <- 3-bytes array)
"""
def load24_littleendian(state0, state1, state2):
    return (state2<<16) | (state1<<8) | (state0)

"""
    For cbd2
    (32-bit integer <- 4-bytes array)
"""
def load32_littleendian(state0, state1, state2, state3):
    return (state3<<24) | (state2<<16) | (state1<<8) | (state0)
    
"""
    Sample polynomial vector
    
    'res'      <- prf(seed in B^32, nonce) (using shake256)
    (192-bytes <- 32-bytes seed)
"""
def prf(outlen, seed, nonce):
    
    strseed = "".join(format(byte, '08b') for byte in seed).encode() # str <- bytes
    strseed += format(nonce, '08b').encode() # concate nonce 
    
    state = hashlib.new('shake_256') # Using python hashlib module
    state.update(strseed)
    
    state = state.hexdigest(outlen)
    state = bytearray.fromhex(state)
    
    # res in B^128 (or B^192) <- prf(seed in B^32, nonce)
    res = [0]*outlen
    for i in range(outlen):
        res[i] = state[i]

    return res
    
""" 
    Sample a polynomial vector 'e1', 'e2', etc.
    (Case : ETA == 2)
    (duplicate process : poly <- 128-byte array)
"""
def cbd2(buf):
    res = poly()
    
    for i in range(N//8):
        t = load32_littleendian(buf[4*i], buf[4*i+1], buf[4*i+2], buf[4*i+3])       #   uint32_t t,d;
        d = t & 0x55555555                                                          #   int16_t a,b;
        d += (t>>1) & 0x55555555
        for j in range(8):
            a = (d >> (4*j+0)) & 0x3
            b = (d >> (4*j+2)) & 0x3
            res.coeffs[8*i+j] = a - b
    return res

""" 
    Sample a polynomial vector 's', 'e', 'r', etc.
    (Case : ETA == 3)
    (duplicate process : poly <- 192-byte array)
"""
def cbd1(buf):
    res = poly()
    for i in range(N//4):
        t = load24_littleendian(buf[3*i], buf[3*i+1], buf[3*i+2])
        d = t & 0x00249249
        d += (t>>1) & 0x00249249
        d += (t>>2) & 0x00249249
        
        for j in range(4):
            a = (d >> (6*j+0)) & 0x7
            b = (d >> (6*j+3)) & 0x7
            res.coeffs[4*i+j] = a - b
    return res

"""
    Sample a matrix 'A'
    (Kyber function : rej_uniform)
    
    Parameters :
                - coeffs : output which is element of R_Q (entry of matrix 'A')
                - len    : length want to sample element (use fixed N)
                - buf    : input buffer used to sampling
                - buflen : length of input buffer
"""
def parse(coeffs, len, buf, buflen): # mat.vec.coeffs, N, buf[506](xof의출력), 506
    ctr = pos = 0
    while ctr < len and pos + 3 <= buflen:
        val0 = ((buf[pos+0] >> 0) | (buf[pos+1] << 8)) & 0xFFF # 12비트에서 샘플링
        val1 = ((buf[pos+1] >> 4) | (buf[pos+2] << 4)) & 0xFFF
        pos += 3
        
        
        if val0 < Q:
            coeffs[ctr] = val0
            ctr+=1
        if ctr < len and val1 < Q:
            coeffs[ctr] = val1
            ctr+=1
    return ctr

"""
    For sample a matrix 'A'
    (Extend function)
"""
def xof(seed, x, y):
    
    extseed = "".join(format(byte, '08b') for byte in seed).encode() # str <- bytes
    extseed += format(x, '08b').encode()
    extseed += format(y, '08b').encode()
    
    state = hashlib.new('shake_128') # Using python hashlib module
    state.update(extseed)
    
    state = state.hexdigest(506) # fixed xof_state structure : len = 25 
    state = bytearray.fromhex(state)
    
    res = [0]*506
    for i in range(506):
        res[i] = state[i]
        
    return res

"""
    Class : Print matrix, vector and element
    ( Print array[ARRAYSIZE], 
      polymat == polyvec[K], 
      polyvec == poly[K], 
      poly == coeffs[N] )
"""
class printf:
        
    @staticmethod
    def printarray(str, state):
        if str != None:
            print(' * %s \n   [ ' %str, end=' ')
        
        if DEBUGTYPE == HEX:
            for i in range(len(state)):
                print('%02x' %state[i], end=' ')
                if len(state) > 32 and i % 32 == 31:
                    print('\n', end='      ')
        elif DEBUGTYPE == DECIMAL:
            for i in range(len(state)):
                print('%03d' %state[i], end=' ')
                if len(state) > 32 and i % 32 == 31:
                    print('\n', end='      ')
                
        print(']')
    
    @staticmethod
    def printpoly(str, state):
        if str != None:
            print(' * %s \n   [ ' %str, end=' ')
        else:
            print('   ', end='')
        
        if DEBUGTYPE == HEX:
            for i in range(N):
                print('%3x ' %state.coeffs[i], end=' ')
                if i % 16 == 15:
                    print('\n   ', end='')
        elif DEBUGTYPE == DECIMAL:
            for i in range(N):
                print('%4d ' %state.coeffs[i], end=' ')
                if i % 16 == 15:
                    print('\n   ', end='')
                
                
        print(']')
    
    @staticmethod
    def printvector(str, state):
        if str != None:
            print(' * %s \n   [ ' %str, end=' ')
        
        for i in range(K):
            printf.printpoly(None, state.vec[i])
        print('   ]')
    
    @staticmethod
    def printmatrix(str, state):
        if str != None:
            print(' * %s \n   [ ' %str, end=' ')
            
        for i in range(K):
            printf.printvector(None, state.mat[i])
        print(']')

"""
    Get random byte array
    (publicseed, noiseseed, coin which is used in sampling of {'A',   's, e',   'e1, e2, r'})
"""
class randombytes:

    def get_randomarray():
        return [random.randint(0, 256) for _ in range(SYMBYTES)] # Using python random module
    
"""
    Treating byte array : concate two byte arrays
"""
class byte:
    
    def concate(a, b):
        res = a[:]
        
        for item in b:
            res.append(item)
        
        return res

"""
    Define poly (== R_Q) and operations
"""
class poly:
    
    def __init__(self):
        self.Q = Q
        self.N = N
        self.coeffs = [0 for _ in range(N)] # Representation of polynomial is signed 16-bit array in Kyber 

    """
        Operator overloading : polynomial addition over R_Q
    """
    def __add__(self, other):
        res = poly()
        for i in range(self.N):
            res.coeffs[i] = (self.coeffs[i] + other.coeffs[i]) % self.Q
        return res
    
    """
        Operator overloading : polynomial subtraction over R_Q
    """
    def __sub__(self, other):
        res = poly()
        for i in range(self.N):
            if self.coeffs[i] >= other.coeffs[i]:
                res.coeffs[i] = self.coeffs[i] - other.coeffs[i]
            else:
                res.coeffs[i] = self.coeffs[i] + self.Q - other.coeffs[i]
        return res

    """
        Operator overloading : polynomial multiplication over R_Q
    """
    def __mul__(self, other):
        res = poly()
        tempcoef = 0
        temppoly = [0 for _ in range(2*N)]
        
        for i in range(N):
            for j in range(N):
                tempcoef = self.coeffs[i] * other.coeffs[j]
                tempcoef = tempcoef + temppoly[i + j]
                temppoly[i + j] = tempcoef % Q
        for i in range(2*N):
            if temppoly[i-N] >= temppoly[i]:
                res.coeffs[i-N] = temppoly[i-N] - temppoly[i]
            else:
                res.coeffs[i-N] = temppoly[i-N] + Q - temppoly[i]
        return res
    
    """
        Modulor reduction for coefficients over Z_Q
    """
    def reduce(self):
        for i in range(self.N):
            self.coeffs[i] %= Q

    """
        Function : Sampling
            Used to sample polynomial depended on ETA1 or ETA2   (use cbd1 or cbd2)
    """
    @staticmethod
    def sample(seed, nonce, eta):
        if eta == 2:
            return cbd2(prf(eta*N//4, seed, nonce)) # outlen = 128, nonce = 0, 1
        elif eta == 3:
            return cbd1(prf(eta*N//4, seed, nonce)) # outlen = 128, nonce = 0, 1
        
    """
        Decompress message 'm'
        'res' (poly) <- byte array 'm' (128-bits) 
        
        If   message bit is '0'
            then coeffs is '0'
        elif message bit is '1' 
            then coeffs is 'rounding(Q/2)' (= 1665)
        
    """
    @staticmethod
    def frommsg(m):
        res = poly()
        
        for i in range(N//8):
            for j in range(8):
                if (m[i]>>j)&1 == 1:
                    res.coeffs[8*i+j] = 1665
        return res
    
    """
        Function : poly.tomsg
        
        Used in compress 'ct_v - Tr(s)*ct_u' which is decryption scheme in Kyber 
            (It means decrypt ciphertext) (recover process!)
        
        If 'ct_v - Tr(s)*ct_u' is close to 0 than 1665
            then message bit is 0
        elif   'ct_v - Tr(s)*ct_u' is close to 1665 than 0
            then message bit is 1
        
    """
    @staticmethod
    def tomsg(p):
        res = [0]*MSGBYTES
        
        for i in range(N//8):
            for j in range(8):
                t = (((p.coeffs[8*i+j]<<1) + Q//2) //Q) &1
                res[i] |= t << j
        return res
    
    """
        Function 'Encode_12' 
            used to pack public key 'pk' and secret key 'sk' forms of entry of polynomial vector 'b', 's'
    """
    @staticmethod
    def tobytes(p):
        res = [0]*((N * 3) // 2) # 256 * (3/2)
        for i in range(N//2):
            t0 = p.coeffs[2*i] % Q
            t1 = p.coeffs[2*i+1] % Q
            res[3*i+0] = (t0 >> 0) & 0xff # 
            res[3*i+1] = ((t0 >> 8) & 0x0f) | ((t1 << 4) & 0xf0)
            res[3*i+2] = (t1 >> 4) & 0xff
        
        # p03 p02 p01 p00    || p13 p12 p11 p10 (each 4-bits)
        # 
        #    p01 p00    || p10 p02 ||     p12 p11 (= output byte array)
        
        return res

    """
        Function 'Decode_12' 
            used to unpack public key 'b' and secret key 's' forms of entry of byte array 'pk', 'sk'
    """
    @staticmethod  
    def frombytes(a):
        res = poly()
        for i in range(N//2):
            res.coeffs[2*i] = (a[3*i] | ((a[3*i+1]<<8)&0xf00))&0xfff  ; res.coeffs[2*i] %= Q
            res.coeffs[2*i+1] = ((a[3*i+1]>>4) | (a[3*i+2]<<4))&0xfff ; res.coeffs[2*i+1] %= Q
        return res
    
    """
        Function poly.compress
        
        Used in compress 'ct_v'
        byte array 'res' (128-bytes) <- poly 'ct_v' (512-bytes) (in Kyber512)
        
        c2 <- Encode_4( Compress_Q(ct_v, 4) )                   (for Kyber512, 'd_v' is 4)
            (4-bits entropy <- 12-bits entropy)
    """
    @staticmethod 
    def compress(p):
        res = [0]*POLYCOMPRESSEDBYTES

        t = [0]*8
        cnt = 0
        for i in range(32): # 
            for j in range(8):
                t[j] = (((p.coeffs[8*i+j]<<4) + Q//2)//Q) & 15
            res[cnt]   = t[0] | (t[1] << 4)
            res[cnt+1] = t[2] | (t[3] << 4)
            res[cnt+2] = t[4] | (t[5] << 4)
            res[cnt+3] = t[6] | (t[7] << 4)
            cnt+=4

        return res
    
    """
        Function poly.decompress
        
        Used in decompress ct_v
        poly ct_v (512-bytes) <- byte array (128-bytes) (in Kyber512)
        
        ct_v <- Decompress_Q( Decode_4(c2, 4) )         (for Kyber512, 'd_v' is 4)
    """
    @staticmethod 
    def decompress(a):
        res = poly()
        cnt = 0
        for i in range(128):
            res.coeffs[2*i]   = ((a[i] & 15)*Q + 8) >> 4
            res.coeffs[2*i+1] = ((a[i] >> 4)*Q + 8) >> 4
        return res

"""
    Define polyvec (== R_Q^K) and operations
"""
class polyvec:
    def __init__(self):
        self.Q = Q
        self.N = N
        self.K = K
        self.vec = [poly() for _ in range(K)]   # Representation of vector is polynomial array

    """
        Operator overloading : polynomial vector addition over R_Q
    """
    def __add__(self, other):
        res = polyvec()
        for i in range(self.K):
            res.vec[i] = self.vec[i] + other.vec[i]
        return res

    """
        Operator overloading : polynomial vector subtraction over R_Q
    """
    def __sub__(self, other):
        res = polyvec()
        for i in range(self.K):
            res.vec[i] = self.vec[i] - other.vec[i]
        return res

    """
        Operator overloading : polynomial vector multiplication over R_Q
    """
    def __mul__(self, other):
        res = poly()
        temp = poly()
        for i in range(self.K):
            temp = self.vec[i] * other.vec[i]
            res = res + temp
        return res

    """
        Modulor reduction over polynomial vector R_Q^K is dup polynomial modulo reduction
    """
    def reduce(self):
        for i in range(self.K):
            (self.vec[i]).reduce()

    """
        Function : Sampling
            Used to sample vector depended on ETA1 or ETA2 (cbd1 or cbd2)
    """
    def sample(self, seed, nonce, eta):
        if eta == 2:
            for i in range(K):
                self.vec[i] = cbd2(prf(eta*N//4, seed, nonce)) # outlen = 128, nonce = 0, 1
                nonce += 1
        elif eta == 3:
            for i in range(K):
                self.vec[i] = cbd1(prf(eta*N//4, seed, nonce)) # outlen = 192, nonce = 0, 1
                nonce += 1

    """
        Function 'Encode_12' 
            used to pack public key 'pk' and secret key 'sk' forms of entry of polynomial vector 'b', 's'
    """
    @staticmethod
    def tobytes(pv): # byte array[384+384] <- polyvec
        res = [0]*(POLYVECBYTES)
        resarray = []
        for i in range(K):
            resarray.append(poly.tobytes(pv.vec[i]))
        for i in range(K):
            for j in range(POLYBYTES):
                res[POLYBYTES*i+j] = resarray[i][j]
        return res

    """
        Function 'Decode_12' 
            used to unpack public key 'b' and secret key 's' forms of entry of byte array 'pk', 'sk'
    """
    @staticmethod 
    def frombytes(a): # polyvec <- byte array
        temp = [0 for _ in range(POLYBYTES)]*K
        for i in range(K):
            temp[i] = a[i*(POLYBYTES) : (i+1)*(POLYBYTES)] # temp[i] = b[i*384 : i*384 + 384]
        
        res = polyvec()
        for i in range(K):
            res.vec[i] = poly.frombytes(temp[i])
        return res
    
    """
        Function polyvec.compress
        
        Used in compress ct_u
        byte array C1 (640-bytes) <- polyvec ct_u (1024-bytes)
        
        C1 <- Encode_d_u( Compress_Q( ct_u, d_u ) )
            (d_u-bits entropy <- 12-bits entropy)
    """
    @staticmethod 
    def compress(ct_u):
        res = [0]*POLYVECCOMPRESSEDBYTES
        t = [0]*4
        cnt = 0
        for i in range(K):
            for j in range(N//4):
                for k in range(4):
                    t[k] = ((ct_u.vec[i].coeffs[4*j+k] << 10) + Q//2)//Q & 0x3ff 
                    
                res[cnt]   = (t[0] >> 0) & 0xff
                res[cnt+1] = ((t[0] >> 8) | (t[1] << 2)) & 0xff
                res[cnt+2] = ((t[1] >> 6) | (t[2] << 4)) & 0xff
                res[cnt+3] = ((t[2] >> 4) | (t[3] << 6)) & 0xff
                res[cnt+4] = (t[3] >> 2) & 0xff
                cnt += 5
        return res
    
    """
        Function polyvec.decompress
    
        Used in decompress C1
        polyvec ct_u (1024-bytes) <- byte array C1 (640-bytes in kyber512) 
        
        ct_u <- Decompress_Q( Decode_10(C1), 10 ) (kyber512 : d_u = 10)
    """
    @staticmethod
    def decompress(c1):
        res = polyvec()
        t = [0]*4
        cnt = 0
        for i in range(K):
            for j in range(N//4):
                t[0] = (c1[cnt+0] >> 0) | (c1[cnt+1] << 8)
                t[1] = (c1[cnt+1] >> 2) | (c1[cnt+2] << 6)
                t[2] = (c1[cnt+2] >> 4) | (c1[cnt+3] << 4)
                t[3] = (c1[cnt+3] >> 6) | (c1[cnt+4] << 2)
                cnt += 5
                
                for k in range(4):
                    res.vec[i].coeffs[4*j+k] = ((t[k] & 0x3FF)*Q + 512) >> 10
            
        return res
    
"""
    Define polymat (== R_Q^(K*K)) and operations
"""
class polymat:
    def __init__(self):
        self.Q = Q
        self.N = N
        self.K = K
        self.mat = [polyvec() for _ in range(K)] # Representation of matrix is vector array

    """
        Operator overloading : multiplication of matrix and vector over R_Q
    """
    def __mul__(self, other): # matrix * vector
        res = polyvec()
        for i in range(K):
            res.vec[i] = self.mat[i] * other
        return res
    
    """
        Finction : transpose matrix (unused)
    """
    @staticmethod
    def trans(A):
        res = polymat()
        for i in range(K):
            for j in range(K):
                res.mat[i].vec[j].coeffs[:] = A.mat[j].vec[i].coeffs[:]
        return res

    """
        Finction : transpose matrix (unused)
    """
    def transmatrix(self):
        res = polymat()
        for i in range(K):
            for j in range(K):
                res.mat[i].vec[j].coeffs[:] = self.mat[j].vec[i].coeffs[:]
        return res
    
    """
        Function : Sampling (= rej_uniform & XOF)
            Used to sample matrix 'A' or 'Tr(A)' which entry is uniformly reduction
            
            Parameter 'transpose' is 0 : do not transpose 
            Parameter 'transpose' is 1 : do transpose
    """
    def sample(self, seed, transpose):
            
        for i in range(K):
            for j in range(K):
                if transpose == True:                                         
                    buf = xof(seed, j, i)                                     # xof : SHAKE-128
                else:
                    buf = xof(seed, i, j)                                     # xof : SHAKE-128
                ctr = parse(self.mat[i].vec[j].coeffs, self.N, buf, len(buf)) # function parse : rejection sampling
                
                if ctr != self.N:
                    self.sample(self, seed)                                   # Difference to kyber reference code : Kyber에서는 xof(keccak_squeeze)함수 반복
            
        
########################################################################################################################
# Kyber class
########################################################################################################################

"""
    Kyber512 class
"""
class kyber512:
    def __init__(self, K, Q, N, ETA1, ETA2):
        self.K = K
        self.Q = Q
        self.N = N
    
    # 1. Key generation part
    def keygen(pk, sk):
        nonce = 0                                                   # Declare nonce used in vector sampling
        publicseed = randombytes.get_randomarray()                  # Get 32-bytes uniformly random array
        noiseseed = randombytes.get_randomarray()                   # Get 32-bytes uniformly random array
        A, s, e = polymat(), polyvec(), polyvec()                   # Declare public matrix 'A', secret key vector 's', secret vector 'e'

        A.sample(publicseed, False)                                 # Sample matrix 'A' uniformly random
        s.sample(noiseseed, nonce, ETA1)                            # Sample vector 's' cbd with ETA1 and nonce
        nonce += K                      
        e.sample(noiseseed, nonce, ETA1)                            # Sample vector 'e' cbd with ETA1 and nonce

        b = A * s + e                                               # Compute public key vector 'b'

        pk_temp = byte.concate(polyvec.tobytes(b), publicseed)      # pk <- ( Encode_12(b) || publicseed )
        for i in range(PUBLICKEYBYTES):
            pk[i] = pk_temp[i] 
        sk_temp = polyvec.tobytes(s)                                # sk <- ( Encode_12(s) )
        for i in range(SECRETKEYBYTES):
            sk[i] = sk_temp[i]
        
    # 2. Encryption part
    def enc(pk, m, coin, ct):               
        nonce = 0                                                   # Declare nonce used in vector sampling
        AT , r, e1, e2 = polymat(), polyvec(), polyvec(), poly()    # Declare parameters used in encryption part
        
        b = polyvec.frombytes(pk[:PUBLICKEYBYTES])                  # b <- Decode_12(pk)
        
        publicseed = pk[POLYVECBYTES:]                              # publicseed <- pk[768:800]
        
        AT.sample(publicseed, True)                                 # Sample matrix 'AT' which is transpose of 'A' in keygen part : True means transpose
        r.sample(coin, nonce, ETA1); nonce += K                     # Sample vector r,  Increase nonce
        e1.sample(coin, nonce, ETA2); nonce += K                    # Sample vector e1, Increase nonce
        e2 = poly.sample(coin, nonce, ETA2)                         # Sample poly e2,   Increase nonce


        mp = poly.frommsg(m)                                        # mp <- Decompress_Q( Decode_1(m), 1 ) 
        ct_u = AT * r + e1                                          # Compute vector forms of first ciphertext 'ct_u'
        ct_v = b * r + e2 + mp                                      # Compute polynomial forms of second ciphertext 'ct_v'
        
        
        ct1 = polyvec.compress(ct_u)                                # c1 <- Encode_d_u( Compress_q(ct_u, d_u) )
        ct2 = poly.compress(ct_v)                                   # c2 <- Encode_d_v( Compress_q(ct_v, d_v) ) ( 128-bytes array <- poly() )
        ct_temp = byte.concate(ct1, ct2)                            # ct <- (c1 || c2)
        
        for i in range(len(ct_temp)):
            ct[i] = ct_temp[i]                     
        
    # 3. Decryption part
    def dec(sk, ct, recovered):
        ct1 = ct[:POLYVECCOMPRESSEDBYTES]                           # ct1 <- ct[0  :640]  (640-bytes in Kyber512/768)
        ct2 = ct[POLYVECCOMPRESSEDBYTES:]                           # ct2 <- ct[640:768]  (128-bytes in Kyber512/768)

        ct_u = polyvec.decompress(ct1)                              # ct_u <- Decompress_Q( Decode_d_u(ct1), d_u )
        ct_v = poly.decompress(ct2)                                 # ct_v <- Decompress_Q( Decode_d_v(ct2), d_v )
        
        s = polyvec.frombytes(sk)                                   # s <- Decode_12(sk)
        
        mp_witherror = ct_v - (s * ct_u)                            # Compute ct_v - s^T * ct_u
        
        recovered_temp = poly.tomsg(mp_witherror)                   # recovered <- Encode_1( Compress_Q(mp with error, 1) )
        for i in range(len(recovered_temp)):
            recovered[i] = recovered_temp[i]


########################################################################################################################
# Kyber512 Test code
# poly.tomsg,       poly.frommsg        : used in 'm'
# polyvec.tobytes,  polyvec.frombytes   : used in 'sk', 'pk'
# poly.compress,    poly.decompress     : used in 'ct'
# polyvec.compress, polyvec.decompress  : used in 'ct'
########################################################################################################################

# 1. Key generation part 
pk = [0]*PUBLICKEYBYTES
sk = [0]*SECRETKEYBYTES
kyber512.keygen(pk, sk)

# Gen message
if MESSAGETYPE == TESTVECTOR:
    m = [0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f]
    # m = [0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]
elif MESSAGETYPE == STDIN:
    buf = input(' * Input your message > ')
    m = [ord(item) for item in buf]
    for i in range(MSGBYTES-len(m)):
        m.append(0)
else:
    print(' Message type error!')
    exit()
    
# Gen random coins
coin = randombytes.get_randomarray()


# 2. Encryption part
ct = [0]*(POLYVECCOMPRESSEDBYTES + POLYCOMPRESSEDBYTES)
kyber512.enc(pk, m, coin, ct)


# 3. Decryption part
recovered = [0]*MSGBYTES
kyber512.dec(sk, ct, recovered)


if MESSAGETYPE == TESTVECTOR:
    printf.printarray('m = ', m)
    printf.printarray('recovered = ', recovered)
    pass
elif MESSAGETYPE == STDIN:
    recoveredstr = ""
    for item in recovered:
        if item == 0:
            break
        recoveredstr += chr(item)
        
    print(recoveredstr)
    
else:
    print(' Message type error!')
    exit()


########################################################################################################################
# Check random
########################################################################################################################

def check_cbd1(vec):
    res = [0]*(2*ETA1+1)
    for i in range(N):
        res[vec.coeffs[i]+ETA1] += 1
    print(res)
    
def check_cbd2(vec):
    res = [0]*(2*ETA2+1)
    for i in range(N):
        res[vec.coeffs[i]+ETA2] += 1
    print(res)

def check_uniform(A): # Number of A_entry is K*K*N(= 2*2*256), A_entry in R_Q
    res = [0]*Q
    for i in range(K):
        for j in range(K):
            for k in range(N):
                res[A.mat[i].vec[j].coeffs[k]] += 1
    for i in range(Q):
        print(res[i], end=' ')
        if i % 64 == 63:
            print()
    print()
    print('( sum( count of A_entry in Z_Q ) = ', sum(res), ' )')
    

""" Key generation """
# printf.printarray('publicseed = ', publicseed)
# printf.printarray('noiseseed = ', noiseseed)
# printf.printmatrix('A = ', A)
# printf.printvector('s = ', s)
# printf.printvector('e = ', e)
# printf.printvector('b = ', b)

# printf.printarray('m = ', m)
# printf.printarray('coin = ', coin)

""" Encryption """
# printf.printmatrix('Tr(A) = ', AT)
# printf.printvector('r = ', r)
# printf.printvector('e1 = ', e1)
# printf.printpoly('e2 = ', e2)
# printf.printvector('ct_u = ', ct_u)
# printf.printpoly('mp = ', mp)
# printf.printpoly('ct_v = ', ct_v)

""" Decryption """
# printf.printpoly('mp(with error) = ', mp_witherror)

""" Check function """
# check_uniform(A)
# check_cbd1(e.vec[0])
# check_cbd1(r.vec[0])
# check_cbd2(e1.vec[1])
# check_cbd2(e2)


# EOF