import random
import hashlib

"""
    For keccak_absorb
"""
def load64_littleendian(state0, state1, state2, state3, state4, state5, state6, state7):
    return (state7 << 56) | (state6 << 48) | (state5 << 40) | (state4 << 32) | (state3 << 24) | (state2<<16) | (state1<<8) | (state0)

def keccak_absorb(r, m, mlen, p):
    t = [0]*200
    
    s = [0]*25
    inx = 0
    while mlen >= r:
        for i in range(r//8):
            s[i] ^= load64_littleendian(m[inx + 8*i], m[inx + 8*i+1], inx + m[8*i+2], inx + m[8*i+3], inx + m[8*i+4], inx + m[8*i+5], inx + m[8*i+6], inx + m[8*i+7])
    
        # keccakf1600_statepermute(s)
        mlen -= r
        inx += r
        
    for i in range(mlen):
        t[i] = m[i]
    t[mlen] = p
    t[r-1] |= 128
    for i in range(r//8):
        s[i] ^= load64_littleendian(t[8*i], t[8*i+1], t[8*i+2], t[8*i+3], t[8*i+4], t[8*i+5], t[8*i+6], t[8*i+7])
        
    return s # s = []*25


########################################################################################################################
"""
    Kyber parameters
"""
global Q, N, K, ETA1, ETA2
Q           = 3329
N           = 256
K           = 2
ETA1        = 3
ETA2        = 2

global ARRAYSIZE
ARRAYSIZE   = 32 # fixed

"""
    Debug options
"""
global DEBUGTYPE, HEX, DECIMAL
HEX         = 16
DECIMAL     = 10
DEBUGTYPE   = HEX

""" 
    Define type of message
        - If you want to use an already defined test vector as a message, define this value(MESSAGETYPE) as TESTVECTOR, 
          and if you want a message as standard input, define this value(MESSAGETYPE) as STDIN
"""
global MESSAGETYPE, TESTVECTOR, STDIN
STDIN       = 1
TESTVECTOR  = 2
MESSAGETYPE = STDIN

"""
    For cbd1
    (24-bit integer <- 3-byte array)
"""
def load24_littleendian(state0, state1, state2):
    return (state2<<16) | (state1<<8) | (state0)

"""
    For cbd2
    (32-bit integer <- 4-byte array)
"""
def load32_littleendian(state0, state1, state2, state3):
    return (state3<<24) | (state2<<16) | (state1<<8) | (state0)
    

"""
    For sample vector
    (res in B^192 <- prf(seed in B^32, nonce) using shake256)
    (192byte <- seed : 32byte)
"""
def prf(outlen, seed, nonce):
    
    strseed = "".join(format(byte, '08b') for byte in seed).encode() # str <- bytes
    strseed += format(nonce, '08b').encode() # concate nonce 
    
    state = hashlib.new('shake_256')
    state.update(strseed)
    
    state = state.hexdigest(outlen)
    state = bytearray.fromhex(state)
    
    # res in B^128 (or B^192) <- prf(seed in B^32, nonce)
    res = [0]*outlen
    for i in range(outlen):
        res[i] = state[i]

    return res
    

""" 
    For sample a polyvec 's', 'e', 'r'
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
    For sample a polyvec 's', 'e', 'r'
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
    For sample a matrix 'A'
    (Kyber function : rej_uniform)
    
    Parameters :
                - coeffs : output which is element of R_Q (entry of matrix A)
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
    
    state = hashlib.new('shake_128')
    state.update(extseed)
    
    state = state.hexdigest(506) # fixed xof_state structure : len = 25 
    state = bytearray.fromhex(state)
    
    res = [0]*506
    for i in range(506):
        res[i] = state[i]
        
    return res

"""
    Print matrix, vector and element
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
            for i in range(ARRAYSIZE):
                print('%02x' %state[i], end=' ')
        elif DEBUGTYPE == DECIMAL:
            for i in range(ARRAYSIZE):
                print('%03d' %state[i], end=' ')
                
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
    Get random byte-array
    (publicseed, noiseseed, coin) 
    (Used to sample  'A',   's, e',   'e1, e2, r')
"""
class randombytes:

    def get_randomarray():
        return [random.randint(0, 256) for _ in range(ARRAYSIZE)]

"""
    Define poly (== R_Q) and operations
"""
class poly:
    def __init__(self):
        self.Q = Q
        self.N = N
        self.coeffs = [0 for _ in range(N)]

    def __add__(self, other):
        res = poly()
        for i in range(self.N):
            res.coeffs[i] = (self.coeffs[i] + other.coeffs[i]) % self.Q
        return res
    
    def __sub__(self, other):
        res = poly()
        for i in range(self.N):
            if self.coeffs[i] >= other.coeffs[i]:
                res.coeffs[i] = self.coeffs[i] - other.coeffs[i]
            else:
                res.coeffs[i] = self.coeffs[i] + self.Q - other.coeffs[i]
        return res

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
    
    @staticmethod
    def sample(seed, nonce, eta):
        if eta == 2:
            return cbd2(prf(eta*N//4, seed, nonce)) # outlen = 128, nonce = 0, 1
        elif eta == 3:
            return cbd1(prf(eta*N//4, seed, nonce)) # outlen = 128, nonce = 0, 1
        
    
    @staticmethod
    def frommsg(m):
        res = poly()
        
        for i in range(N//8):
            for j in range(8):
                if (m[i]>>j)&1 == 1:
                    res.coeffs[8*i+j] = 1665
        return res

    @staticmethod
    def tomsg(mp):
        res = [0]*ARRAYSIZE
        
        for i in range(N//8):
            for j in range(8):
                t = (((mp.coeffs[8*i+j]<<1) + Q//2) //Q) &1
                res[i] |= t << j
        return res
    
"""
    Define polyvec (== R_Q^K) and operations
"""
class polyvec:
    def __init__(self):
        self.Q = Q
        self.N = N
        self.K = K
        self.vec = [poly() for _ in range(K)]
    
    def __add__(self, other):
        res = polyvec()
        for i in range(self.K):
            res.vec[i] = self.vec[i] + other.vec[i]
        return res
    
    def __sub__(self, other):
        res = polyvec()
        for i in range(self.K):
            res.vec[i] = self.vec[i] - other.vec[i]
        return res

    def __mul__(self, other):
        res = poly()
        temp = poly()
        for i in range(self.K):
            temp = self.vec[i] * other.vec[i]
            res = res + temp
        return res
    
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
    Define polymat (== R_Q^(K*K)) and operations
"""
class polymat:
    def __init__(self):
        self.Q = Q
        self.N = N
        self.K = K
        self.mat = [polyvec() for _ in range(K)]

    def __mul__(self, other): # matrix * vector
        res = polyvec()
        for i in range(K):
            res.vec[i] = self.mat[i] * other
        return res
    
    @staticmethod
    def trans(A):
        res = polymat()
        for i in range(K):
            for j in range(K):
                res.mat[i].vec[j].coeffs[:] = A.mat[j].vec[i].coeffs[:]
        return res

    def transmatrix(self):
        res = polymat()
        for i in range(K):
            for j in range(K):
                res.mat[i].vec[j].coeffs[:] = self.mat[j].vec[i].coeffs[:]
        return res
    
    def sample(self, seed, transpose):
            
        for i in range(K):
            for j in range(K):
                if transpose == True:
                    buf = xof(seed, j, i)
                else:
                    buf = xof(seed, i, j)
                ctr = parse(self.mat[i].vec[j].coeffs, self.N, buf, len(buf))
                
                if ctr != self.N:
                    self.sample(self, seed)
            
        
########################################################################################################################
# 1. Key generate part
publicseed = randombytes.get_randomarray()
noiseseed = randombytes.get_randomarray()

A, s, e = polymat(), polyvec(), polyvec()
nonce = 0

A.sample(publicseed, False)
s.sample(noiseseed, nonce, ETA1)
nonce+=K
e.sample(noiseseed, nonce, ETA1)

b = A * s + e

# Gen message
if MESSAGETYPE == TESTVECTOR:
    m = [0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f]
elif MESSAGETYPE == STDIN:
    buf = input(' * Input your message > ')
    m = [ord(item) for item in buf]
    for i in range(32-len(m)):
        m.append(0)
else:
    print(' Message type error!')
    exit()

# Gen coin
coin = randombytes.get_randomarray()


# 2. Encryption part
AT , r, e1, e2 = polymat(), polyvec(), polyvec(), poly()
nonce = 0

AT.sample(publicseed, True)
r.sample(coin, nonce, ETA1); nonce += K
e1.sample(coin, nonce, ETA2); nonce += K
e2 = poly.sample(coin, nonce, ETA2)

mp = poly.frommsg(m)
ct_u = AT * r + e1
ct_v = b * r + e2 + mp


# 3. Decryption part
mp_witherror = ct_v - (s * ct_u)
recovered = poly.tomsg(mp_witherror)


########################################################################################################################



printf.printarray('publicseed = ', publicseed)
printf.printarray('noiseseed = ', noiseseed)
printf.printmatrix('A = ', A)
printf.printvector('s = ', s)
printf.printvector('e = ', e)
printf.printvector('b = ', b)

printf.printarray('m = ', m)
printf.printarray('coin = ', coin)

printf.printmatrix('Tr(A) = ', AT)
printf.printvector('r = ', r)
printf.printvector('e1 = ', e1)
printf.printpoly('e2 = ', e2)
printf.printvector('ct_u = ', ct_u)
printf.printpoly('mp = ', mp)
printf.printpoly('ct_v = ', ct_v)

printf.printpoly('mp(with error) = ', mp_witherror)

if MESSAGETYPE == TESTVECTOR:
    printf.printarray('recovered = ', recovered)
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
    

# check_uniform(A)
# check_cbd1(e.vec[0])
# check_cbd1(r.vec[0])
# check_cbd2(e1.vec[1])
# check_cbd2(e2)


# EOF