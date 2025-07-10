
import sympy
import random
from sympy import nextprime
from fpylll import IntegerMatrix, LLL
import math

# Steg 1: skapa två 64-bitars primtal p och q
bit_size = 64
suffix_bits = 16  # Simulerad DRMM-nivå (antal kända bitar)

# Slumpmässiga primtal p och q
def generate_prime_with_suffix(bits, suffix):
    while True:
        high = random.getrandbits(bits - suffix.bit_length())
        candidate = (high << suffix.bit_length()) + suffix
        if sympy.isprime(candidate):
            return candidate

# Generera primtal p, q med samma suffix (t.ex. 0b1010110101011101)
suffix = random.getrandbits(suffix_bits)
p = generate_prime_with_suffix(bit_size, suffix)
q = generate_prime_with_suffix(bit_size, suffix)
N = p * q

# Steg 2: definiera kända bitar
b = p & ((1 << suffix_bits) - 1)
d = q & ((1 << suffix_bits) - 1)

# Steg 3: beräkna det justerade värdet
K = N - b * d  # K = ac + ad + bc

# Steg 4: bygg polynomet f(a, c) = ac + ad + bc - K
# Vi approximerar ad + bc ≈ (a + c) * suffix → gör en 1D-approximation:
# f(x) = x^2 + 2*b*x - K ≡ 0

# Vi försöker hitta små rötter till f(x) = x^2 + 2*b*x - K (modulo N)
# Genom lattice-teknik

# Bygg lattice-basen
X = 1 << (bit_size - suffix_bits)  # Gissning för hur stor x är (~ökända bitar)
f0 = [1, 2 * b, -K]

# Skapa en 3x3 lattice för x^2 + 2bx - K ≡ 0
B = IntegerMatrix(3, 3)
B[0, 0] = X * X
B[0, 1] = 0
B[0, 2] = 0
B[1, 0] = X * f0[1]
B[1, 1] = X
B[1, 2] = 0
B[2, 0] = f0[2]
B[2, 1] = 0
B[2, 2] = 1

# Kör LLL-reduktion
LLL.reduction(B)

# Extrahera koefficienterna från första vektorn
a, b_, c = B[0, 0], B[0, 1], B[0, 2]

# Lös den approximativa polynomekvationen: a*x^2 + b_*x + c = 0
roots = sympy.solve(a * sympy.Symbol('x')**2 + b_ * sympy.Symbol('x') + c, sympy.Symbol('x'))

# Filtrera heltalsrötter
int_roots = [int(r) for r in roots if r.is_real and abs(r - int(r)) < 1e-6]

# Presentera
print("Suffix (simulerat från DRMM):", bin(suffix))
print("p =", p)
print("q =", q)
print("N =", N)
print("Approximerade rötter (troliga övre bitar):", int_roots)
print("Verkliga övre bitar p:", p >> suffix_bits)
print("Verkliga övre bitar q:", q >> suffix_bits)
