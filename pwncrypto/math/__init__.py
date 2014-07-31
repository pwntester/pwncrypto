'''Mathematical building blocks used for basic cryptographic theory and
operations.'''

from foreign import *
from gmpy2 import mpz
import random

def check_prime(p):
    """Miller-Rabin prime test"""
    if p & 1 == 0:
        return False

    m = p - 1
    s = 0
    while m & 1 == 0:
        m >>= 1
        s += 1

    for j in range(100):
        a = random.randint(2, p - 2)
        if gcd(a, p) != 1:
            return False

        b = pow(a, m * (1 << s), p)
        if b in (0, 1, p - 1):
            continue

        for i in range(s):
            b = pow(b, 2, p)

            if b == 1:
                return False

            if b == p - 1:
                if i < s - 1:
                    break
                else:
                    return False
        else:
            return False
    return True

def generate_prime(keysize=1024):
    """ Return a random prime number of keysize bits in size. """
    while True:
        num = random.randrange(2 ** (keysize - 1), 2 ** (keysize))
        if check_prime(num):
            return num

def crt(ml, al):
    """
    Chinese Remainder Theorem:
    ms = list of pairwise relatively prime integers
    as = remainders when x is divided by ms
    (ai is 'each in as', mi 'each in ms')

    The solution for x modulo M (M = product of ms) will be:
    x = a1*M1*y1 + a2*M2*y2 + ... + ar*Mr*yr (mod M),
    where Mi = M/mi and yi = (Mi)^-1 (mod mi) for 1 <= i <= r.
    """

    M = reduce(lambda x, y: x * y, ml)        # multiply ml together
    Ms = [M / mi for mi in ml]   # list of all M/mi

    ys = [invert(Mi, mi) for Mi, mi in zip(Ms, ml)]  # uses inverse,eea
    return reduce(lambda x, y: x + y, [mpz(ai) * mpz(Mi) * yi for ai, Mi, yi in zip(al, Ms, ys)]) % M
