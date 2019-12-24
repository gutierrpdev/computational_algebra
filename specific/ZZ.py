from random import random, randint

from generic.EuclideanDomain import EuclideanDomain


class ZZ(EuclideanDomain):
    def rem(self, elem1, elem2):
        return elem1 % elem2

    def add_id(self):
        return 0

    def mul_id(self):
        return 1

    def add_inv(self, elem):
        return -elem

    def add(self, elem1, elem2):
        return elem1 + elem2

    def div_mod(self, elem1, elem2):
        return elem1 // elem2, elem1 % elem2

    def mul(self, elem1, elem2):
        return elem1 * elem2

    def mul_scalar(self, elem, scalar):
        return elem*scalar

    @staticmethod
    def miller_rabin_prime_check(n, k=4):
        if n <= 3 or n % 2 == 0:
            return False

        r, d = ZZ.rd_form(n)

        for _ in range(k):
            if not ZZ.miller_test(r, d, n):
                return False
        return True

    @staticmethod
    def miller_test(r, d, n):
        a = randint(2, n - 2)
        x = pow(a, d, n)
        if x == 1 or x == (n - 1):
            return True
        for _ in range(r - 1):
            x = (x * x) % n
            if x == 1:
                return False
            elif x == (n - 1):
                return True
        return False

    '''
    :returns: input expressed as 2^r * d + 1
    '''
    @staticmethod
    def rd_form(n):
        r = 0
        d = n - 1
        while d % 2 == 0:
            d //= 2
            r += 1
        return r, d


if __name__ == "__main__":
    ZZ = ZZ()
    print("check if numbers from AKS test are 'probably' prime")
    print(ZZ.miller_rabin_prime_check(8))
    print(ZZ.miller_rabin_prime_check(27))
    print(ZZ.miller_rabin_prime_check(79))
    print(ZZ.miller_rabin_prime_check(83))
    print(ZZ.miller_rabin_prime_check(89))
    print(ZZ.miller_rabin_prime_check(97))
    print(ZZ.miller_rabin_prime_check(379))
    print(ZZ.miller_rabin_prime_check(381))
    print(ZZ.miller_rabin_prime_check(6))
    print(ZZ.miller_rabin_prime_check(81))
