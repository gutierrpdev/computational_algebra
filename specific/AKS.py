from math import sqrt, floor

from generic.Polynomial import Polynomial
from specific.ZX import ZX
from specific.ZZ import ZZ


class AKS:

    @staticmethod
    def compute_primality(n):
        if n <= 1:
            return False

        if AKS.is_power(n):
            return False

        r = AKS.compute_r(n)

        if r == n:
            return True

        if ZZ().gcd(n, r) > 1:
            return False

        top = 2 * (len(bin(n)) - 2) * floor(sqrt(r)) + 1
        for j in range(1, top + 1):
            _ZX = ZX()
            # (X + j)^n
            elem1 = _ZX.exp(Polynomial([j, 1], ZZ()), n)
            # X^n + j
            aux = [0] * (n + 1)
            aux[n] = 1
            aux[0] = j
            elem2 = Polynomial(aux, ZZ())

            # (X^r - 1)
            aux2 = [0] * (r + 1)
            aux2[r] = 1
            aux2[0] = -1
            divisor = Polynomial(aux2, ZZ())

            # mod (X^r - 1)
            res1 = Polynomial(_ZX.extended_synthetic_division(elem1.coefs, divisor.coefs)[1], ZZ())
            res2 = Polynomial(_ZX.extended_synthetic_division(elem2.coefs, divisor.coefs)[1], ZZ())

            # convert to Zn[X]
            res1 = Polynomial([elem % n for elem in res1.coefs], ZZ())
            res2 = Polynomial([elem % n for elem in res2.coefs], ZZ())

            if res1 != res2:
                return False
        return True

    @staticmethod
    def is_power(n):
        # 2 ^ k
        if (- n & n) == n:
            return True
        lgn = 1 + (len(bin(abs(n))) - 2)
        for b in range(2, lgn):
            # b lg a = lg n
            low_a = 1
            high_a = 1 << (lgn // b + 1)
            while low_a < high_a - 1:
                mid_a = (low_a + high_a) >> 1
                ab = (mid_a ** b)
                if ab > n:
                    high_a = mid_a
                elif ab < n:
                    low_a = mid_a
                else:
                    # mid_a ^ b
                    return True
        return False

    @staticmethod
    def compute_r(n):
        r = 2
        while r <= n:
            if ZZ().gcd(n, r) != 1:
                return r
            else:
                mul_order = 1
                nc = n
                # Multiply by n modulo r until nc becomes 1
                while nc != 1:
                    nc *= n
                    nc %= r
                    mul_order += 1
                if mul_order > 4 * (len(bin(n)) - 2) ** 2:
                    return r
            r += 1


if __name__ == "__main__":
    print(AKS.is_power(8))
    print(AKS.is_power(27))
    print(AKS.is_power(97))
    print(AKS.is_power(3021377))
    print(AKS.is_power(3021373))
    print(AKS.is_power(6))
    print(AKS.is_power(81))

    print("...............")
    print("check if previous numbers are prime")
    print(AKS.compute_primality(8))
    print(AKS.compute_primality(27))
    print(AKS.compute_primality(97))
    print(AKS.compute_primality(379))
    print(AKS.compute_primality(381))
    print(AKS.compute_primality(6))
    print(AKS.compute_primality(81))

