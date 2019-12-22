"""
Pollard Rho Algorithm
-> for computing the discrete logarithm of an element in the multiplicative group of a finite field.
:param g: generator of the multiplicative cyclic group of Fq (in Polynomial form)
:param h: element in Fq whose discrete logarithm base g we are interested in finding (in Polynomial form)
:param finite_field: finite field Fq
"""
from generic.FiniteField import FiniteField
from generic.Polynomial import Polynomial
from specific.ZZ import ZZ
from specific.Zp import Zp


class DiscreteLogarithm:

    def __init__(self, g, h, finite_field):
        self.g, self.h, self.finite_field = g, h, finite_field
        # p is the order of the multiplicative group of Fq.
        self.p = self.finite_field.cardinality() - 1
        self.q = self.p // 2

    # given elem in polynomial form, we define subset(elem) as the result of converting the polynomial coefficients of
    # elem as a number base p to an integer base 10 mod 3 (we require 3 different subsets of approximately equal size
    # in order for the algorithm to work as expected).
    def subset(self, elem):
        res = 0
        acc = 1
        # compute number representation of polynomial as an integer (base 10).
        for x in elem.coefs:
            res += x * acc
            acc *= self.finite_field.p
        return res % 3

    def _f(self, elem):
        subset = self.subset(elem)
        if subset == 0:
            return self.finite_field.mul(self.h, elem)
        elif subset == 1:
            return self.finite_field.mul(elem, elem)
        else:
            return self.finite_field.mul(self.g, elem)

    def _h(self, elem, n):
        subset = self.subset(elem)
        if subset == 0:
            return n
        elif subset == 1:
            return (2 * n) % self.q
        else:
            return (n + 1) % self.q

    def _g(self, elem, n):
        subset = self.subset(elem)
        if subset == 0:
            return (n + 1) % self.q
        elif subset == 1:
            return (2 * n) % self.q
        else:
            return n

    def pollard(self):
        ai, bi, xi = 1, 1, self.finite_field.mul(self.g, self.h)
        a2i, b2i, x2i = 1, 1, self.finite_field.mul(self.g, self.h)
        print("subset x:", self.subset(xi))

        while True:
            xi, ai, bi = self._f(xi), self._g(xi, ai), self._h(xi, bi)
            x2i, a2i, b2i = self._f(self._f(x2i)), self._g(self._f(x2i), a2i), self._h(self._f(x2i), b2i)

            print(xi, ai, bi)
            print(x2i, a2i, b2i)

            if xi == x2i:
                r = bi - b2i
                if r == 0:
                    return None
                x = ZZ().extended_gcd(r, self.p)[1] * (ai - a2i) % self.q
                return x


if __name__ == "__main__":
    F4 = FiniteField(2, 2, Polynomial([1, 1, 1], Zp(2)))
    F4Alpha = Polynomial([0, 1], Zp(2))
    F4AlphaPlusOne = Polynomial([1, 1], Zp(2))
    F4One = Polynomial([1], Zp(2))
    _x = DiscreteLogarithm(F4Alpha, F4AlphaPlusOne, F4)
    print(_x.pollard())
