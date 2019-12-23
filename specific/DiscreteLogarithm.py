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
    # in order for the algorithm to work as expected). Need mul_id not to be in set 1, therefore, we need to add
    # one to the result.
    def subset(self, elem):
        res = 0
        acc = 1
        # compute number representation of polynomial as an integer (base 10).
        for x in elem.coefs:
            res += x * acc
            acc *= self.finite_field.p
        return (res + 1) % 3

    # f is defined such a way that given an element, it can multiply by either g, h or itself. _g and _h are
    # defined to return how many times we have multiplied by g and h.
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
            return (2 * n) % self.p
        else:
            return (n + 1) % self.p

    def _g(self, elem, n):
        subset = self.subset(elem)
        if subset == 0:
            return (n + 1) % self.p
        elif subset == 1:
            return (2 * n) % self.p
        else:
            return n

    def pollard(self):
        ai, bi, xi = 0, 0, self.finite_field.mul_id()
        a2i, b2i, x2i = 0, 0, self.finite_field.mul_id()

        while True:
            xi, ai, bi = self._f(xi), self._h(xi, ai), self._g(xi, bi)
            x2i, a2i, b2i = self._f(self._f(x2i)), self._h(self._f(x2i), self._h(x2i, a2i)), \
                            self._g(self._f(x2i), self._g(x2i, b2i))

            if xi == x2i:
                r = bi - b2i
                if r == 0:
                    return None
                x = (ZZ().extended_gcd(r, self.p)[1] * (a2i - ai)) % self.p
                return x


if __name__ == "__main__":
    F4 = FiniteField(2, 2, Polynomial([1, 1, 1], Zp(2)))
    F4Alpha = Polynomial([0, 1], Zp(2))
    F4AlphaPlusOne = Polynomial([1, 1], Zp(2))
    F4One = Polynomial([1], Zp(2))
    _x = DiscreteLogarithm(F4AlphaPlusOne, F4One, F4)
    print(_x.pollard())
