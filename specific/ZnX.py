from generic.EuclideanDomain import EuclideanDomain
from generic.Ring import Ring
from generic.Polynomial import Polynomial
from specific.Zn import Zn


class ZnX(EuclideanDomain):

    def div_mod(self, elem1, elem2):
        pass

    def __init__(self, n):
        self.n = n

    def add_id(self):
        return Polynomial([0], Zn(self.n))

    def mul_id(self):
        return Polynomial([1], Zn(self.n))

    def add_inv(self, elem):
        lis = list(map(lambda x: Zn(self.n).add_inv(x), elem.coefs))
        return Polynomial(lis, Zn(self.n))

    def add(self, elem1, elem2):
        lis = list(map(lambda x1, x2: Zn(self.n).add(x1, x2), elem1.coefs, elem2.coefs))
        lis.extend(elem1.coefs[len(lis):])
        lis.extend(elem2.coefs[len(lis):])
        return Polynomial(lis, Zn(self.n))

    def mul(self, elem1, elem2):
        res = [0] * (elem1.degree() + elem2.degree() + 1)
        for in1, e1 in enumerate(elem1.coefs):
            for in2, e2 in enumerate(elem2.coefs):
                res[in1 + in2] += Zn(self.n).mul(e1, e2)
        return Polynomial(res, Zn(self.n))

    # Fast polynomial division by using Extended Synthetic Division.
    def extended_synthetic_division(self, dividend, divisor):
        dividend.reverse()
        divisor.reverse()
        out = list(dividend)  # Copy the dividend
        normalizer = divisor[0]
        for i in range(len(dividend) - len(divisor) + 1):
            out[i] = Zn(self.n).div(out[i], normalizer)
            # for general polynomial division (when polynomials are non-monic),
            # we need to normalize by dividing the coefficient with the divisor's first coefficient
            coef = out[i]
            if coef != Zn(self.n).add_id():  # useless to multiply if coef is 0
                for j in range(1, len(divisor)):
                    # in synthetic division, we always skip the first coefficient of the divisor,
                    # because it is only used to normalize the dividend coefficients
                    aux = Zn(self.n).add_inv(Zn(self.n).mul(divisor[j], coef))
                    out[i + j] = Zn(self.n).add(out[i + j], aux)
        separator = len(dividend) - len(divisor) + 1
        quotient = out[:separator]
        quotient.reverse()
        remainder = out[separator:]
        remainder.reverse()
        dividend.reverse()
        divisor.reverse()
        if len(quotient) == 0:
            quotient = [Zn(self.n).add_id()]
        if len(remainder) == 0:
            remainder = [Zn(self.n).add_id()]
        return quotient, remainder  # return quotient, remainder.

