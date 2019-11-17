from generic.Ring import Ring
from generic.Polynomial import Polynomial
from specific.ZZ import ZZ


class ZX(Ring):
    def add_id(self):
        return Polynomial([0], ZZ())

    def mul_id(self):
        return Polynomial([1], ZZ())

    def add_inv(self, elem):
        lis = list(map(lambda x: ZZ().add_inv(x), elem.coefs))
        return Polynomial(lis, ZZ())

    def add(self, elem1, elem2):
        lis = list(map(lambda x1, x2: ZZ().add(x1, x2), elem1.coefs, elem2.coefs))
        lis.extend(elem1.coefs[len(lis):])
        lis.extend(elem2.coefs[len(lis):])
        return Polynomial(lis, ZZ())

    def mul(self, elem1, elem2):
        res = [0] * (elem1.degree() + elem2.degree() + 1)
        for in1, e1 in enumerate(elem1.coefs):
            for in2, e2 in enumerate(elem2.coefs):
                res[in1 + in2] += ZZ().mul(e1, e2)
        return Polynomial(res, ZZ())

    # Here we make use of the following theorem:
    # mcd(a, b) = mcd(b, primitive_part(r)), if a, b in Z[X] and r is such that a = bq + r, with deg(r) < deg(b)
    # this can be done under the assumption that the coefficients of a are "convenient", that is,
    # at each step of the synthetic division we are able to perform a division valid in Z[X].
    # Note that, in order to force this situation, a little trick may be performed on a:
    # it suffices to multiply the whole polynomial a by the leading coefficient of b to guarantee that
    # the previously described condition is always fulfilled.
    # NOTE: Z[X] is NOT an euclidean domain (think of a = x + 1 and b = 2x + 1)
    def gcd(self, a, b):
        c = self.primitive_part(a)
        d = self.primitive_part(b)

        if c.degree() < d.degree():
            c, d = d, c

        while d != self.add_id():
            lc = c.get_leading_coef()
            n = c.degree() - d.degree() + 1
            aux = Polynomial([lc ** n], ZZ())
            aux = self.mul(aux, c)
            r = Polynomial(self.extended_synthetic_division(aux.coefs, d.coefs)[1], ZZ())
            c = d
            d = self.primitive_part(r)
        gcd_ab = ZZ().gcd(self.mcd_polynomial_coefs(a), self.mcd_polynomial_coefs(b))
        return self.mul(Polynomial([gcd_ab], ZZ()), c)

    # polynomial divided by gcd of all its coefficients
    def primitive_part(self, elem):
        if elem == self.add_id():
            return self.add_id()
        else:
            # divide all coefficients by common gcd to obtain primitive part
            return Polynomial(list(map(lambda x: ZZ().div(x, self.mcd_polynomial_coefs(elem)), elem.coefs)), ZZ())

    @staticmethod
    def mcd_polynomial_coefs(elem):
        # gcd not defined for zero values: remove them all before proceeding
        without_zeros = list(filter(lambda x: x != 0, elem.coefs))
        # calculate gcd of all non-zero values in list
        return ZZ().multiple_gcd(without_zeros)

    # Fast polynomial division by using Extended Synthetic Division.
    @staticmethod
    def extended_synthetic_division(dividend, divisor):
        dividend.reverse()
        divisor.reverse()
        out = list(dividend)  # Copy the dividend
        normalizer = divisor[0]
        for i in range(len(dividend) - len(divisor) + 1):
            out[i] = ZZ().div(out[i], normalizer)
            # for general polynomial division (when polynomials are non-monic),
            # we need to normalize by dividing the coefficient with the divisor's first coefficient
            coef = out[i]
            if coef != ZZ().add_id():  # useless to multiply if coef is 0
                for j in range(1, len(divisor)):
                    # in synthetic division, we always skip the first coefficient of the divisor,
                    # because it is only used to normalize the dividend coefficients
                    aux = ZZ().add_inv(ZZ().mul(divisor[j], coef))
                    out[i + j] = ZZ().add(out[i + j], aux)
        separator = len(dividend) - len(divisor) + 1
        quotient = out[:separator]
        quotient.reverse()
        remainder = out[separator:]
        remainder.reverse()
        dividend.reverse()
        divisor.reverse()
        if len(quotient) == 0:
            quotient = [ZZ().add_id()]
        if len(remainder) == 0:
            remainder = [ZZ().add_id()]
        return quotient, remainder  # return quotient, remainder.


if __name__ == "__main__":
    ZX = ZX()
    _a = Polynomial([1, 2, 1], ZZ())
    _b = Polynomial([1, 0, 2], ZZ())
    _c = Polynomial([1, 2, 2, 1], ZZ())
    _d = Polynomial([27, 9, 18, 0, 9], ZZ())
    _e = Polynomial([3, 6, 3], ZZ())
    _f = Polynomial([3, 6, 6, 3], ZZ())
    _g = Polynomial([1, 0, 1], ZZ())
    _h = Polynomial([5, 3], ZZ())
    _i = Polynomial([1, 1], ZZ())
    _j = Polynomial([1, 2], ZZ())
    # print(ZX.extended_synthetic_division(_c.coefs, _a.coefs))
    # print(ZX.mcd_polynomial_coefs(_d))
    # print(ZX.primitive_part(_d).coefs)
    print(ZX.gcd(_c, _a).coefs)
    print(ZX.gcd(_a, _c).coefs)
    print(ZX.gcd(_d, _b).coefs)
    print(ZX.gcd(_e, _f).coefs)
    print(ZX.gcd(_g, _h).coefs)
    print(ZX.gcd(_i, _j).coefs)
