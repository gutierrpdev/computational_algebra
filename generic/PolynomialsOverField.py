from generic.EuclideanDomain import EuclideanDomain
from generic.Polynomial import Polynomial


# Basic implementation of a Polynomial RING over a Field.
# Note that the Ring K[X], where K is a Field, is in fact a Euclidean Domain and thus the Euclidean Algorithms
# may be used within our construction.
class PolynomialsOverField(EuclideanDomain):

    def div_mod(self, elem1, elem2):
        if elem1.degree() < elem2.degree():
            return Polynomial([self.base_field.add_id()], self.base_field), elem1

        a = list(elem1.coefs)
        a.reverse()
        b = list(elem2.coefs)
        b.reverse()
        quotient, remainder = self.extended_synthetic_division(a, b)
        quotient.reverse()
        remainder.reverse()
        return Polynomial(quotient, self.base_field), Polynomial(remainder, self.base_field)

    def __init__(self, base_field):
        self.base_field = base_field

    def add_id(self):
        return Polynomial([self.base_field.add_id()], self.base_field)

    def mul_id(self):
        return Polynomial([self.base_field.mul_id()], self.base_field)

    def add_inv(self, elem):
        lis = list(map(lambda x: self.base_field.add_inv(x), elem.coefs))
        return Polynomial(lis, self.base_field)

    def add(self, elem1, elem2):
        lis = list(map(lambda x1, x2: self.base_field.add(x1, x2), elem1.coefs, elem2.coefs))
        lis.extend(elem1.coefs[len(lis):])
        lis.extend(elem2.coefs[len(lis):])
        return Polynomial(lis, self.base_field)

    def mul(self, elem1, elem2):
        res = [self.base_field.add_id()] * (elem1.degree() + elem2.degree() + 1)
        for in1, e1 in enumerate(elem1.coefs):
            for in2, e2 in enumerate(elem2.coefs):
                res[in1 + in2] = self.base_field.add(res[in1 + in2], self.base_field.mul(e1, e2))
        return Polynomial(res, self.base_field)

    # Fast polynomial division by using Extended Synthetic Division.
    def extended_synthetic_division(self, dividend, divisor):

        out = list(dividend)  # Copy the dividend
        normalizer = divisor[0]
        for i in range(len(dividend) - len(divisor) + 1):
            out[i] = self.base_field.div(out[i], normalizer)  # for general polynomial division (when polynomials are non-monic),
            # we need to normalize by dividing the coefficient with the divisor's first coefficient
            coef = out[i]
            if coef != self.base_field.add_id():  # useless to multiply if coef is 0
                for j in range(1, len(divisor)):  # in synthetic division, we always skip the first coefficient of the divisor,
                    # because it is only used to normalize the dividend coefficients
                    aux = self.base_field.add_inv(self.base_field.mul(divisor[j], coef))
                    out[i + j] = self.base_field.add(out[i + j], aux)
        separator = len(dividend) - len(divisor) + 1
        quotient = out[:separator]
        remainder = out[separator:]

        if len(quotient) == 0:
            quotient = [self.base_field.add_id()]
        if len(remainder) == 0:
            remainder = [self.base_field.add_id()]
        return quotient,  remainder  # return quotient, remainder.
