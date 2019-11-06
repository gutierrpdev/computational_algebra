from generic.Field import Field
from generic.Polynomial import Polynomial
from generic.Polynomials import Polynomials
from specific.Zp import Zp


class FiniteField(Field):

    def __init__(self, p, k, f):
        self.p, self.k, self.f = p, k, f
        self.field = Polynomials(Zp(p))

    def mul_inv(self, elem):
        _, d, _ = self.extended_gcd(elem, self.f)
        return d

    def add_id(self):
        return Polynomial([0], self.field)

    def mul_id(self):
        return Polynomial([1], self.field)

    def add_inv(self, elem):
        pass

    def add(self, elem1, elem2):
        pass

    def mul(self, elem1, elem2):
        pass
