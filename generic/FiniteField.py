from generic.Field import Field
from generic.Polynomial import Polynomial
from generic.PolynomialsOverField import PolynomialsOverField
from specific.Zp import Zp


class FiniteField(Field):

    def __init__(self, p, k, f):
        self.p, self.k, self.f = p, k, f
        self.field = PolynomialsOverField(Zp(p))

    def cardinality(self):
        return self.p ** self.k

    def div_mod(self, elem1, elem2):
        # make both polynomials of degree lower than f's
        elem1 = self.rep(elem1)
        elem2 = self.rep(elem2)

        # make sure that elem1 has a higher degree than elem2's.
        # note that adding f is essentially no different from adding 0 in this field.
        elem1 = self.field.add(self.f, elem1)

        return self.field.div_mod(elem1, elem2)

    def add_id(self):
        return Polynomial([0], self.field)

    def mul_id(self):
        return Polynomial([1], self.field)

    def add_inv(self, elem):
        return self.rep(self.field.add_inv(elem))

    def add(self, elem1, elem2):
        return self.rep(self.field.add(elem1, elem2))

    def mul(self, elem1, elem2):
        return self.rep(self.field.mul(elem1, elem2))

    # elem's equivalence class
    def rep(self, elem):
        return self.field.div_mod(elem, self.f)[1]


if __name__ == "__main__":
    F9 = FiniteField(3, 2, Polynomial([2, 1, 1], Zp(3)))
    aux = F9.div(F9.mul_id(), Polynomial([0, 1], Zp(3)))
    aux2 = F9.mul_inv(Polynomial([0, 1], Zp(3)))
    assert(aux == aux2)
    print(F9.mul(aux, Polynomial([0, 1], Zp(3))).coefs)
