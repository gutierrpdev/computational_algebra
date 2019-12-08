from generic.Field import Field
from generic.Polynomial import Polynomial
from generic.PolynomialsOverField import PolynomialsOverField
from specific.Zp import Zp
import copy


class FiniteField(Field):

    def __init__(self, p, k, f):
        self.p, self.k, self.f = p, k, f
        self.field = PolynomialsOverField(Zp(p))

    def cardinality(self):
        return self.p ** self.k

    def div_mod(self, elem1, elem2):

        e1 = copy.deepcopy(elem1)
        e2 = copy.deepcopy(elem2)

        # make both polynomials of degree lower than f's
        e1 = self.rep(e1)
        e2 = self.rep(e2)

        # make sure that elem1 has a higher degree than elem2's.
        # note that adding f is essentially no different from adding 0 in this field.
        # e1 = self.field.add(self.f, e1)

        # result = self.field.div_mod(e1, e2)
        result = self.mul(e1, self.mul_inv(e2))
        # make sure result is returned in proper representation.
        return self.rep(result), self.add_id()

    def add_id(self):
        return Polynomial([0], self.field)

    def mul_id(self):
        return Polynomial([1], self.field)

    def add_inv(self, elem):
        return self.rep(self.field.add_inv(elem))

    def mul_inv(self, elem):
        res = self.field.extended_gcd(elem, self.f)
        return res[1]

    def add(self, elem1, elem2):
        return self.rep(self.field.add(elem1, elem2))

    def mul(self, elem1, elem2):
        return self.rep(self.field.mul(elem1, elem2))

    # elem's equivalence class
    def rep(self, elem):
        return self.field.div_mod(elem, self.f)[1]

    def __eq__(self, other):
        if isinstance(other, FiniteField):
            return self.p == other.p and self.k == other.k
        return NotImplemented


if __name__ == "__main__":

    # test to check if rep works correctly for GF(4)
    Z2 = Zp(2)
    F4 = FiniteField(2, 2, Polynomial([1, 1, 1], Z2))
    F4One = F4.mul_id()
    F4Zero = F4.add_id()
    F4Alpha = Polynomial([0, 1], Z2)
    F4AlphaPlusOne = Polynomial([1, 1], Z2)
    # Both of the following elements should have rep [1]
    print(F4.rep(F4One), F4.rep(F4.add(F4.mul(F4Alpha, F4Alpha), F4Alpha)))
    # Multiplication table for F4.
    print("Multiplication table for F4")
    print(F4.mul(F4One, F4One))
    print(F4.mul(F4One, F4Alpha))
    print(F4.mul(F4One, F4AlphaPlusOne))
    print(F4.mul(F4Alpha, F4Alpha))
    print(F4.mul(F4Alpha, F4AlphaPlusOne))
    print(F4.mul(F4AlphaPlusOne, F4AlphaPlusOne))

    print("Inverses table for F4")
    print(F4.mul_inv(F4One))
    print(F4.mul_inv(F4Alpha))
    print(F4.mul_inv(F4AlphaPlusOne))
