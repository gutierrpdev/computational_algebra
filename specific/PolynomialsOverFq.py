from generic.PolynomialsOverField import PolynomialsOverField
from generic.Polynomial import Polynomial
from specific.Zp import Zp


class PolynomialsOverFq(PolynomialsOverField):
    # checks whether element is irreducible over Fq[X] using the IPT algorithm.
    def is_irreducible(self, f):
        x = Polynomial([0, 1], self.base_field)
        h = self.rem(x, f)
        for k in range(1, f.degree()//2):
            aux = h
            # calculate h = h ** card, where card = cardinality of Fq
            for j in range(1, self.base_field.cardinality()):
                h = self.mul(h, aux)

            if self.gcd(self.add(h, self.add_inv(x)), f) != self.mul_id():
                return False
        return True


if __name__ == "__main__":
    Z2 = Zp(2)
    Z2X = PolynomialsOverFq(Z2)
    print(Z2X.is_irreducible(Polynomial([1, 0, 1], Z2)))
