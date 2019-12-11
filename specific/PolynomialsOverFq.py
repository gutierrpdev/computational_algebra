from generic.FiniteField import FiniteField
from generic.PolynomialsOverField import PolynomialsOverField
from generic.Polynomial import Polynomial
from specific.Zp import Zp


class PolynomialsOverFq(PolynomialsOverField):
    # checks whether element is irreducible over Fq[X] using the IPT algorithm.
    def is_irreducible(self, f):
        x = Polynomial([self.base_field.add_id(), self.base_field.mul_id()], self.base_field)
        h = self.rem(x, f)
        for k in range(1, f.degree()//2+1):
            aux = h
            # calculate h = h ** card, where card = cardinality of Fq
            for j in range(1, self.base_field.cardinality()):
                h = self.mul(h, aux)
            h = self.rem(h, f)
            gcd = self.gcd(self.add(h, self.add_inv(x)), f)
            if self.div(gcd, Polynomial([gcd.get_leading_coef()], self.base_field)) != self.mul_id():
                return False
        return True

    # perform the SFD algorithm to compute a square free decomposition of a monic polynomial f.
    def square_free_decomposition(self, f):
        _L = []
        s = 1
        while True:
            j = 1
            g = self.div(f, self.gcd(f, f.derivative()))
            print("g: ", g)
            while g != self.mul_id():
                print("g_loop: ", g)
                f = self.div(f, g)
                h = self.gcd(f, g)
                m = self.div(g, h)
                if m != self.mul_id():
                    _L.append((m, j*s))
                g = h
                j += 1
            if f != self.mul_id():
                print("f_loop: ", f)
                f = self.pth_root(f)
                print("p-th root: ", f)
                s = self.base_field.p * s

            if f == self.mul_id():
                break
        return _L

    def pth_root(self, f):
        # f = g(X^p), with g in PolynomialsOverFq
        g = []
        i = 0
        # x^4 + x^2 + 1 -> x^2 + x + 1 (p = 2)
        while i <= f.degree():
            g.append(f.coefs[i])
            i += self.base_field.p

        print("g*: ", g)
        h = [self.base_field.exp(a, self.base_field.p**(self.base_field.k - 1)) for a in g]
        return Polynomial(h, f.base_ring)


if __name__ == "__main__":
    Z2 = Zp(2)
    Z2X = PolynomialsOverFq(Z2)
    print(Z2X.is_irreducible(Polynomial([1, 0, 1], Z2)))
    # false
    print(Z2X.is_irreducible(Polynomial([1, 1, 1, 1, 1], Z2)))
    # true
    print(Z2X.is_irreducible(Polynomial([0, 1, 1, 1, 1], Z2)))
    # false
    print(Z2X.is_irreducible(Polynomial([1, 1, 0, 1], Z2)))
    # true
    print(Z2X.is_irreducible(Polynomial([1, 0, 1, 1], Z2)))
    # true
    print(Z2X.is_irreducible(Polynomial([1, 0, 0, 0, 1, 1], Z2)))
    # false
    F16 = FiniteField(2, 4, Polynomial([1, 1, 0, 0, 1], Z2))
    F16X = PolynomialsOverFq(F16)
    One = F16.mul_id()
    Zero = F16.add_id()
    Two = F16.add(One, One)
    # print(F16X.is_irreducible(Polynomial([One, Zero, One, Zero, One], F16)))
    # false
    F4 = FiniteField(2, 2, Polynomial([1, 1, 1], Z2))
    F4X = PolynomialsOverFq(F4)
    F4One = F4.mul_id()
    F4Zero = F4.add_id()
    F4Alpha = Polynomial([0, 1], Z2)
    F4AlphaPlusOne = Polynomial([1, 1], Z2)

    # For every polynomial with degree 2 in F4, check if algorithm works
    print("Check if Polynomials over F4 works...")
    print("----------------------------------------")
    print(F4X.is_irreducible(Polynomial([F4Zero, F4Zero, F4One], F4)))
    print(F4X.is_irreducible(Polynomial([F4One, F4Zero, F4One], F4)))
    print(F4X.is_irreducible(Polynomial([F4Alpha, F4Zero, F4One], F4)))
    print(F4X.is_irreducible(Polynomial([F4AlphaPlusOne, F4Zero, F4One], F4)))
    print(F4X.is_irreducible(Polynomial([F4Zero, F4One, F4One], F4)))
    print(F4X.is_irreducible(Polynomial([F4One, F4One, F4One], F4)))
    print(F4X.is_irreducible(Polynomial([F4Alpha, F4One, F4One], F4)))
    print(F4X.is_irreducible(Polynomial([F4AlphaPlusOne, F4One, F4One], F4)))
    print(F4X.is_irreducible(Polynomial([F4Zero, F4Alpha, F4One], F4)))
    print(F4X.is_irreducible(Polynomial([F4One, F4Alpha, F4One], F4)))
    print(F4X.is_irreducible(Polynomial([F4Alpha, F4Alpha, F4One], F4)))
    print(F4X.is_irreducible(Polynomial([F4AlphaPlusOne, F4Alpha, F4One], F4)))
    print(F4X.is_irreducible(Polynomial([F4Zero, F4AlphaPlusOne, F4One], F4)))
    print(F4X.is_irreducible(Polynomial([F4One, F4AlphaPlusOne, F4One], F4)))
    print(F4X.is_irreducible(Polynomial([F4Alpha, F4AlphaPlusOne, F4One], F4)))
    print(F4X.is_irreducible(Polynomial([F4AlphaPlusOne, F4AlphaPlusOne, F4One], F4)))
    # False
    # False
    # False
    # False
    # False
    # False
    # True
    # True
    # False
    # True
    # True
    # False
    # False
    # True
    # False
    # True

    print("-------------------------------------")
    print("Check if SFD Algorithm works properly...")
    F4X = PolynomialsOverFq(F4)
    li = F4X.square_free_decomposition(Polynomial([F4One, F4Zero, F4One], F4))
    print(li)
    li = F4X.square_free_decomposition(Polynomial([F4Zero, F4Zero, F4One], F4))
    print(li)
    li = F4X.square_free_decomposition(Polynomial([F4AlphaPlusOne, F4Zero, F4One], F4))
    print(li)

