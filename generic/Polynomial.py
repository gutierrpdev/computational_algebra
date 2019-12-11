class Polynomial:

    def __init__(self, coefs, base_ring):
        i = len(coefs) - 1
        while i > 0 and coefs[i] == base_ring.add_id():
            i -= 1
        self.coefs = coefs[:i+1]
        self.base_ring = base_ring

    def degree(self):
        return len(self.coefs) - 1

    def get_coefs(self):
        return self.coefs

    def get_base_ring(self):
        return self.base_ring

    def get_leading_coef(self):
        return self.coefs[self.degree()]

    def derivative(self):
        res = [self.base_ring.mul_scalar(self.coefs[i], i) for i in range(1, self.degree()+1)]
        i = self.degree() - 1
        while i > 0 and self.coefs[i] == self.base_ring.add_id():
            i -= 1
        return Polynomial(res[:i+1], self.base_ring)

    def __eq__(self, other):
        if isinstance(other, Polynomial):
            if len(self.coefs) != len(other.coefs):
                return False
            is_same = True
            for i in range(0, len(self.coefs)):
                is_same = is_same and self.coefs[i] == other.coefs[i]
            return is_same
        return NotImplemented

    def __str__(self):
        return str(self.coefs)

    def __repr__(self):
        # aux = ""
        # for i in range(0, len(self.coefs)):
        #     if self.coefs[i] != self.base_ring.add_id():
        #         aux += (self.coefs[i]).__str__() + "*X^" + str(i) + " + "
        # return aux
        return str(self.coefs)


if __name__ == "__main__":
    from specific.Zp import Zp
    from generic.FiniteField import FiniteField
    f = Polynomial([1, 1, 1, 1, 1], Zp(2))
    print(f.derivative())
    # should be X^2 + 1
    F4 = FiniteField(2, 2, Polynomial([1, 1, 1], Zp(2)))
    F4One = F4.mul_id()
    f = Polynomial([F4One, F4One, F4One, F4One, F4One], F4)
    print(f.derivative())
