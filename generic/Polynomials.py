from generic.Polynomial import Polynomial
from generic.Ring import Ring


class Polynomials(Ring):

    def __init__(self, base_ring):
        self.base_ring = base_ring

    def add_id(self):
        return Polynomial([self.base_ring.add_id()], self.base_ring)

    def mul_id(self):
        return Polynomial([self.base_ring.mul_id()], self.base_ring)

    def add_inv(self, elem):
        lis = list(map(lambda x: self.base_ring.add_inv(x), elem.coefs))
        return Polynomial(lis, self.base_ring)

    def add(self, elem1, elem2):
        lis = list(map(lambda x1, x2: self.base_ring.add(x1, x2), elem1.coefs, elem2.coefs))
        lis.extend(elem1.coefs[len(lis):])
        lis.extend(elem2.coefs[len(lis):])
        return Polynomial(lis, self.base_ring)

    def mul(self, elem1, elem2):
        res = [0] * (elem1.degree() + elem2.degree() + 1)
        for in1, e1 in enumerate(elem1.coefs):
            for in2, e2 in enumerate(elem2.coefs):
                res[in1 + in2] += self.base_ring.mul(e1, e2)
        return Polynomial(res, self.base_ring)
