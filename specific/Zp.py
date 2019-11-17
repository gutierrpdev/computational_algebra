from generic.Field import Field
from specific.ZZ import ZZ


class Zp(Field):

    def __init__(self, p):
        self.p = p

    def mul_inv(self, elem):
        _, d, _ = ZZ().extended_gcd(elem, self.p)
        return d % self.p

    def add_id(self):
        return 0

    def mul_id(self):
        return 1

    def add_inv(self, elem):
        return self.p - (elem % self.p)

    def add(self, elem1, elem2):
        return (elem1 + elem2) % self.p

    def mul(self, elem1, elem2):
        return (elem1 * elem2) % self.p
