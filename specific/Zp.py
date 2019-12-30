import random

from generic.Field import Field
from specific.ZZ import ZZ


class Zp(Field):

    def __init__(self, p):
        self.p = p
        self.k = 1

    def mul_inv(self, elem):
        _, d, _ = ZZ().extended_gcd(elem, self.p)
        return d % self.p

    def add_id(self):
        return 0

    def mul_id(self):
        return 1

    def add_inv(self, elem):
        return (self.p - (elem % self.p)) % self.p

    def add(self, elem1, elem2):
        return (elem1 + elem2) % self.p

    def mul(self, elem1, elem2):
        return (elem1 * elem2) % self.p

    def mul_scalar(self, elem, scalar):
        return (elem * scalar) % self.p

    def cardinality(self):
        return self.p

    def random_elem(self):
        return random.randint(0, self.p-1)

    def __eq__(self, other):
        if isinstance(other, Zp):
            return self.p == other.p
        return NotImplemented
