from generic.Field import Field
from generic.Ring import Ring
from specific.ZZ import ZZ


class Zn(Ring):

    def __init__(self, n):
        self.n = n

    def add_id(self):
        return 0

    def mul_id(self):
        return 1

    def add_inv(self, elem):
        return (self.n - (elem % self.n)) % self.n

    def add(self, elem1, elem2):
        return (elem1 + elem2) % self.n

    def mul(self, elem1, elem2):
        return (elem1 * elem2) % self.n

    def mul_scalar(self, elem, scalar):
        return (elem * scalar) % self.n

    def cardinality(self):
        return self.p

    def __eq__(self, other):
        if isinstance(other, Zn):
            return self.n == other.n
        return NotImplemented
