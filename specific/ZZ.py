from generic.EuclideanDomain import EuclideanDomain


class ZZ(EuclideanDomain):
    def rem(self, elem1, elem2):
        return elem1 % elem2

    def add_id(self):
        return 0

    def mul_id(self):
        return 1

    def add_inv(self, elem):
        return -elem

    def add(self, elem1, elem2):
        return elem1 + elem2

    def div_mod(self, elem1, elem2):
        return elem1 // elem2, elem1 % elem2

    def mul(self, elem1, elem2):
        return elem1 * elem2

    def mul_scalar(self, elem, scalar):
        return elem*scalar
