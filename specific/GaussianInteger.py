from generic.EuclideanDomain import EuclideanDomain


class GaussianIntegers(EuclideanDomain):
    def div_mod(self, elem1, elem2):
        pass

    def add_id(self):
        return 0, 0

    def mul_id(self):
        return 1, 0

    def add_inv(self, elem):
        a, b = elem
        return -a, -b

    def add(self, elem1, elem2):
        a1, b1 = elem1
        a2, b2 = elem2
        return a1 + a2, b1 + b2

    def mul(self, elem1, elem2):
        a1, b1 = elem1
        a2, b2 = elem2
        return a1 * a2 - b1 * b2, a1 * b2 + a2 * b1
