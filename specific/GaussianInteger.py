from generic.EuclideanDomain import EuclideanDomain


class GaussianIntegers(EuclideanDomain):

    # Proof of degree's validity
    # http://fermatslasttheorem.blogspot.com/2005/06/division-algorithm-for-gaussian.html
    def div_mod(self, elem1, elem2):
        a, b = elem1
        c, d = elem2
        norm2 = c ** 2 + d ** 2
        r, s = round((a * c + b * d) / norm2), round((b * c - a * d) / norm2)
        coc = (r, s)
        rem = self.add(elem1, self.add_inv(self.mul(coc, elem2)))
        return coc, rem

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
