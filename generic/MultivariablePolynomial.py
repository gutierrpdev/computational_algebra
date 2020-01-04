

class MultivariablePolynomial:

    def __init__(self, coefs_map, base_field):
        self.coefs = coefs_map
        self.base_field = base_field
        self.md = self.multidegree()

    # returns (lexicographic)
    @staticmethod
    def monomial_less(m1, m2):
        i = 0
        while i < len(m1) and m1[i] - m2[i] == 0:
            i += 1

        if m1[i] < m2[i]:
            return True
        else:
            return False

    def multidegree(self):
        i = 0
        max_m = None
        for elem in self.coefs:
            if max_m is None:
                max_m = elem
            elif self.monomial_less(max_m, elem):
                max_m = elem
        return max_m

    def leading_coef(self):
        return self.coefs[self.md]

    def leading_monomial(self):
        return self.md

    def leading_term(self):
        return MultivariablePolynomial({self.md: self.leading_coef()}, self.base_field)

    def __repr__(self):
        return str(self.coefs)


if __name__ == "__main__":
    from specific.ZZ import ZZ
    dic = {(1, 1, 2): 4, (3, 0, 0): 4, (0, 4, 0): -5, (1, 2, 1): 7}
    p = MultivariablePolynomial(dic, ZZ())
    print("Multidegree:", p.md)
    print("L.coef:", p.leading_coef())
    print("L.monomial:", p.leading_monomial())
    print("L.term:", p.leading_term())
