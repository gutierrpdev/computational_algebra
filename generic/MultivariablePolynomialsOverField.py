import copy

from generic.MultivariablePolynomial import MultivariablePolynomial


class MultivariablePolynomialsOverField:

    def __init__(self, field):
        self.field = field

    def add(self, p1, p2):
        res = copy.deepcopy(p1.coefs)
        for elem in p2.coefs:
            if elem in res:
                suma = self.field.add(res[elem], p2.coefs[elem])
                if suma == self.field.add_id():
                    res.pop(elem)
                else:
                    res[elem] = suma
            else:
                res[elem] = p2.coefs[elem]
        return MultivariablePolynomial(res, self.field)

    # takes two Multivariable Polynomials with only one term each and returns the division of the two, if possible.
    # otherwise, returns None.
    def div_term(self, p1, p2):
        # only allow division of single terms
        if len(p1.coefs) != len(p2.coefs) or len(p1.coefs) != 1:
            return None
        p1_item = list(p1.coefs.items())[0]
        p2_item = list(p2.coefs.items())[0]
        res = [0] * len(p1_item[0])

        for i in range(len(p1_item[0])):
            if p1_item[0][i] >= p2_item[0][i]:
                res[i] = p1_item[0][i] - p2_item[0][i]
            else:
                return None
        return MultivariablePolynomial({tuple(res): self.field.div(p1_item[1], p2_item[1])}, self.field)

    # takes two single-term Multivariable Polynomials and returns the product of the two
    def mul_term(self, p1, p2):
        # only allow single-term Polynomials here
        if len(p1.coefs) != len(p2.coefs) or len(p1.coefs) != 1:
            return None
        p1_item = list(p1.coefs.items())[0]
        p2_item = list(p2.coefs.items())[0]
        res = [0] * len(p1_item[0])

        for i in range(len(p1_item[0])):
            res[i] = p1_item[0][i] + p2_item[0][i]
        return MultivariablePolynomial({tuple(res): self.field.mul(p1_item[1], p2_item[1])}, self.field)

    # takes two Multivariable-Polynomials and returns the product of the two.
    def mul(self, p1, p2):
        # result of the multiplication.
        res = MultivariablePolynomial({}, self.field)
        for elem1 in p1.coefs:
            for elem2 in p2.coefs:
                term1 = MultivariablePolynomial({elem1: p1.coefs[elem1]}, self.field)
                term2 = MultivariablePolynomial({elem2: p2.coefs[elem2]}, self.field)
                res = self.add(res, self.mul_term(term1, term2))
        return res

    # division algorithm from Modern Computer Algebra, page 599.
    def div(self, f, fs):
        r = MultivariablePolynomial({}, self.field)
        p = copy.deepcopy(f)
        qs = [MultivariablePolynomial({}, self.field)] * len(fs)

        while p != MultivariablePolynomial({}, self.field):
            i = 0
            quot = None
            while i < len(fs) and (quot is None):
                quot = self.div_term(p.leading_term(), fs[i].leading_term())
                i += 1

            # nothing divides lt(p)
            if quot is None:
                r = self.add(r, p.leading_term())
                p = self.add(p, self.add_inv(p.leading_term()))
            else:
                i -= 1
                qs[i] = self.add(qs[i], self.div_term(p.leading_term(), fs[i].leading_term()))
                aux = self.mul(self.div_term(p.leading_term(), fs[i].leading_term()), fs[i])
                p = self.add(p, self.add_inv(aux))
        return qs, r

    def add_inv(self, pol):
        res = copy.deepcopy(pol)
        for elem in res.coefs:
            res.coefs[elem] = self.field.add_inv(res.coefs[elem])
        return res

    def s_poly(self, g, h):
        alpha = g.multidegree()
        beta = h.multidegree()
        sigma = tuple([max(alpha[i], beta[i]) for i in range(len(alpha))])
        x = MultivariablePolynomial({sigma: self.field.mul_id()}, self.field)
        t1 = self.mul(self.div_term(x, g.leading_term()), g)
        t2 = self.mul(self.div_term(x, h.leading_term()), h)
        return self.add(t1, self.add_inv(t2))

    # Buchberger Algorithm from Modern Computer Algebra page 610
    def groebner_basis(self, fs):
        _g = fs

        while True:
            _s = []
            for i in range(len(_g)-1):
                for j in range(i + 1, len(_g)):
                    aux = self.s_poly(_g[i], _g[j])
                    r = self.div(aux, _g)[1]
                    if len(_g) == 5 and r != MultivariablePolynomial({}, self.field):
                        print(self.div(aux, _g), "r", i, j)
                    if r != MultivariablePolynomial({}, self.field):
                        _s.append(r)
            if not _s:
                return _g
            else:
                _g.extend(_s)

    # check whether f is in <f1, ..., fs>. Theorem 21.28, page 605, MCA
    def in_ideal(self, f, fs):
        _g = self.groebner_basis(fs)
        return self.div(f, fs)[1] == MultivariablePolynomial({}, self.field)


if __name__ == "__main__":
    from specific.ZZ import ZZ
    dic1 = {(1, 1, 2): 4, (3, 0, 0): 4, (0, 4, 0): -5, (1, 2, 1): 7}
    p1 = MultivariablePolynomial(dic1, ZZ())
    dic2 = {(1, 1, 2): 4, (3, 0, 0): 4, (0, 5, 0): -3, (1, 2, 1): -7}
    p2 = MultivariablePolynomial(dic2, ZZ())
    print("p1 = ", p1)
    print("p2 = ", p2)
    print("p1 + p2 = ", MultivariablePolynomialsOverField(ZZ()).add(p1, p2))
    print("p1 * p2 = ", MultivariablePolynomialsOverField(ZZ()).mul(p1, p2))

    dic3 = {(2, 4, 6, 8): 25}
    p3 = MultivariablePolynomial(dic3, ZZ())
    dic4 = {(2, 2, 3, 2): 5}
    p4 = MultivariablePolynomial(dic4, ZZ())
    print("p3 = ", p3)
    print("p4 = ", p4)
    print("p3 / p4 = ", MultivariablePolynomialsOverField(ZZ()).div_term(p3, p4))
    print("p3 * p4 (term) = ", MultivariablePolynomialsOverField(ZZ()).mul_term(p3, p4))
    print("p3 * p4 (mul) = ", MultivariablePolynomialsOverField(ZZ()).mul(p3, p4))

    dic5 = {(1, 0, 0): 2}
    p5 = MultivariablePolynomial(dic5, ZZ())
    dic6 = {(1, 0, 0): 1, (0, 1, 0): 1, (0, 0, 1): 1}
    p6 = MultivariablePolynomial(dic6, ZZ())
    print("p5 = ", p5)
    print("p6 = ", p6)
    print("p5 * p6 (mul) = ", MultivariablePolynomialsOverField(ZZ()).mul(p5, p6))

    # example from page 599, modern computer algebra
    from specific.QQ import QQ
    dic7 = {(2, 1): 1, (1, 2): 1, (0, 2): 1}
    dic8 = {(1, 1): 1, (0, 0): -1}
    dic9 = {(0, 2): 1, (0, 0): -1}
    f = MultivariablePolynomial(dic7, QQ())
    f1 = MultivariablePolynomial(dic8, QQ())
    f2 = MultivariablePolynomial(dic9, QQ())

    print("f = ", f)
    print("f1 = ", f1)
    print("f2 = ", f2)
    print("f / f1, f2 = ", MultivariablePolynomialsOverField(QQ()).div(f, [f1, f2]))

    dic_g = {(3, 0): 1, (1, 1): -2}
    dic_h = {(2, 1): 1, (0, 2): -2, (1, 0): 1}
    g = MultivariablePolynomial(dic_g, QQ())
    h = MultivariablePolynomial(dic_h, QQ())
    print(MultivariablePolynomialsOverField(QQ()).groebner_basis([g, h]))
