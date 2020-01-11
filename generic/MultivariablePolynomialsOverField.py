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
            if p1_item[0][i] % p2_item[0][i] == 0:
                res[i] = p1_item[0][i] // p2_item[0][i]
            else:
                return None
        return MultivariablePolynomial({tuple(res): self.field.div(p1_item[1], p2_item[1])}, self.field)


if __name__ == "__main__":
    from specific.ZZ import ZZ
    dic1 = {(1, 1, 2): 4, (3, 0, 0): 4, (0, 4, 0): -5, (1, 2, 1): 7}
    p1 = MultivariablePolynomial(dic1, ZZ())
    dic2 = {(1, 1, 2): 4, (3, 0, 0): 4, (0, 5, 0): -3, (1, 2, 1): -7}
    p2 = MultivariablePolynomial(dic2, ZZ())
    print("p1 = ", p1)
    print("p2 = ", p2)
    print("p1 + p2 = ", MultivariablePolynomialsOverField(ZZ()).add(p1, p2))

    dic3 = {(2, 4, 6, 8): 25}
    p3 = MultivariablePolynomial(dic3, ZZ())
    dic4 = {(2, 2, 3, 2): 5}
    p4 = MultivariablePolynomial(dic4, ZZ())
    print("p3 = ", p3)
    print("p4 = ", p4)
    print("p3 / p4 = ", MultivariablePolynomialsOverField(ZZ()).div_term(p3, p4))
