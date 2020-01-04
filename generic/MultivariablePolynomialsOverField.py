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


if __name__ == "__main__":
    from specific.ZZ import ZZ
    dic1 = {(1, 1, 2): 4, (3, 0, 0): 4, (0, 4, 0): -5, (1, 2, 1): 7}
    p1 = MultivariablePolynomial(dic1, ZZ())
    dic2 = {(1, 1, 2): 4, (3, 0, 0): 4, (0, 5, 0): -3, (1, 2, 1): -7}
    p2 = MultivariablePolynomial(dic2, ZZ())
    print("p1 = ", p1)
    print("p2 = ", p2)
    print("p1 + p2 = ", MultivariablePolynomialsOverField(ZZ()).add(p1, p2))
