import copy

from generic.Polynomial import Polynomial
from generic.PolynomialsOverField import PolynomialsOverField


class CantorZassenhaus:

    def __init__(self, f):
        # field = Fq, f = f, fx = F[X]
        self.field = f.base_ring
        self.f = f
        self.fx = PolynomialsOverField(f.base_ring)
        self._l = f.degree()

    def distinct_degree_factorization(self):
        _L = []
        _X = Polynomial([self.field.add_id(), self.field.mul_id()], self.field)
        h = self.fx.rem(_X, self.f)
        k = 0
        fc = copy.deepcopy(self.f)

        while fc != self.fx.mul_id():
            h = self.fx.rem(self.fx.exp(h, self.field.cardinality()), fc)
            k += 1
            aux = self.fx.add(h, self.fx.add_inv(_X))
            g = self.fx.gcd(aux, fc)
            if g.degree() > 1:
                _L.append((g, k))
                fc = self.fx.div(fc, g)
                h = self.fx.rem(h, fc)
        return _L

    # generate a random polynomial in F[X]/(h)
    def random_polynomial_h(self, h):
        res = []
        for i in range(h.degree()):
            res.append(self.field.random_elem())
        pol = Polynomial(res, self.field)
        return self.fx.rem(pol, h)

    def mk(self, k):
        w = self.field.k
        res = [self.field.add_id()] * (2 ** (w * k - 1) + 1)
        acc = 1
        for i in range(w * k):
            res[acc] = self.field.mul_id()
            acc *= 2
        return Polynomial(res, self.field)

    # f is a monic polynomial of degree l and k the number of times
    def equal_degree_factorization(self, pol, k):
        mk = self.mk(k)
        fc = copy.deepcopy(pol)
        r = pol.degree() // k
        # print("r: ", r)
        _H = [fc]
        # print("_H: ", _H)
        while len(_H) < r:
            _Hp = []
            for elem in _H:
                alpha = self.random_polynomial_h(elem)
                eval_mk = self.fx.evaluate_polynomial(mk, alpha)
                eval_mk = self.fx.rem(eval_mk, fc)
                d = self.fx.gcd(eval_mk, elem)
                if d.degree() == 0 or d == elem:
                    _Hp.append(elem)
                else:
                    # print("d: ", d)
                    # print("elem: ", elem)
                    _Hp.append(d)
                    # print("Div :", self.fx.div(elem, d))
                    _Hp.append(self.fx.div(elem, d))
            _H = _Hp
        return _H

    def compute(self):
        step1 = self.distinct_degree_factorization()
        # print("Step 1 result: ", step1)
        # for pol in step1:
        #    print(self.fx.obtain_monic_representation(pol[0]))
        res = []
        for elem in step1:
            # print("elem for: ", elem)
            sol = self.equal_degree_factorization(elem[0], elem[1])
            # print(sol)
            res.extend(sol)
        return res


if __name__ == "__main__":
    from generic.FiniteField import FiniteField
    from specific.Zp import Zp

    F4 = FiniteField(2, 2, Polynomial([1, 1, 1], Zp(2)))
    F4One = F4.mul_id()
    F4Zero = F4.add_id()
    F4Alpha = Polynomial([0, 1], Zp(2))
    F4AlphaPlusOne = Polynomial([1, 1], Zp(2))
    f = Polynomial([F4Alpha, F4AlphaPlusOne, F4One], F4)
    print(".............")
    print("Compute factorization through Berkelamp's algorithm")
    print(CantorZassenhaus(f).compute())
