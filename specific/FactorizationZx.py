import copy
from itertools import combinations
from math import ceil, sqrt

from generic.Polynomial import Polynomial
from generic.PolynomialsOverField import PolynomialsOverField
from specific.CantorZassenhaus import CantorZassenhaus
from specific.ZX import ZX
from specific.ZZ import ZZ
from specific.Zp import Zp


class FactorizationZx:

    def __init__(self, f):
        self.f = f
        self.n = f.degree()
        self.A = ZX().maxnorm(self.f)

    @staticmethod
    def symmetrical_representation_number(number, p):
        if number > (p - 1) // 2:
            return number - p
        else:
            return number

    def symmetric_form(self, pol, p):
        return Polynomial([self.symmetrical_representation_number(x, p) for x in pol.coefs], pol.base_ring)

    def modular_factorization(self, f_bar, p, fx):
        # f_bar_decomposition = Berlekamp(f_bar).compute()
        f_bar_decomposition = CantorZassenhaus(f_bar).compute()
        f_bar_monic_decomposition = [fx.obtain_monic_representation(factor) for factor in f_bar_decomposition]
        f_bar_symmetric_decomposition = [self.symmetric_form(factor, p) for factor in f_bar_monic_decomposition]
        return f_bar_symmetric_decomposition

    def factorize_big_prime(self):
        if self.n == 1:
            return [self.f]
        b = self.f.get_leading_coef()
        B = ceil(sqrt(self.n + 1) * pow(2, self.n) * self.A * b)

        p = ZZ().obtain_random_prime_in_range(2 * B + 1, 4 * B - 1)
        f_bar = Polynomial(self.f.coefs, Zp(p))
        f_bar_der = f_bar.derivative()
        fx = PolynomialsOverField(Zp(p))

        while fx.gcd(f_bar, f_bar_der).degree() > 1:
            p = ZZ().obtain_random_prime_in_range(2 * B + 1, 4 * B - 1)
            f_bar = Polynomial(self.f.coefs, Zp(p))
            f_bar_der = f_bar.derivative()

        g_i = self.modular_factorization(f_bar, p, fx)
        r = len(g_i)
        f_star = copy.deepcopy(self.f)
        T = set(range(r))
        s = 1
        G = []
        while 2 * s <= len(T):
            goto = False
            combination = list(combinations(T, s))
            for subset_tuple in combination:
                subset = set(subset_tuple)
                g_star = self.obtain_factors(subset, g_i, fx, b, p)
                h_star = self.obtain_factors(T.difference(subset), g_i, fx, b, p)
                if ZX().onenorm(g_star) * ZX().onenorm(h_star) <= B:
                    T = T.difference(subset)
                    G.append(ZX().primitive_part(g_star))
                    f_star = ZX().primitive_part(h_star)
                    b = f_star.get_leading_coef()
                    goto = True
                    break
            if not goto:
                s += 1
        G.append(f_star)
        return G

    def obtain_factors(self, T, g_i, fx, b, p):
        g = fx.mul_id()
        for index in T:
            g = fx.mul(g_i[index], g)
        g_star = fx.mul_field_scalar(g, b)
        g_star_module_p = Polynomial([x % p for x in g_star.coefs], ZX())
        return self.symmetric_form(g_star_module_p, p)


if __name__ == "__main__":
    g = Polynomial([-15, -17, -16, -1, 1], ZZ())
    fact = FactorizationZx(g)
    print(fact.factorize_big_prime())
