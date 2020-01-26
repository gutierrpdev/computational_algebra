import copy
from functools import reduce
from itertools import combinations
from math import ceil, sqrt, floor, log2, log

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
                g_star = self.obtain_factors(subset, g_i, b, p)
                h_star = self.obtain_factors(T.difference(subset), g_i, b, p)
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

    def obtain_factors(self, T, g_i, b, p):
        zx = ZX()
        g = zx.mul_id()
        for index in T:
            g = zx.mul(g_i[index], g)
        g_star = zx.mul_scalar(g, b)
        g_star_module_p = Polynomial([x % p for x in g_star.coefs], ZX())
        return self.symmetric_form(g_star_module_p, p)

    # Modern Computer Algebra, page 446
    @staticmethod
    def hensel_step(m, f, g, h, s, t):
        # step 1
        zx = ZX()
        e = zx.add(f, zx.add_inv(zx.mul(g, h)))
        e = zx.mod_elem(e, m * m)
        q = zx.mul(s, e)
        quot, rem = zx.extended_synthetic_division(q.coefs, h.coefs)
        q = zx.mod_elem(Polynomial(quot, ZZ()), m * m)
        r = zx.mod_elem(Polynomial(rem, ZZ()), m * m)
        g_star = zx.add(g, zx.mul(t, e))
        g_star = zx.add(g_star, zx.mul(q, g))
        g_star = zx.mod_elem(g_star, m * m)
        h_star = zx.mod_elem(zx.add(h, r), m * m)

        # step 2
        b = zx.mul(s, g_star)
        b = zx.add(b, zx.mul(t, h_star))
        b = zx.add(b, zx.add_inv(zx.mul_id()))
        b = zx.mod_elem(b, m * m)
        c = zx.mul(s, b)
        quot, rem = zx.extended_synthetic_division(c.coefs, h_star.coefs)
        c = zx.mod_elem(Polynomial(quot, ZZ()), m * m)
        d = zx.mod_elem(Polynomial(rem, ZZ()), m * m)
        s_star = zx.add(s, zx.add_inv(d))
        s_star = zx.mod_elem(s_star, m * m)
        t_star = zx.add(t, zx.add_inv(zx.mul(t, b)))
        t_star = zx.add(t_star, zx.add_inv(zx.mul(c, g_star)))
        t_star = zx.mod_elem(t_star, m * m)
        return g_star, h_star, s_star, t_star

    # Modern Computer Algebra, page 450
    def multifactor_hensel_lifting(self, p, f, fr, l):
        zx = ZX()
        r = len(fr)
        # step 1
        if r == 1:
            # extended gcd in ZZ
            if f.get_leading_coef() == 1:
                return [f]

            pl = p ** l
            _, w, _ = ZZ().extended_gcd(f.get_leading_coef(), pl)
            return [zx.mod_elem(zx.mul_scalar(f, w), pl)]

        # step 2
        k = floor(r / 2)
        d = ceil(log2(l))

        # step 3
        g = reduce(zx.mul, fr[0:k])
        g = zx.mul_scalar(g, f.get_leading_coef())
        g = zx.mod_elem(g, p)
        h = reduce(zx.mul, fr[k:r])
        h = zx.mod_elem(h, p)

        # step 4: ZZ/<p> is always a field (namely Zp), since p is prime
        zpx = PolynomialsOverField(Zp(p))
        coef, s, t = zpx.extended_gcd(g, h)
        s = zpx.mul_scalar(s, Zp(p).mul_inv(coef.coefs[0]))
        t = zpx.mul_scalar(t, Zp(p).mul_inv(coef.coefs[0]))

        # step 5
        m = p
        for j in range(d):
            # step 6
            g, h, s, t = self.hensel_step(m, f, g, h, s, t)
            m **= 2

        _k1 = self.multifactor_hensel_lifting(p, g, fr[0:k], l)
        _k2 = self.multifactor_hensel_lifting(p, h, fr[k:r], l)
        _k1.extend(_k2)
        return _k1

    # Page 453, Modern Computer Algebra
    def factorize_prime_power(self):
        # step 1
        if self.n == 1:
            return [self.f]
        b = self.f.get_leading_coef()
        B = ceil(sqrt(self.n + 1) * pow(2, self.n) * self.A * b)
        C = (self.n + 1) ** (2 * self.n) * (self.A ** (2 * self.n - 1))
        gamma = ceil(2 * log2(C))

        # step 2
        p = ZZ().obtain_random_prime_in_range(3, ceil(2 * gamma * log(gamma)))
        f_bar = Polynomial(self.f.coefs, Zp(p))
        f_bar_der = f_bar.derivative()
        fx = PolynomialsOverField(Zp(p))

        while b % p == 0 or fx.gcd(f_bar, f_bar_der).degree() > 1:
            p = ZZ().obtain_random_prime_in_range(2, ceil(2 * gamma * log(gamma)))
            f_bar = Polynomial(self.f.coefs, Zp(p))
            f_bar_der = f_bar.derivative()
        l = ceil(log(2 * B + 1, p))

        # step 3
        h_i = self.modular_factorization(f_bar, p, fx)

        # step 4
        fact = self.multifactor_hensel_lifting(p, self.f, h_i, l)
        g_i = [self.symmetric_form(factor, p**l) for factor in fact]

        # step 5
        r = len(g_i)
        f_star = copy.deepcopy(self.f)
        T = set(range(r))
        s = 1
        G = []

        # step 6
        while 2 * s <= len(T):
            goto = False
            combination = list(combinations(T, s))
            for subset_tuple in combination:
                subset = set(subset_tuple)
                g_star = self.obtain_factors(subset, g_i, b, p**l)
                h_star = self.obtain_factors(T.difference(subset), g_i, b, p**l)
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


if __name__ == "__main__":
    g = Polynomial([-15, -17, -16, -1, 1], ZZ())
    fact = FactorizationZx(g)
    print(fact.factorize_big_prime())
    print(fact.factorize_prime_power())

    p_big = Polynomial([-12, -24, -18, -2, 8, 3, 0, 2, 1], ZZ())
    # print(FactorizationZx(p_big).factorize_big_prime())
    print(FactorizationZx(p_big).factorize_prime_power())

    m = 5
    f = Polynomial([-1, 0, 0, 0, 1], ZZ())
    g = Polynomial([-2, -1, 2, 1], ZZ())
    h = Polynomial([-2, 1], ZZ())
    s = Polynomial([-2], ZZ())
    t = Polynomial([-1, -2, 2], ZZ())

    print(fact.hensel_step(m, f, g, h, s, t))

    f1 = Polynomial([-1, 1], ZZ())
    f2 = Polynomial([-2, 1], ZZ())
    f3 = Polynomial([2, 1], ZZ())
    f4 = Polynomial([1, 1], ZZ())
    print(fact.multifactor_hensel_lifting(5, f, [f1, f2, f3, f4], 4))
