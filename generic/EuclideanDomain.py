from abc import abstractmethod
from generic.UFD import UFD
from functools import reduce


class EuclideanDomain(UFD):

    def div(self, elem1, elem2):
        return self.div_mod(elem1, elem2)[0]

    def rem(self, elem1, elem2):
        return self.div_mod(elem1, elem2)[1]

    @abstractmethod
    def div_mod(self, elem1, elem2):
        pass

    # Euclidean algorithm for greatest common divisor
    def gcd(self, a, b):
        while b != self.add_id():
            r = self.rem(a, b)
            a = b
            b = r
        return a

    # Extended euclidean algorithm for gcd + Bezout's coefficients
    # a*s + b*t = d
    def extended_gcd(self, a, b):
        r, rp = a, b
        s, sp = self.mul_id(), self.add_id()
        t, tp = self.add_id(), self.mul_id()
        while rp != self.add_id():
            q, rpp = self.div(r, rp), self.rem(r, rp)
            r, s, t, rp, sp, tp = rp, sp, tp, rpp, self.add(s, self.add_inv(self.mul(sp, q))), self.add(t, self.add_inv(self.mul(tp, q)))
        d = r
        return d, s, t

    # Chinese Remainder Theorem solver for a system of type: x congruent to a_i mod n_i, 0 <= i < k
    # under the assumption that for every i, 0 <= deg (a_i) < deg (n_i)
    # Note that this is always possible since deg(z) is well defined for any z within a euclidean domain.
    def chinese_remainder_theorem(self, a_vec, n_vec):
        # multiply all elements in n_vec
        n = reduce(self.mul, n_vec)
        e_vec = [None] * len(n_vec)
        for i in range(0, len(n_vec)):
            n_star = self.div(n, n_vec[i])
            b_i = self.rem(n_star, n_vec[i])
            t_i = self.extended_gcd(b_i, n_vec[i])[1]
            e_vec[i] = self.mul(n_star, t_i)
        return self.rem(reduce(self.add, map(lambda t: self.mul(t[0], t[1]), zip(a_vec, e_vec))), n)

    # return gcd of a list of elements in Euclidean Domain.
    def multiple_gcd(self, elem_list):
        if len(elem_list) == 1:
            return elem_list[0]
        else:
            return self.gcd(elem_list[0], self.multiple_gcd(elem_list[1:]))
