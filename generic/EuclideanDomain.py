from abc import abstractmethod
from generic.UFD import UFD


class EuclideanDomain(UFD):
    def div(self, elem1, elem2):
        return divmod(elem1, elem2)[0]

    def rem(self, elem1, elem2):
        return divmod(elem1, elem2)[1]

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
    def extended_gcd(self, a, b):
        r, rp = a, b
        s, sp = self.mul_id(), self.add_id()
        t, tp = self.add_id(), self.mul_id()
        while rp != self.add_id():
            q, rpp = self.div(r, rp), self.rem(r, rp)
            r, s, t, rp, sp, tp = rp, sp, tp, rpp, self.add(s, self.add_inv(self.mul(sp, q))), self.add(t, self.add_inv(self.mul(tp, q)))
        d = r
        return d, s, t
