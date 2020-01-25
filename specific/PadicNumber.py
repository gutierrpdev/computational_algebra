from generic.EuclideanDomain import EuclideanDomain
from specific.ZZ import ZZ

# p-adic numbers are represented as tuples of the form (a,b), where n = a/b. This
# is the fraction representation.
class PadicNumber(EuclideanDomain):

    def __init__(self, p):
        self.p = p

    def mul_id(self):
        return 1, 1

    # Given two elements of the form (a,b), returns two lists containing the lists of
    # quotient and remainders when performing the Euclidean Algorithm. This elements are of the
    # form (a,b,c), with n = a/b * p**c.
    def p_adicEuclideanAlgorithm(self, elem1, elem2):
        dividendo = self.obtain_padic_representation(elem1)
        divisor = self.obtain_padic_representation(elem2)
        cociente, resto = self.div_mod(dividendo, divisor)
        qi = []
        nui = []

        # resto != 0
        while resto[0] != 0:
            dividendo, divisor = divisor, resto
            qi.append(cociente)
            nui.append(resto)
            cociente, resto = self.div_mod(dividendo, divisor)

        return qi, nui

    # Uses Euclidean Algorithm method to return the gcd.
    def gcd(self, elem1, elem2):
        qi, nui = self.p_adicEuclideanAlgorithm(elem1, elem2)
        return abs(nui[-1][0]), abs(nui[-1][1])

    def add_id(self):
        return 0, 1

    # Obtains primitive form of a fraction
    def simplify(self, elem):
        a, b = elem
        if a == 0:
            return 0, 1
        gcd = ZZ().gcd(a,b)
        return a // gcd, b // gcd

    def symmetric_representation(self, elem, mod):
        if abs(elem) > (mod-1) // 2:
            if elem < 0:
                return mod + elem
            else:
                return elem - mod
        else:
            return elem

    # Obtains division following the p-adic division algorithm. Input can be in fraction form or
    # s*p^v form.
    def div_mod(self, dividendo, divisor):
        if len(dividendo) == 2:
            a, b, c1 = self.obtain_padic_representation(dividendo)
        else:
            a, b, c1 = dividendo
        if len(divisor) == 2:
            c, d, c2 = self.obtain_padic_representation(divisor)
        else:
            c, d, c2 = divisor

        # |sigma|p < |tau|p
        if c2 < c1:
            return divisor, (0,1,0)

        # Otherwise
        s = (a*d) % self.p **(c2-c1+1)
        t = (b*c) % self.p **(c2-c1+1)

        _, t_inv, _ = ZZ().extended_gcd(t, self.p**(c2-c1+1))
        quotient = self.symmetric_representation((s*t_inv) % self.p**(c2-c1+1), self.p**(c2-c1+1))

        # nu = sigma - q*tau
        remainder = self.add(self.obtain_fraction_representation(dividendo),
                             self.add_inv(self.mul(self.obtain_fraction_representation(divisor),
                                                   (quotient, self.p**(c2-c1)))))

        # In case input was in fraction form, output is given in the same format
        if len(dividendo) == 2:
            return (quotient, self.p**(c2-c1)), remainder
        return (quotient, 1, c1-c2), self.obtain_padic_representation(remainder)

    # Given a p-adic number of the form (a,b,c), returns a/b * p^c
    def obtain_fraction_representation(self, elem):
        if len(elem) == 2:
            return elem
        a, b, c = elem
        if c < 0:
            b *= self.p**(-c)
        elif c > 0:
            a *= self.p**c
        return a, b

    def add(self, elem1, elem2):
        a, b = elem1
        c, d = elem2
        num = a*d + b*c
        div = b*d
        return self.simplify((num,div))

    def add_inv(self, elem):
        return -elem[0], elem[1]

    def mul(self, elem1, elem2):
        a, b = elem1
        c, d = elem2
        num = a*c
        div = b*d
        return self.simplify((num,div))

    # Given a fraction of form (a,b), with n = a/b; obtains p-adic representation of the form
    # n = c / d * p^v, where gcd(c,p) = 1 and gcd(d,p) = 1.
    def obtain_padic_representation(self, number):
        if len(number) == 3:
            return number
        a,b = number
        c1, c2 = 0,0
        while a != 0 and a % self.p == 0:
            c1 += 1
            a //= self.p
        while b != 0 and b % self.p == 0:
            c2 += 1
            b //= self.p
        return a, b, c1 - c2

if __name__ == "__main__":
    p7 = PadicNumber(7)
    n1 = 181625, 11
    n2 = 10555, 2
    print(p7.p_adicEuclideanAlgorithm(n1, n2))
    print(p7.gcd(n1, n2))