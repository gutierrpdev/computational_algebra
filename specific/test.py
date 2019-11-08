from generic.FiniteField import FiniteField
from specific.ZZ import ZZ
from generic.Polynomial import Polynomial
from generic.Polynomials import Polynomials
from generic.Ring import Ring
from specific.Zp import Zp


# ZZ = ZZ()
# print(ZZ.rem(5, 3))
# print(ZZ.add_id())
# print(ZZ.mul_id())
# print(ZZ.mul(7, -2))
# print(ZZ.gcd(178, 4))
# print(ZZ.gcd(7, -5))
# print(ZZ.extended_gcd(15, 25))
# print(ZZ.extended_gcd(25, 15))
# print(ZZ.extended_gcd(7, 25))
# PZ = Polynomials(ZZ)
# print(PZ.add_inv(Polynomial([1, 2], ZZ)).coefs)
# print(PZ.add(Polynomial([1, 2], ZZ), Polynomial([1, 2, 3], ZZ)).coefs)
# print(PZ.mul(Polynomial([1, 2], ZZ), Polynomial([1, 2, 3], ZZ)).coefs)

# Z3 = Zp(3)
# PZ3 = Polynomials(Z3)
# a, b = PZ3.div_mod(Polynomial([1, 2, 1], Z3), Polynomial([1, 1], Z3))
# print(a.coefs, b.coefs)
# a, b = PZ3.div_mod(Polynomial([1, 2, 1], Z3), Polynomial([0, 1], Z3))
# print(a.coefs, b.coefs)
# a, b = PZ3.div_mod(Polynomial([1, 2, 1], Z3), Polynomial([1], Z3))
# print(a.coefs, b.coefs)

Z2 = Zp(2)
F4 = FiniteField(2, 1, Polynomial([1, 1, 1], Z2))
# print(F4.rep(Polynomial([1, 0, 1], Z2)).coefs)
# print(F4.rep(Polynomial([0, 1], Z2)).coefs)
f = Polynomial([1, 0, 1], Z2)
print(F4.div_mod(Polynomial([0, 1], Z2), f)[0].coefs)
