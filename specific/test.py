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

Z11 = Zp(11)
print(Z11.mul_inv(5))
print(Z11.mul_inv(1))
Z4 = FiniteField(2, 2, Polynomial([1, 1, 1], Zp(2)))
print(Z4.mul_inv(Polynomial([1, 1], Zp(2))))

