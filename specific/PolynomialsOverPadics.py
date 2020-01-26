from generic.PolynomialsOverField import PolynomialsOverField
from specific.PadicNumber import PadicNumber


# Not much to do here, this class is mostly intended for testing purposes.
class PolynomialsOverPadics(PolynomialsOverField):

    def __init__(self, p):
        self.base_field = PadicNumber(p)


if __name__ == "__main__":
    from generic.Polynomial import Polynomial
    pp2 = PolynomialsOverPadics(2)
    p2 = PadicNumber(2)
    f1 = Polynomial([p2.mul_id(), p2.mul_id()], p2)
    f2 = Polynomial([p2.mul_id(), p2.mul_id(), p2.mul_id()], p2)
    f12 = pp2.mul(f1, f2)
    f11 = pp2.mul(f1, f1)
    f22 = pp2.mul(f2, f2)
    print("f1: ", f1)
    print("f2: ", f2)
    print("f11: ", f11)
    print("f12: ", f12)
    print("f22: ", f22)
    print("-------------------------------")
    res1 = pp2.gcd(f11, f12)
    print("gcd f11, f12: ", res1)
    print("divided by f1: ", pp2.div(res1, f1))
    print("divided by f2: ", pp2.div(res1, f2))
    print("-------------------------------")
    res2 = pp2.gcd(f11, f22)
    print("gcd f11, f22: ", res2)
    print("divided by f1: ", pp2.div(res2, f1))
    print("divided by f2: ", pp2.div(res2, f2))
    print("-------------------------------")
    res3 = pp2.gcd(f12, f2)
    print("gcd f12, f2: ", res3)
    print("divided by f1: ", pp2.div(res3, f1))
    print("divided by f2: ", pp2.div(res3, f2))
    print("-------------------------------")
    res4 = pp2.gcd(f22, pp2.mul(f22, f2))
    print("gcd f22, f222: ", res4)
    print("divided by f1: ", pp2.div(res4, f1))
    print("divided by f2: ", pp2.div(res4, f2))
    print("divided by f22: ", pp2.div(res4, f22))
    print("-------------------------------")
