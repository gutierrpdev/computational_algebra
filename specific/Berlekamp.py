import copy
from functools import reduce

from generic.FiniteField import FiniteField
from generic.Polynomial import Polynomial
from generic.PolynomialsOverField import PolynomialsOverField
from specific.PolynomialsOverFq import PolynomialsOverFq
from specific.Zp import Zp


class Berlekamp:

    def __init__(self, f):
        # field = Fq, f = f, fx = F[X]
        self.field = f.base_ring
        self.f = f
        self.fx = PolynomialsOverField(f.base_ring)
        self._l = f.degree()

    def compute(self):
        print("Computando base...")
        base = self.step_b1()
        print("Buscando elemento aleatorio")
        return self.step_b2(base)

    # elem is an element in F[X]
    def rep(self, elem):
        return self.fx.rem(elem, self.f)

    # Compute a Berlekamp basis for E := F[X]/(f), where f is a polynomial over F
    def step_b1(self):
        # xi = X * 1_F
        xi = Polynomial([self.field.add_id(), self.field.mul_id()], self.field)
        # l x l matrix with coefficients in F
        q_matrix = []
        # compute alpha = xi^q through repeated squaring.
        # we do this by making all calculations on F[X] and then computing the associated elements in E via rem with f.
        alpha = self.fx.exp(xi, self.field.cardinality())
        alpha = self.rep(alpha)

        print("Tengo el representante")

        # beta = 1_E, no quotient needed in this case.
        beta = self.fx.mul_id()

        for i in range(self._l):
            # Row_i(Q) = Vec_S(beta)
            q_matrix.append(self.vec_s(copy.deepcopy(beta.coefs)))
            # Q[i][i] = Q[i][i] - 1
            q_matrix[i][i] = self.field.add(q_matrix[i][i], self.field.add_inv(self.field.mul_id()))
            # beta = beta * alpha
            beta = self.rep(self.fx.mul(beta, alpha))
        print("Q = ", q_matrix)

        # compute a basis of the row null space of Q
        _, x_matrix = self.gaussian_elimination(q_matrix)
        res = [Polynomial(elem, self.field) for elem in x_matrix]
        print("base = ", res)
        return res

    def step_b2(self, base):
        h_list = [self.f]
        base_list = [self.rep(elem) for elem in base]
        m1_computed = self.m1()
        h_p = []

        while len(h_list) < len(base_list):
            g_aux = [self.fx.mul_field_scalar(gi, self.field.random_elem()) for gi in base_list]
            # g = c1 * g1 + ... + cr * gr (in F[X])
            g = reduce(self.fx.add, g_aux)
            h_p = []
            for h in h_list:
                # beta = [g]_h (in F[X]/(h))
                beta = self.fx.rem(g, h)
                # compute rep(M1(beta))
                eval_m1_beta = self.rep(self.fx.evaluate_polynomial(m1_computed, beta))
                d = self.fx.gcd(eval_m1_beta, h)
                if d == self.fx.mul_id() or d == h or d.degree() < 1:
                    h_p.append(h)
                else:
                    h_p.append(d)
                    h_p.append(self.fx.div(h, d))
            h_list = h_p
        return h_p

    def m1(self):
        coefs = []
        if self.field.p > 2:
            coefs = [self.field.add_id()] * ((self.field.cardinality() - 1) // 2 + 1)
            coefs[(self.field.cardinality() - 1) // 2] = self.field.mul_id()
            coefs[0] = self.field.add_inv(self.field.add_id())
        else:
            # 1 + X^2 + X^4 + ...
            coefs = [self.field.mul_id() if i % 2 == 0 else self.field.add_id() for i in range(2 * self.field.k - 1)]
        return Polynomial(coefs, self.field)

    # coordinates of beta in standard basis S = {1, X, X^2,..., X^l}
    def vec_s(self, beta):
        res = [self.field.add_id()] * self._l
        res[:len(beta)] = beta
        return res

    # returns (gaussian elimination of matrix, basis for the row null space of matrix)
    def gaussian_elimination(self, _matrix):
        q_matrix = copy.deepcopy(_matrix)
        f_one = self.field.mul_id()
        f_zero = self.field.add_id()

        # identity matrix
        x_matrix = [[copy.deepcopy(f_one) if row == col else copy.deepcopy(f_zero) for col in range(len(q_matrix[0]))] for row in range(len(q_matrix))]

        r = -1
        for j in range(len(q_matrix[0])):
            print("Q matrix: ", q_matrix)
            print("X matrix: ", x_matrix)
            _l = -1
            i = r
            while _l == -1 and i < len(q_matrix)-1:
                i += 1
                if q_matrix[i][j] != f_zero:
                    _l = i
            if _l != -1:
                r += 1

                # swap rows r and l of X
                aux = x_matrix[r]
                x_matrix[r] = x_matrix[_l]
                x_matrix[_l] = aux
                # swap rows r and l of Q
                aux = q_matrix[r]
                q_matrix[r] = q_matrix[_l]
                q_matrix[_l] = aux

                # B(r, j)^(-1)
                row_divisor = self.field.mul_inv(q_matrix[r][j])

                # Row_r(X) = B(r, j)^(-1) * Row_r(X)
                x_matrix[r] = [self.field.mul(row_divisor, x) for x in x_matrix[r]]
                # Row_r(Q) = B(r, j)^(-1) * Row_r(Q)
                q_matrix[r] = [self.field.mul(row_divisor, x) for x in q_matrix[r]]

                for _i in range(len(q_matrix)):
                    if _i != r:
                        # Row_i(X) = Row_i(X) - B(i, j) * Row_r(X)
                        aux = [self.field.add_inv(self.field.mul(q_matrix[_i][j], x)) for x in x_matrix[r]]
                        x_matrix[_i] = [self.field.add(x, y) for x, y in zip(x_matrix[_i], aux)]
                        # Row_i(B) = Row_i(B) - B(i, j) * Row_r(B)
                        aux = [self.field.add_inv(self.field.mul(q_matrix[_i][j], x)) for x in q_matrix[r]]
                        q_matrix[_i] = [self.field.add(x, y) for x, y in zip(q_matrix[_i], aux)]
        r += 1
        return q_matrix, x_matrix[-(len(q_matrix)-r):]


if __name__ == "__main__":

    # Z3 = Zp(3)
    # matrix = [[0, 1, 1], [2, 1, 2], [2, 2, 0]]
    # elimination, change = Berlekamp(Polynomial([1], Z3)).gaussian_elimination(matrix)
    # print(".............")
    # print(change)
    # print(elimination)
    #
    # F4 = FiniteField(2, 2, Polynomial([1, 1, 1], Zp(2)))
    # F4One = F4.mul_id()
    # F4Zero = F4.add_id()
    # F4Alpha = Polynomial([0, 1], Zp(2))
    # F4AlphaPlusOne = Polynomial([1, 1], Zp(2))
    # f = Polynomial([F4Alpha, F4AlphaPlusOne, F4One], F4)
    # print(".............")
    # print(Berlekamp(f).step_b1())
    #
    # print(".............")
    # print("Compute factorization through Berkelamp's algorithm")
    # print(Berlekamp(f).step_b2([Polynomial([F4One], F4), Polynomial([F4Zero, F4One], F4)]))
    # print(Berlekamp(f).compute())
    # factors = Berlekamp(f).compute()
    # for factor in factors:
    #     print("Factor", PolynomialsOverField(F4).obtain_monic_representation(factor))
    #
    # Z7 = Zp(7)
    # g = Polynomial([6,11,6,1], Z7)
    # factors = Berlekamp(g).compute()
    # for factor in factors:
    #     print("Factor", PolynomialsOverField(Z7).obtain_monic_representation(factor))

    # Z6473 = Zp(1201)
    # g = Polynomial([-40, -62, -21, 2,1], Z6473)
    # factors = Berlekamp(g).compute()
    # for factor in factors:
    #     print("Factor", PolynomialsOverField(Z6473).obtain_monic_representation(factor))

    Z3 = Zp(3)
    g = Polynomial([1, 1, 2, 3, 1], Z3)
    factors = Berlekamp(g).compute()
    for factor in factors:
        print("Factor", PolynomialsOverField(Z3).obtain_monic_representation(factor))
