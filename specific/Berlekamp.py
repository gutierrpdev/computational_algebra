import copy

from generic.Polynomial import Polynomial
from specific.PolynomialsOverFq import PolynomialsOverFq
from specific.Zp import Zp


class Berlekamp:

    def __init__(self, f):
        # field = Fq, f = f, fx = F[X]
        self.field = f.base_ring
        self.f = f
        self.fx = PolynomialsOverFq(f.base_ring)
        self._l = f.degree()

    def compute(self):
        self.step_b1()

    # Compute a Berlekamp basis for E := F[X]/(f), where f is a polynomial over F
    def step_b1(self):
        # xi = X * 1_F
        xi = Polynomial([self.field.add_id(), self.field.mul_id()], self.field)
        # l x l matrix with coefficients in F
        q_matrix = []
        # compute alpha = xi^q through repeated squaring.
        # we do this by making all calculations on F[X] and then computing the associated elements in E via rem with f.
        alpha = self.fx.exp(xi, self.field.cardinality())
        alpha = self.fx.rem(alpha, self.f)

        # beta = 1_E, no quotient needed in this case.
        beta = self.fx.mul_id()

        for i in range(self.l):
            # Row_i(Q) = Vec_S(beta)
            q_matrix.append(self.vec_s(copy.deepcopy(beta.coefs)))
            # Q[i][i] = Q[i][i] - 1
            q_matrix[i][i] = self.field.add(q_matrix[i][i], self.field.add_inv(self.field.mul_id()))
            # beta = beta * alpha
            beta = self.fx.rem(self.fx.mul(beta, alpha), self.f)

        # compute a basis of of the row null space of Q
        self.gaussian_elimination(q_matrix)

    # coordinates of beta in standard basis S = {1, X, X^2,..., X^l}
    def vec_s(self, beta):
        res = [self.field.add_id()] * self.l
        res[:len(beta)] = beta
        return res

    # returns gaussian elimination of matrix, matrix whose m-n last rows form a basis for the row null space of matrix
    def gaussian_elimination(self, _matrix):
        q_matrix = copy.deepcopy(_matrix)
        f_one = self.field.mul_id()
        f_zero = self.field.add_id()

        # identity matrix
        x_matrix = [[f_one if row == col else f_zero for col in range(len(q_matrix[0]))] for row in range(len(q_matrix))]

        r = -1
        for j in range(len(q_matrix[0])):
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
        return q_matrix, x_matrix


if __name__ == "__main__":
    Z3 = Zp(3)
    matrix = [[0, 1, 1], [2, 1, 2], [2, 2, 0]]
    elimination, change = Berlekamp(Polynomial([1], Z3)).gaussian_elimination(matrix)
    print(".............")
    print(change)
    print(elimination)
