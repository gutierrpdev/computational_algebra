class Polynomial:

    def __init__(self, coefs, base_ring):
        i = len(coefs) - 1
        while i > 0 and coefs[i] == base_ring.add_id():
            i -= 1
        self.coefs = coefs[:i+1]
        self.base_ring = base_ring

    def degree(self):
        return len(self.coefs) - 1

    def get_coefs(self):
        return self.coefs

    def get_base_ring(self):
        return self.base_ring

    def get_leading_coef(self):
        return self.coefs[self.degree()]

    def __eq__(self, other):
        if isinstance(other, Polynomial):
            return self.coefs == other.coefs and (type(self.base_ring) is type(other.base_ring))
        return NotImplemented
