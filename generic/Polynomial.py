from generic.Ring import Ring


class Polynomial:

    def __init__(self, coefs, base_ring):
        self.coefs = coefs
        self.base_ring = base_ring

    def degree(self):
        return len(self.coefs) - 1

    def get_coefs(self):
        return self.coefs

    def get_base_ring(self):
        return self.base_ring

    def get_leading_coef(self):
        return self.coefs[self.size()]
