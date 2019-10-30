from generic.Field import Field


class FiniteField(Field):

    def __init__(self, p, k, f):
        self.p, self.k, self.f = p, k, f
