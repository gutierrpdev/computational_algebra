from generic.Field import Field


class Zp(Field):

    def __init__(self, p):
        self.p = p

    def mul_inv(self, elem):
        pass

    def div_mod(self, elem1, elem2):
        pass

    def add_id(self):
        return 0

    def mul_id(self):
        return 1

    def add_inv(self, elem):
        pass

    def add(self, elem1, elem2):
        return (elem1 + elem2) % self.p

    def mul(self, elem1, elem2):
        pass
