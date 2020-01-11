from generic.Field import Field


# QQ is implemented as a float here
class QQ(Field):
    def add_id(self):
        return 0

    def mul_id(self):
        return 1

    def add_inv(self, elem):
        return -elem

    def add(self, elem1, elem2):
        return elem1 + elem2

    def mul(self, elem1, elem2):
        return elem1 * elem2

    def div(self, elem1, elem2):
        return elem1 / elem2
