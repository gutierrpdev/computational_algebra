from abc import ABC, abstractmethod

from generic.EuclideanDomain import EuclideanDomain


# Would rather implement mul_inv and div_mod interdependently and let the user redefine one of the methods
class Field(EuclideanDomain, ABC):
    @abstractmethod
    def mul_inv(self, elem):
        return

    def div_mod(self, elem1, elem2):
        return self.mul(elem1, self.mul_inv(elem2)), self.add_id()
