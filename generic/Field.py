from abc import ABC, abstractmethod

from generic.EuclideanDomain import EuclideanDomain


class Field(EuclideanDomain, ABC):
    @abstractmethod
    def mul_inv(self, elem):
        pass

    def div_mod(self, elem1, elem2):
        return self.mul(elem1, self.mul_inv(elem2)), self.add_id()
