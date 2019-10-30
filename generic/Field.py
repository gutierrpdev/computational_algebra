from abc import ABC, abstractmethod

from generic.EuclideanDomain import EuclideanDomain


class Field(EuclideanDomain, ABC):
    @abstractmethod
    def mul_inv(self, elem):
        pass
