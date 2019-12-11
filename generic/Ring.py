import copy
from abc import ABC, abstractmethod


class Ring(ABC):
    @abstractmethod
    def add_id(self):
        pass

    @abstractmethod
    def mul_id(self):
        pass

    @abstractmethod
    def add_inv(self, elem):
        pass

    @abstractmethod
    def add(self, elem1, elem2):
        pass

    @abstractmethod
    def mul(self, elem1, elem2):
        pass

    def mul_scalar(self, elem, scalar):
        res = self.add_id()
        for i in range(abs(scalar)):
            res = self.add(res, elem)
        if scalar < 0:
            res = self.add_inv(res)
        return res

    # compute elem**scalar (scalar must be >= 0)
    def exp(self, elem, scalar):
        res = copy.deepcopy(elem)
        if scalar == 0:
            return self.mul_id()
        elif scalar > 0:
            for i in range(scalar-1):
                res = self.mul(res, elem)
        return res
