from generic.Ring import Ring
from abc import abstractmethod


class UFD(Ring):
    @abstractmethod
    def gcd(self, elem1, elem2):
        pass
