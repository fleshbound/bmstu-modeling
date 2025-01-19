import math
from math import erf, sqrt

from numpy import log as ln
from numpy.random import normal
from random import random


class EvenDistribution:
    def __init__(self, a: float, b: float):
        self.a = a
        self.b = b

    def generate(self):
        return self.a + (self.b - self.a) * random()


class PoissonDistribution:
    """
    Класс, представляющий распределение Пуассона.

    Attributes:
        lambd (float): Параметр lambda (среднее значение) распределения.

    Methods:
        generate(): Генерирует случайное значение Пуассоновской случайной величины.
    """

    def __init__(self, lambd):
        """
        Конструктор класса PoissonDistribution.

        Args:
            lambd (float): Параметр lambda (среднее значение) распределения.
        """
        if lambd <= 0:
            raise ValueError("Lambda must be a positive value.")
        self.lambd = lambd

    def generate(self):
        """
        Генерирует случайное значение Пуассоновской случайной величины.
        Использует алгоритм, основанный на равномерном распределении.

        Returns:
            int: Случайное значение Пуассоновской случайной величины.
        """
        L = math.exp(-self.lambd)
        k = 0
        p = 1
        while p > L:
            k += 1
            p *= random()
        return k - 1


class ErlangDistribution:
    def __init__(self, k: int, l: float):
        self.k = k
        self.l = l

    def generate(self):
        return 1 / (self.k * self.l) * \
            sum([ln(1 - random()) for _ in range(self.k)])


class NormalDistribution:
    def __init__(self, m: float, sigma: float):
        self.sigma = sigma
        self.m = m

    def generate(self):
        return normal(self.m, self.sigma)
        #return self.m + self.sigma * random()