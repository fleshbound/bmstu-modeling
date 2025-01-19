import math
import random

from abc import ABC, abstractmethod

class Distribution(ABC):
    
    @abstractmethod
    def __init__(self) -> None:
        raise NotImplementedError("Метод должен быть реализован в классах-потомках!")
    

    @abstractmethod
    def get_value(self) -> float:
        raise NotImplementedError("Метод должен быть реализован в классах-потомках!")

class PoissonDistribution(Distribution):
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

    def get_value(self):
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
            p *= random.random()
        return k - 1


class UniformDistribution(Distribution):

    def __init__(self, lower: float, upper: float) -> None:
        self.lower = lower
        self.upper = upper
    
    
    def get_value(self) -> float:
        return random.uniform(self.lower, self.upper)
    
    
class ErlangDistribution(Distribution):
    
    def __init__(self, shape: int, rate: float) -> None:
        if not isinstance(shape, int):
            raise ValueError("Параметр 'порядок' распределения Эрланга должен быть целым числом!")
        if shape < 1:
            raise ValueError("Параметр 'порядок' распределения Эрланга должен быть больше 0!")
        
        if rate <= 0.00:
            raise ValueError("Параметр 'интенсивность' распределения Эрланга должен быть больше 0!")
        
        self.rate = rate
        self.shape = shape


    def get_value(self) -> float:
        return random.gammavariate(self.shape, self.rate)   