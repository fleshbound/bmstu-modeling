import numpy as np
from collections import Counter

def get_probability(num):
    # Проверяем, что число является целым и положительным
    #if not isinstance(num, int) or num <= 0:
    #    raise ValueError("Введите положительное целое число.")
    
    # Определяем количество цифр в числе
    num_digits = len(str(num))
    
    # Определяем диапазон в зависимости от количества цифр
    if num_digits == 1:
        lower_bound = 0
        upper_bound = 9
    elif num_digits == 2:
        lower_bound = 10
        upper_bound = 99
    elif num_digits == 3:
        lower_bound = 100
        upper_bound = 999
    elif num_digits == 4:
        lower_bound = 1000
        upper_bound = 9999
    else:
        raise ValueError("Функция поддерживает только числа с 1-4 цифрами.")
    
    # Общее количество чисел в диапазоне
    total_numbers = upper_bound - lower_bound + 1
    
    # Вероятность появления данного числа
    probability = 1 / total_numbers
    
    return probability

def shannon_entropy(data):
    counter = Counter(data)
    total_count = len(data)
    
    probabilities = [count / total_count for count in counter.values()]

    #probabilities = [get_probability(num) for num in data]
    
    entropy = -sum(p * np.log2(p) for p in probabilities if p > 0)
    return entropy
