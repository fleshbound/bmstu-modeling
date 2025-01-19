from scipy.stats import poisson
import numpy as np

from lab_01.src.draw import draw_universal, draw_poisson
from lab_01.src.values import ud_function, ud_density


def main():
    run = True
    choice = 1
    while run:
        choice = input('Вариант распределения (1 - равномерное, 2 - Пуассона): ')

        try:
            choice = int(choice)
            run = False
        except Exception:
            run = True

    if choice == 1:
        a = float(input("Введите a: "))
        b = float(input("Введите b: "))
        # a = 0
        # b = 5
        delta = b - a
        x = np.arange(a - delta / 2, b + delta / 2, 0.001)
        y_function = [ud_function(a, b, _x) for _x in x]
        y_density = [ud_density(a, b, _x) for _x in x]
        draw_universal(x, y_function, y_density, 'Равномерное распределение')

    # mu = 5
    # x = np.linspace(poisson.ppf(0.01, mu), poisson.ppf(0.95, mu))
    # y_function = poisson_function(x, mu)
    # y_density = poisson_density(x, mu)
    # lin = np.arange(poisson.ppf(0.0001, mu), poisson.ppf(0.99, mu), 0.01)
    elif choice == 2:
        mu = float(input("Введите lambda: "))
        lin = np.arange(0,18, 0.01)
        y_function = poisson.cdf(lin, mu)
        y_density = poisson.pmf(lin, mu)
        draw_poisson(lin, y_function, y_density, 'Распределение Пуассона')


if __name__ == '__main__':
    main()
