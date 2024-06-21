import numpy as np
from matplotlib import pyplot as plt


def task1():
    min_x = 0
    max_x = 1.2
    n_x = 50
    x = np.linspace(min_x, max_x, n_x)
    h = (max_x - min_x) / n_x
    picard_y = np.empty_like(x)
    exp_y = np.empty_like(x)
    euler_y = np.empty_like(x)
    euler_y1 = np.empty_like(x)
    euler_y[0] = 1
    euler_y1[0] = 2

    for i in range(len(x)):
        exp_y[i] = (1 + 2 * x[i]
                    - 0.7 * (x[i] ** 2)
                    - 0.77 / 3 * (x[i] ** 3)
                    + 0.051 * (x[i] ** 4))

        picard_y[i] = (1 + 2 * x[i]
                       - 0.7 * (x[i] ** 2)
                       - 0.77 / 3 * (x[i] ** 3)
                       - 47 / 1500 * (x[i] ** 4)
                       - 0.0007 * (x[i] ** 5)
                       - 1 / (1.2 * (10 ** 5)) * (x[i] ** 6))

        if i > 0:
            euler_y1[i] = euler_y1[i - 1] + h * (-0.1 * (euler_y1[i - 1] ** 2) - (1 + 0.1 * x[i - 1]) * euler_y[i - 1])
            euler_y[i] = euler_y[i - 1] + h * euler_y1[i - 1]

    plt.plot(x, exp_y, 'o', label='Разложение в ряд')
    plt.plot(x, picard_y, 'x', label='Метод Пикара')
    plt.plot(x, euler_y, '*', label='Метод Эйлера')
    plt.legend(fontsize=10)
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def task2():
    min_y = -1.5
    max_y = 1.5
    n_y = 90
    y = np.linspace(min_y, max_y, n_y)
    picard_x1 = np.empty_like(y)
    picard_x2 = np.empty_like(y)
    picard_x3 = np.empty_like(y)
    picard_x4 = np.empty_like(y)
    x = np.empty_like(y)

    for i in range(len(x)):
        x[i] = - (y[i] ** 2) / 2 - 0.5 + np.exp(y[i] ** 2)
        picard_x1[i] = 0.5 + (y[i] ** 4) / 4 + (y[i] ** 2) / 2
        picard_x2[i] = 0.5 + (y[i] ** 4) / 2 + (y[i] ** 2) / 2 + (y[i] ** 6) / 12
        picard_x3[i] = 0.5 + (y[i] ** 4) / 2 + (y[i] ** 2) / 2 + (y[i] ** 6) / 6 + (y[i] ** 8) / 48
        picard_x4[i] = 0.5 + (y[i] ** 4) / 2 + (y[i] ** 2) / 2 + (y[i] ** 6) / 6 + (y[i] ** 8) / 24 + (y[i] ** 10) / 240

    plt.plot(x, y, 'o', label='Аналитический метод')
    plt.plot(picard_x1, y, 'x', label='Метод Пикара (1-е приближение)')
    plt.plot(picard_x2, y, '*', label='Метод Пикара (2-е приближение)')
    plt.plot(picard_x3, y, '+', label='Метод Пикара (3-е приближение)')
    plt.plot(picard_x4, y, 'v', label='Метод Пикара (4-е приближение)')
    plt.legend(fontsize=10)
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def task3():
    # x_min = 0, x_max ~ 1, h=1e-4 ???
    # overflow: x ~
    # u(2) = 317.25
    # min_x = 0
    # max_x = 2
    # step_x = 0.0000005
    min_x = 0
    max_x = 2
    step_x = 0.0000005
    x = np.arange(min_x, max_x + step_x, step_x)
    picard_y1 = np.empty_like(x)
    picard_y2 = np.empty_like(x)
    picard_y3 = np.empty_like(x)
    picard_y4 = np.empty_like(x)
    euler_y = np.empty_like(x)
    euler_y[0] = 0

    for i in range(len(x)):
        if i > 0:
            euler_y[i] = euler_y[i - 1] + step_x * (x[i - 1] ** 2 + euler_y[i - 1] ** 2)

        picard_y1[i] = (x[i] ** 3) / 3
        picard_y2[i] = (x[i] ** 3) / 3 + (x[i] ** 7) / 63
        picard_y3[i] = (x[i] ** 3) / 3 + (x[i] ** 7) / 63 + (x[i] ** 11) * 2 / 2079 + (x[i] ** 15) / 59535
        picard_y4[i] = ((x[i] ** 3) / 3 + (x[i] ** 7) / 63 + (x[i] ** 11) * 2 / 2079 + (x[i] ** 15) / 59535
                        + (x[i] ** 15) * 4 / 93555 + (x[i] ** 19) * 2 / 3393495 + (x[i] ** 19) * 4 / 2488563
                        + (x[i] ** 23) * 4 / 99411543 + (x[i] ** 23) * 2 / 86266215 + (x[i] ** 27) * 4 / 3341878155
                        + (x[i] ** 31) / 109876902875)

    print(f'xmin: {min_x}\nxmax: {max_x}\nh: {step_x}')
    print(f'|  x   | picard1 | picard2 | picard3 | picard4 |  euler  |')

    for i in range(0, len(x), 100):
        print(f'| {x[i]:.6f} | {picard_y1[i]:.5f} | {picard_y2[i]:.5f} | {picard_y3[i]:.5f} '
              f'| {picard_y4[i]:.5f} | {euler_y[i]:.5f} |')

    plt.plot(x, euler_y, 'o', label='Метод Эйлера')
    plt.plot(x, picard_y1, 'x', label='Метод Пикара (1-е приближение)')
    plt.plot(x, picard_y2, '*', label='Метод Пикара (2-е приближение)')
    plt.plot(x, picard_y3, '+', label='Метод Пикара (3-е приближение)')
    plt.plot(x, picard_y4, 'v', label='Метод Пикара (4-е приближение)')
    plt.legend(fontsize=10)
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def main():
    # task1()
    # task2()
    task3()


if __name__ == '__main__':
    main()
