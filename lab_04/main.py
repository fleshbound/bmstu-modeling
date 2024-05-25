import numpy as np

a1 = 0.0134
b1 = 1
c1 = 4.35e-4
m1 = 1

a2 = 2.049
b2 = 0.563e-3
c2 = 0.528e5
m2 = 1

R = 0.5  # см
l = 10  # см
T0 = 300
t_max = 60  # с
F_max = 50  # Вт/см2

k0 = 1.0  # см-1


def c_fun(T):
    """коэффициент теплоемкости, Дж/см3 К"""
    return a2 + b2 * (T ** m2) - c2 / (T ** 2)


def lambda_fun(T):
    """коэффициент теплопроводности, Вт/см К"""
    return a1 * (b1 + c1 * (T ** m1))


def alpha(x):
    return 0.0125 * l / (x + 0.25 * l)


def p_fun(x):
    return 2 / R * alpha(x)


def kappa(T1, T2):
    return (lambda_fun(T1) + lambda_fun(T2)) / 2


def k_fun(T):
    return k0 * ((T / 300) ** 2)


def F0(t):
    return F_max / t_max * t * np.exp(-(t / t_max - 1))


def f_fun(x, t, T):
    """t = time, T = Temperature"""
    return k_fun(T) * F0(t) * np.exp(-k_fun(T) * x) + 2 * T0 / R * alpha(x)


def left_condition(x0, t, h, t_step, y0, y1):
    """t = time, t_step = time_step; y0, y1 - temperatures"""
    c_fun_half = (c_fun(x0) + c_fun(x0 + h)) / 2
    p_fun_half = (p_fun(x0) + p_fun(x0 + h)) / 2
    f_fun_half = (f_fun(x0, t, y0) + f_fun(x0 + h, t, y1)) / 2
    a = c_fun_half * h / 8
    b = kappa(x0, x0 + h) * t_step / h
    c = p_fun_half * h * t_step / 8
    d = c_fun(x0) * h / 4

    k0 = a + b + c + d + h * t_step * p_fun(x0) / 4
    m0 = a - b + c
    p0 = a * (y0 + y1) + t_step * F0(t) + h * t_step * (f_fun_half + f_fun(x0, t, y0)) / 4 + d * y0

    return k0, m0, p0


def right_condition(xN, t, h, t_step, yN_1, yN):
    """t = time, t_step = time_step; y0, y1 - temperatures"""
    c_fun_half = (c_fun(xN - h) + c_fun(xN)) / 2
    p_fun_half = (p_fun(xN - h) + p_fun(xN)) / 2
    f_fun_half = (f_fun(xN - h, t, yN_1) + f_fun(xN, t, yN)) / 2
    a = c_fun_half * h / 8
    b = kappa(xN - h, xN) * t_step / h
    c = p_fun_half * h * t_step / 8
    d = c_fun(xN) * h / 4

    k0 = a - b + c
    m0 = d + a + b + c + h * t_step * p_fun(xN) / 4
    p0 = d * yN + a * (yN_1 + yN) + h * t_step * (f_fun_half + f_fun(xN, t, yN)) / 4 - t_step * F

    return k0, m0, p0


def right_sweep(a, b, h):
    k0, m0, p0 = 
