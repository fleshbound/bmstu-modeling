from lsm import LeastSquaresMethodLine

import numpy as np
from matplotlib import pyplot as plt

r = 0.35
t_w = 2000
t_0 = 1e4
p = 4

c = 3e10

t_k = [2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
k_1 = [8.2E-3, 2.768E-02, 6.560E-02, 1.281E-01, 2.214E-01, 3.516E-01, 5.248E-01, 7.472E-01, 1.025E+00]
k_2 = [1.600E+00, 5.400E+00, 1.280E+01, 2.500E+01, 4.320E+01, 6.860E+01, 1.024E+02, 1.458E+02, 2.000E+02]

a1, b1 = LeastSquaresMethodLine(x=np.log(t_k), y=np.log(k_1)).get_solve()
a2, b2 = LeastSquaresMethodLine(x=np.log(t_k), y=np.log(k_2)).get_solve()

is_k1 = None

EPS = 1e-7


def t(z: int | float | np.ndarray):
    return (t_w - t_0) * z ** p + t_0


def u_p(z: int | float | np.ndarray):
    return 3.084e-4 / (np.exp(4.799e4 / t(z)) - 1)


def k1(_t: int | float | np.ndarray):
    return np.exp(a1 * np.log(t(_t)) + b1)


def k2(_t: int | float | np.ndarray):
    return np.exp(a2 * np.log(t(_t)) + b2)


def der_u(z, _f):
    if is_k1:
        return -3 * r * k1(t(z)) / c * _f

    return -3 * r * k2(t(z)) / c * _f


def der_f(z, u, _f):
    if is_k1:
        common_part = r * c * k1(t(z)) * (u_p(z) - u)
    else:
        common_part = r * c * k2(t(z)) * (u_p(z) - u)

    res = -_f / z + common_part if abs(z) > EPS else common_part / 2

    return res


def div_flux(z, u):
    if is_k1:
        return c * k1(z) * (u_p(z) - u)

    return c * k2(z) * (u_p(z) - u)


def _lambda(z):
    if is_k1:
        return c / (3 * k1(z))

    return c / (3 * k2(z))


def _lambda_new(z):
    if is_k1:
        return 1 / (3 * k1(z))

    return 1 / (3 * k2(z))


def _p(z):
    if is_k1:
        return c * k1(z)

    return c * k2(z)


def f(z):
    if is_k1:
        return c * k1(z) * u_p(z)

    return c * k2(z) * u_p(z)


def _p_new(z):
    if is_k1:
        return 1 * k1(z)

    return 1 * k2(z)


def f_new(z):
    if is_k1:
        return 1 * k1(z) * u_p(z)

    return 1 * k2(z) * u_p(z)


def v_n(z, h):
    return ((z + h / 2) ** 2 - (z - h / 2) ** 2) / 2


def kappa1(z1, z2):
    return (_lambda(z1) + _lambda(z2)) / 2


def kappa2(z1, z2):
    return (2 * _lambda(z1) * _lambda(z2)) / (_lambda(z1) + _lambda(z2))


def kappa_new(z1, z2):
    return (_lambda_new(z1) + _lambda_new(z2)) / 2


def left_boundary_condition(z0, g0, h):
    m0 = +kappa_new(z0, z0 + h) * (z0 + h / 2) / (r ** 2 * h) + _p_new(z0 + h / 2) * (z0 + h / 2) * h / 8  # правильно
    k0 = -kappa_new(z0, z0 + h) * (z0 + h / 2) / (r ** 2 * h) - _p_new(z0 + h / 2) * (z0 + h / 2) * h / 8
    p0 = +(f_new(z0 + h / 2) * (z0 + h / 2) + f_new(z0) * z0) * h / 4

    return k0, m0, p0


# При x = N
def right_boundary_condition(zn, h):
    mn = -kappa_new(zn - h, zn) * (zn - h / 2) / (r ** 2 * h) + _p_new(zn - h / 2) * (zn - h / 2) * h / 8  # правильно
    kn = +kappa_new(zn - h, zn) * (zn - h / 2) / (r ** 2 * h) + 0.393 * 1 * zn / r + _p_new(zn) * zn * h / 4 + _p_new(zn - h / 2) * (zn - h / 2) * h / 8
    pn = +(f_new(zn) * zn + f_new(zn - h / 2) * (zn - h / 2)) * h / 4

    return kn, mn, pn


def right_sweep(a, b, h):
    k0, m0, p0 = left_boundary_condition(a, 0, h)
    kn, mn, pn = right_boundary_condition(b, h)

    ksi = [0, -k0 / m0]
    eta = [0, p0 / m0]

    z = h
    n = 1

    while z < b + h / 2:
        a_n = (z - h / 2) * (kappa_new(z - h, z)) / (r ** 2 * h)
        c_n = ((z + h / 2) * kappa_new(z, z + h)) / (r ** 2 * h)
        b_n = a_n + c_n + _p_new(z) * v_n(z, h)
        d_n = f_new(z) * v_n(z, h)

        ksi.append(c_n / (b_n - a_n * ksi[n]))
        eta.append((a_n * eta[n] + d_n) / (b_n - a_n * ksi[n]))

        n += 1
        z += h

    u = [0] * n
    u[n - 1] = (pn - mn * eta[n - 1]) / (mn * ksi[n - 1] + kn)

    for i in range(n - 2, -1, -1):
        u[i] = ksi[i + 1] * u[i + 1] + eta[i + 1]

    return u


def left_sweep(a, b, h):
    k0, m0, p0 = left_boundary_condition(a, 0, h)
    kn, mn, pn = right_boundary_condition(b, h)

    ksi = [-mn / kn, 0]
    eta = [pn / kn, 0]

    z = b - h
    n = -2
    cnt = 1

    while z > a - h / 2:
        a_n = (z - h / 2) * (kappa_new(z - h, z)) / (r ** 2 * h)
        c_n = ((z + h / 2) * kappa_new(z, z + h)) / (r ** 2 * h)
        b_n = a_n + c_n + _p_new(z) * v_n(z, h)
        d_n = f_new(z) * v_n(z, h)

        # print(a_n, b_n, c_n, d_n)
        ksi.insert(0, a_n / (b_n - c_n * ksi[n]))
        eta.insert(0, (c_n * eta[n] + d_n) / (b_n - c_n * ksi[n]))

        n -= 1
        z -= h
        cnt += 1

    u = [0] * (cnt)
    u[0] = (p0 - k0 * eta[0]) / (m0 + k0 * ksi[0])

    for i in range(1, cnt):
        u[i] = ksi[i - 1] * u[i - 1] + eta[i - 1]

    return u


def meetings_sweep(a, b, h, n_eq):
    k0, m0, p0 = left_boundary_condition(a, 0, h)
    kn, mn, pn = right_boundary_condition(b, h)

    z = h
    n = 1

    ksi_r = [0, -k0 / m0]
    eta_r = [0, p0 / m0]

    while z < n_eq * h:
        a_n = (z - h / 2) * (kappa_new(z - h, z)) / (r ** 2 * h)
        c_n = ((z + h / 2) * kappa_new(z, z + h)) / (r ** 2 * h)
        b_n = a_n + c_n + _p_new(z) * v_n(z, h)
        d_n = f_new(z) * v_n(z, h)

        ksi_r.append(c_n / (b_n - a_n * ksi_r[n]))
        eta_r.append((a_n * eta_r[n] + d_n) / (b_n - a_n * ksi_r[n]))

        n += 1
        z += h

    # левая часть прогонки
    ksi_l = [-mn / kn, 0]
    eta_l = [pn / kn, 0]

    z = b - h
    n1 = -2
    cnt = 1

    while z > n_eq * h - h / 2:
        a_n = (z - h / 2) * (kappa_new(z - h, z)) / (r ** 2 * h)
        c_n = ((z + h / 2) * kappa_new(z, z + h)) / (r ** 2 * h)
        b_n = a_n + c_n + _p_new(z) * v_n(z, h)
        d_n = f_new(z) * v_n(z, h)

        ksi_l.insert(0, a_n / (b_n - c_n * ksi_l[n1]))
        eta_l.insert(0, (c_n * eta_l[n1] + d_n) / (b_n - c_n * ksi_l[n1]))

        n1 -= 1
        z -= h
        cnt += 1

    u = [0] * (n + cnt)
    u[n_eq] = (ksi_r[-1] * eta_l[0] + eta_r[-1]) / (1 - ksi_r[-1] * ksi_l[0])

    for i in range(n_eq - 1, -1, -1):
        u[i] = ksi_r[i + 1] * u[i + 1] + eta_r[i + 1]

    for i in range(n_eq + 1, n + cnt):
        _i = i - n_eq
        u[i] = ksi_l[_i - 1] * u[i - 1] + eta_l[_i - 1]

    return u


def flux1(u, z, h):
    f_res = [0]

    try:
        for i in range(1, len(u) - 1):
            curr_f = -(_lambda(z[i]) / r) * (u[i + 1] - u[i - 1]) / (2 * h)
            f_res.append(curr_f)

        f_res.append(-(_lambda(z[len(u) - 1]) / r) * (3 * u[-1] - 4 * u[-2] + u[-3]) / (2 * h))
    except IndexError:
        print("i", end="")
        return []

    return f_res


def flux2(z, u, h):
    _f = [0]
    f_res = [0]

    if is_k1:
        k = k1
    else:
        k = k2

    try:
        for i in range(1, len(z)):
            _f.append(k(z[i]) * (u_p(z[i]) - u[i]) * z[i])
            f_res.append((c * r / z[i]) * h * ((_f[0] + _f[i]) / 2 + sum(_f[1:-1])))
    except IndexError:
        print("i", end="")
        return []

    return f_res


def flux3(z, u, h):
    _f = [0]
    f_res = [0]

    if is_k1:
        k = k1
    else:
        k = k2

    for i in range(1, len(z)):
        _f.append(k(z[i]) * (u_p(z[i]) - u[i]) * z[i])

        _sum = 0

        for _k in range(1, len(z[:i]), 2):
            _sum += (_f[_k - 1] + 4 * _f[_k] + _f[_k + 1])

        f_res.append((c * r / z[i]) * (h / 3) * _sum)

    return f_res


def write_result_to_file(filepath, z_res, u_res, f_res, f_res2):
    file = open(filepath, "w", encoding="utf-8")

    file.write(f"Число узлов n = {len(z_res)}\n")
    file.write("-" * 86 + "\n")
    file.write(f'| {"x": ^7} | {"u(z)": ^22} | {"f(z)": ^22} | {"f(z) integral": ^22} |\n')
    file.write("-" * 86 + "\n")

    for i in range(len(z_res)):
        file.write(f"| {z_res[i]: ^7.5f} | {u_res[i]: ^22.10e} | {f_res[i]: ^22.10e} | {f_res2[i]: ^22.10e} |\n")

    file.write("-" * 86)

    file.close()


def cmp_res_by_input_data(a, b, h):
    plt.ioff()

    z_res = np.arange(a, b + h, h)
    _t_0 = np.linspace(100, 1e6, 100)
    u0_t_0_res = []
    f0_t_0_res = []
    global t_0
    for curr_t_0 in _t_0:
        t_0 = curr_t_0
        u_res = right_sweep(a, b, h)
        f_res = flux1(u_res, z_res, h)
        u0_t_0_res.append(u_res[0])
        f0_t_0_res.append(f_res[0])
    t_0 = 1e4
    fig, axs = plt.subplots(2, 3, figsize=(13, 9))
    axs[0, 0].plot(_t_0, u0_t_0_res, 'r')
    axs[0, 0].set_xlabel("T_0")
    axs[0, 0].set_ylabel("u(0)")
    axs[0, 0].grid()
    axs[1, 0].plot(_t_0, f0_t_0_res, 'r')
    axs[1, 0].set_xlabel("T_0")
    axs[1, 0].set_ylabel("F(0)")
    axs[1, 0].grid()

    p_list = list(range(1, 11))
    u0_p_res = []
    f0_p_res = []
    global p
    for curr_p in p_list:
        p = curr_p
        u_res = right_sweep(a, b, h)
        f_res = flux1(u_res, z_res, h)
        u0_p_res.append(u_res[0])
        f0_p_res.append(f_res[0])
    p = 4
    axs[0, 1].plot(p_list, u0_p_res, 'r')
    axs[0, 1].set_xlabel("p")
    axs[0, 1].set_ylabel("u(0)")
    axs[0, 1].grid()
    axs[1, 1].plot(p_list, f0_p_res, 'r')
    axs[1, 1].set_xlabel("p")
    axs[1, 1].set_ylabel("F(0)")
    axs[1, 1].grid()

    r_list = np.arange(0.1, 1 + 0.05, 0.05)
    u0_r_res = []
    f0_r_res = []
    global r
    for curr_r in r_list:
        r = curr_r
        u_res = right_sweep(a, b, h)
        f_res = flux1(u_res, z_res, h)
        u0_r_res.append(u_res[0])
        f0_r_res.append(f_res[0])
    r = 0.35
    axs[0, 2].plot(r_list, u0_r_res, 'r')
    axs[0, 2].set_xlabel("r")
    axs[0, 2].set_ylabel("u(0)")
    axs[0, 2].grid()

    axs[1, 2].plot(r_list, f0_r_res, 'r')
    axs[1, 2].set_xlabel("r")
    axs[1, 2].set_ylabel("u(0)")
    axs[1, 2].grid()

    # Сохраняем графики в файл
    fig.savefig(f"data/cmp_res_by_input_data.png")

    # Включаем интерактивный режим обратно
    plt.ion()


def process(start, stop, step, a, b, sweeps):
    print(f"start data rewriting: z in [{a}, {b}], n in [{start}, {stop}], step = {step}, sweeps = {sweeps}")

    for n in range(start, stop + 1, step):
        print(f".", end="")
        h = (b - a) / n

        for sweep in sweeps:
            if sweep == 'r':
                u_res = right_sweep(a, b, h)
            elif sweep == 'l':
                u_res = left_sweep(a, b, h)
            else:
                u_res = meetings_sweep(a, b, h, n // 2)

            # print(u_res)

            z_res = np.arange(a, b + h, h)
            f_res = flux1(u_res, z_res, h)
            f_res2 = flux2(z_res, u_res, h)

            if len(f_res) == 0:
                print("e1", end="")
                continue

            if len(f_res2) == 0:
                print("e2", end="")
                continue

            # f_res2 = flux3(z_res, u_res, h)
            up_res = [0] * len(z_res)
            div_f = [0] * len(z_res)

            for i in range(len(z_res)):
                up_res[i] = u_p(z_res[i])
                div_f[i] = div_flux(z_res[i], u_res[i])

            write_result_to_file(f"data/{sweep}/{n}.txt", z_res, u_res, f_res, f_res2)

    print("\nsuccess")


def get_data():
    global is_k1
    is_k1 = True
    a, b = 0, 1
    sweeps = ['r', 'l', 'm']
    start, stop, step = 10, 100, 10
    process(start, stop, step, a, b, sweeps)
    start, stop, step = 200, 1000, 100
    process(start, stop, step, a, b, sweeps)
    start, stop, step = 2000, 10000, 1000
    process(start, stop, step, a, b, sweeps)


def main() -> None:
    global is_k1
    var = int(input('Какой вариант данных использовать? (1, 2): ').strip())
    if var == 1:
        is_k1 = True
    elif var == 2:
        is_k1 = False
    else:
        print("Неизвестный вариант")
        exit(1)

    a, b = 0, 1
    sweep = input('Какую прогонку использовать? (r, l, m): ').strip()
    n_inp = input('Введите число узлов (по умолчанию 200): ').strip()
    if n_inp == '':
        n = 200
    else:
        n = int(n_inp)
    h = (b - a) / n

    if sweep == 'r':
        u_res = right_sweep(a, b, h)
    elif sweep == 'l':
        u_res = left_sweep(a, b, h)
    elif sweep == 'm':
        u_res = meetings_sweep(a, b, h, n // 2)
    else:
        print('Неизвестный способ прогонки')
        exit(1)

    z_res = np.arange(a, b + h, h)
    f_res = flux1(u_res, z_res, h)
    f_res2 = flux2(z_res, u_res, h)
    # f_res2 = flux3(z_res, u_res, h)

    up_res = [0] * len(z_res)
    div_f = [0] * len(z_res)

    for i in range(len(z_res)):
        up_res[i] = u_p(z_res[i])
        div_f[i] = div_flux(z_res[i], u_res[i])

    write_result_to_file(f"data/{sweep}/{n}.txt", z_res, u_res, f_res, f_res2)

    plt.figure(figsize=(9, 6))
    plt.subplot(2, 2, 1)
    plt.plot(z_res, u_res, 'r', label='u(z)')
    plt.plot(z_res, up_res, 'b', label='u_p')
    plt.legend()
    plt.grid()

    plt.subplot(2, 2, 2)
    plt.plot(z_res, f_res, 'g', label='F(z)')
    plt.legend()
    plt.grid()

    plt.subplot(2, 2, 3)
    plt.plot(z_res, f_res2, 'g', label='F(z) integral')
    plt.legend()
    plt.grid()

    plt.subplot(2, 2, 4)
    plt.plot(z_res, div_f, 'y', label='divF(z)')
    plt.legend()
    plt.grid()

    plt.show()

    # cmp_res_by_input_data(a, b, h)


if __name__ == '__main__':
    main()
    # get_data()
