import numpy as np
import matplotlib.pyplot as plt

T_TABLE = np.array([2000, 3000, 4000,
                    5000, 6000, 7000,
                    8000, 9000, 10000])

K_TABLE_1 = np.array([8.200E-03, 2.768E-02, 6.560E-02, 1.281E-01,
                    2.214E-01, 3.516E-01, 5.248E-01, 7.472E-01,
                    1.025E+00])

K_TABLE_2 = np.array([1.600e+00, 5.400e+00, 1.280e+01,
                      2.500e+01, 4.320e+01, 6.860e+01,
                      1.024e+02, 1.458e+02, 2.000e+02])

DI_EPS = 1e-4
R = 0.35
T_W = 2000  # 10000
T_0 = 10000
P = 4
C = 3e+10
min_z = 0
max_z = 1
n = 1000
z, h = np.linspace(min_z, max_z, n, retstep=True)


def T(z):
    return (T_W - T_0) * np.power(z, P) + T_0


def u_p(z):
    return 3.084 * 1e-4 / (np.exp(4.799 * 1e+4 / T(z)) - 1)


def k1(T):
    xi = np.log(T_TABLE)
    eta = np.log(K_TABLE_1)
    return np.exp(np.interp(np.log(T), xi, eta))


def k2(T):
    xi = np.log(T_TABLE)
    eta = np.log(K_TABLE_2)
    return np.exp(np.interp(np.log(T), xi, eta))


def runge_kutta(cur_h, z_n, y_n, u_n, f, phi):
    k1 = f(z_n, y_n, u_n)
    q1 = phi(z_n, y_n, u_n)
    k2 = f(z_n + cur_h / 2, y_n + cur_h * k1 / 2, u_n + cur_h * q1 / 2)
    q2 = phi(z_n + cur_h / 2, y_n + cur_h * k1 / 2, u_n + cur_h * q1 / 2)
    k3 = f(z_n + cur_h / 2, y_n + cur_h * k2 / 2, u_n + cur_h * q2 / 2)
    q3 = phi(z_n + cur_h / 2, y_n + cur_h * k2 / 2, u_n + cur_h * q2 / 2)
    k4 = f(z_n + cur_h, y_n + cur_h * k3, u_n + cur_h * q3)
    q4 = phi(z_n + cur_h, y_n + cur_h * k3, u_n + cur_h * q3)

    y = y_n + cur_h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    u = u_n + cur_h / 6 * (q1 + 2 * q2 + 2 * q3 + q4)

    return y, u


def get_solution(k, t0, f0, up0):
    """z in [0; 1]"""

    # U'
    def phi(z, y, t):
        return -3.0 * R * k(T(z)) * y / C / up0

    # F'
    def f(z, y, t):
        if abs(z) < 1e-4:
            return R * C * k(T(z)) * (u_p(z) - t * up0) / 2

        return R * C * k(T(z)) * (u_p(z) - t * up0) - y / z

    y = np.zeros(len(z))  # -> f
    t = np.zeros(len(z))  # -> t = u / u_p(0)
    y[0] = f0
    t[0] = t0
    cur_h = h

    for i in range(1, len(z)):
        z_n = z[i - 1]
        y_n = y[i - 1]
        t_n = t[i - 1]
        k1 = f(z_n, y_n, t_n)
        q1 = phi(z_n, y_n, t_n)
        k2 = f(z_n + h / 2, y_n + h * k1 / 2, t_n + h * q1 / 2)
        q2 = phi(z_n + h / 2, y_n + h * k1 / 2, t_n + h * q1 / 2)
        k3 = f(z_n + h / 2, y_n + h * k2 / 2, t_n + h * q2 / 2)
        q3 = phi(z_n + h / 2, y_n + h * k2 / 2, t_n + h * q2 / 2)
        k4 = f(z_n + h, y_n + h * k3, t_n + h * q3)
        q4 = phi(z_n + h, y_n + h * k3, t_n + h * q3)

        y[i] = y_n + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t[i] = t_n + h / 6 * (q1 + 2 * q2 + 2 * q3 + q4)


    return z, y, t


def solve(k, u0, f0):
    up0 = u_p(0)
    t0 = u0 / up0
    z, f, t = get_solution(k, t0, f0, up0)
    u = t * up0
    # z, f, u = get_solution(k, u0, f0, up0)
    return z, f, u


def psi(f, u):
    return f - 0.393 * C * u


def psi_end(u0, k):
    _, f, u = solve(k, u0, 0)
    return psi(f[-1], u[-1])


def get_u0(k):
    start = 0.1 * u_p(0)
    stop = 1 * u_p(0)
    cur = 0

    while abs((stop - start) / stop) > DI_EPS:
        cur = (start + stop) / 2

        if psi_end(start, k) * psi_end(cur, k) < 0:
            stop = cur
        else:
            start = cur

    return cur, start, stop


def get_solution_2(k, u0, f0, up0):
    """z in [0; 1]"""
    # U'
    def phi(z, y, t):
        return -3.0 * R * k(T(z)) * y / C

    # F'
    def f(z, y, t):
        if abs(z) < 1e-4:
            return R * C * k(T(z)) * (u_p(z) - t) / 2

        return R * C * k(T(z)) * (u_p(z) - t) - y / z

    cur_h = h
    cur_z = min_z
    z_n = max_z
    cur_y = f0
    cur_u = u0
    z = []
    y = []
    u = []
    z.append(cur_z)
    y.append(cur_y)
    u.append(cur_u)

    while abs(cur_z - z_n) > 1e-6:
        y_h, u_h = runge_kutta(cur_h, cur_z, cur_y, cur_u, f, phi)

        new_h = cur_h / 2
        y_2h_1, u_2h_1 = runge_kutta(new_h, cur_z, cur_y, cur_u, f, phi)
        y_2h, u_2h = runge_kutta(new_h, cur_z + new_h, y_2h_1, u_2h_1, f, phi)

        eps_y = (y_2h - y_h) / y_2h
        eps_u = (u_2h - u_h) / u_2h

        new_y = y_h
        new_u = u_h
        print(f"EPS_Y: {eps_y:.5e} EPS_U: {eps_u:.5e}")

        while abs(eps_y) > 1e-4 or abs(eps_u) > 1e-4:
            new_h = cur_h / 2
            tmp_z = cur_z + new_h
            tmp_z_end = cur_z + cur_h
            tmp_y_res = y_h
            tmp_u_res = u_h
            print(f"NEW_H: {new_h:.3e}")

            while abs(tmp_z - (new_h + tmp_z_end)) > 1e-6:
                tmp_y_res, tmp_u_res = runge_kutta(new_h, tmp_z, tmp_y_res, tmp_u_res, f, phi)
                tmp_z += new_h

            eps_y = (tmp_y_res - y_h) / y_2h
            eps_u = (tmp_u_res - u_h) / u_2h
            print(f"EPS_Y: {eps_y:.5e} EPS_U: {eps_u:.5e}")
            # print(eps_y, "---", eps_u)

            cur_h = new_h
            new_y = tmp_y_res
            new_u = tmp_u_res

        cur_y = new_y
        cur_u = new_u
        # print(f"CUR_Y: {cur_y:.5e} CUR_U: {cur_u:.5e}")
        y.append(new_y)
        u.append(new_u)
        cur_z += cur_h
        z.append(cur_z)

    # print("CUR_H: ", cur_h)

    return z, y, u


def solve_2(k, u0, f0):
    up0 = u_p(0)
    z, f, u = get_solution_2(k, u0, f0, up0)
    return z, f, u


def run(n_k):
    def graph(cur_z):
        T_z = np.array(list(map(T, z)))
        k = np.array(list(map(k_t, T_z)))

        fig, axs = plt.subplots(2, 3)
        for i in range(2):
            for j in range(3):
                axs[i][j].grid()

        axs[0][0].plot(cur_z, f, color="red")  # , label="F(z)")
        axs[0][0].set_title("F(z)")
        axs[0][1].plot(cur_z, u, color="red")  # , label="U(z)")
        axs[0][1].set_title("U(z)")
        axs[0][2].plot(cur_z, u_p(z), color="red")  # , label="U_p(z)")
        axs[0][2].set_title("U_p(z)")
        axs[1][0].plot(z, T(z), color="red")  # , label="T(z)")
        axs[1][0].set_title("T(z)")
        axs[1][1].plot(z, k, color="red")  # , label="k(z)")
        axs[1][1].set_title("k1(z)")

        # _p = np.vectorize(psi)(f, u)
        # axs[1][2].plot(cur_z, _p, color="red")  # , label="k(z)")
        # axs[1][2].set_title(f"psi(z), u0 = {u0:.3e}")

        axs[1][2].plot(cur_z, u, color="red")  # , label="U(z)")
        axs[1][2].plot(cur_z, u_p(z), color="blue")  # , label="U_p(z)")
        axs[1][2].set_title("U(z), U_p(z)")

        fig.set_size_inches(12, 7)
        fig.tight_layout()
        plt.show()

    k_t = k1 if n_k == 1 else k2
    # solve_f = solve if n_k == 1 else solve_2

    if n_k == 1:
        u0, u0_min, u0_max = get_u0(k_t)
        print(f"u0 interval: [{u0_min:.4e}; {u0_max:.4e}]\nu0 = {u0:.4e}")
    else:
        u0 = 3.9207e-07

    all_z, f, u = solve(k_t, u0, 0)
    graph(all_z)

    # all_z, f, u = solve_2(k_t, u0, 0)
    # graph(all_z)




if __name__ == '__main__':
    n_k = 1
    run(n_k)
