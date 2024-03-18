import numpy as np
import math
import matplotlib.pyplot as plt

T_TABLE = np.array([2000, 3000, 4000,
                    5000, 6000, 7000,
                    8000, 9000, 10000])

K_TABLE_1 = np.array([8.200E-03, 2.768E-02, 6.560E-02, 1.281E-01,
                    2.214E-01, 3.516E-01, 5.248E-01, 7.472E-01,
                    1.025E+00]) * 28

K_TABLE_2 = np.array([1.600e+00, 5.400e+00, 1.280e+01,
                      2.500e+01, 4.320e+01, 6.860e+01,
                      1.024e+02, 1.458e+02, 2.000e+02])

DI_EPS = 1e-4
R = 0.35
T_W = 2000
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


def get_t_solution(k, t0, f0, up0):
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

    return y, t


def solve(k, u0, f0):
    up0 = u_p(0)
    # up0 = 1
    t0 = u0 / up0
    F, t = get_t_solution(k, t0, f0, up0)
    u = t * up0
    return F, u


def psi_end(u0, k):
    f, u = solve(k, u0, 0)
    return f[-1] - 0.393 * C * u[-1]


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


def run(n_k):
    k_t = k1 if n_k == 1 else k2

    u0, u0_min, u0_max = get_u0(k_t)
    print(f"u0 interval: [{u0_min:.4e}; {u0_max:.4e}]\nu0 = {u0:.4e}")
    
    f, u = solve(k_t, u0, 0)
    T_z = np.array(list(map(T, z)))
    k = np.array(list(map(k_t, T_z)))

    fig, axs = plt.subplots(2, 3)
    for i in range(2):
        for j in range(3):
            if i != 1 or j != 2:
                axs[i][j].grid()

    axs[0][0].plot(z, f, color="red")  # , label="F(z)")
    axs[0][0].set_title("F(z)")
    axs[0][1].plot(z, u, color="red")  # , label="U(z)")
    axs[0][1].set_title("U(z)")
    axs[0][2].plot(z, u_p(z), color="red")  # , label="U_p(z)")
    axs[0][2].set_title("U_p(z)")
    axs[1][0].plot(z, T(z), color="red")  # , label="T(z)")
    axs[1][0].set_title("T(z)")
    axs[1][1].plot(z, k, color="red")  # , label="k(z)")
    axs[1][1].set_title("k1(z)")

    fig.set_size_inches(12, 7)
    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    n_k = 1
    run(n_k)
