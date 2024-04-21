import numpy as np
import matplotlib.pyplot as plt

T_TABLE = np.array([2000, 3000, 4000,
                    5000, 6000, 7000,
                    8000, 9000, 10000])

K_TABLE_1 = np.array([8.200E-03, 2.768E-02, 6.560E-02, 1.281E-01,
                    2.214E-01, 3.516E-01, 5.248E-01, 7.472E-01,
                    1.025E+00]) * 2

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


k = k1


def sweep(z_start, z_end, N, F_0, F_N):
    z, h = np.linspace(z_start, z_end, N, retstep=True)

    def lmbd(z):
        return - C / (3 * k(z))

    def p(z):
        return - C * k(z)

    def f(z):
        return - C * k(z) * u_p(z)

    def kappa(z):
        """n - index in z"""
        a = lmbd(z - h / 2)
        b = lmbd(z - h / 2)
        return (a + b) / 2
        # return 2 * a * b / (a + b)


    def start_MNP():
        z_1_2 = z_start + h / 2
        alpha = z_1_2 * kappa(z_1_2) / (R ** 2  * h)
        beta = h * p(z_1_2) * z_1_2 / 8

        M0 = -alpha - beta - p(z_start) * z_start * h / 4
        K0 = alpha - beta
        # P0 =


    # xi = np.zeros(N + 1)
    # eta = np.zeros(N + 1)
    #
    # kappa_1_2 = kappa(z, 0)
    # z_1_2 = z[0] + h / 2
    # z_0 = z[0]
    # p_1_2 = p(z_1_2)
    # f_1_2 = f(z_1_2)
    # f_0 = f(z_0)
    # p_0 = p(z_0)
    # K_0 = z_1_2 * kappa_1_2 / (R * R * h) + h * p_1_2 * z_1_2 / 8
    # M_0 = - z_1_2 * kappa_1_2 / (R * R * h) - h * p_1_2 * z_1_2 / 8 + p_0 * z_0 * h / 4
    # P_0 = - F_0 * z_0 / R - (f_1_2 * z_1_2 + f_0 * z_0) * h / 4
    #
    # xi[1] = - K_0 / M_0
    # eta[1] = P_0
    #
    # for n in range(2, N + 1):
    #     z_n_p1_2 = z[n] + h / 2
    #     z_n_m1_2 = z[n] - h / 2
    #     V_n = (z_n_p1_2 ** 2 - z_n_m1_2 ** 2) / 2
    #     z_n = z[n]
    #     p_n = p(z_n)
    #     f_n = f(z_n)
    #     kappa_n_p1_2 = kappa_plus(z, n)
    #     kappa_n_m1_2 = kappa_minus(z, n)
    #
    #     C_n = z_n_p1_2 * kappa_n_p1_2 / (R * R * h)
    #     A_n = z_n_m1_2 * kappa_n_m1_2 / (R * R * h)
    #     B_n = A_n + C_n - p_n * V_n
    #     D_n = f_n * V_n
    #
    #     xi[n] = C_n / (B_n - A_n * xi[n - 1])
    #     eta[n] = (A_n * eta[n - 1] + D_n) / (B_n - A_n * xi[n - 1])
    #
    # kappa_N_1_2 = kappa_minus(z, N - 1)
    # z_N_1_2 = z[N - 1] - h / 2
    # z_N = z[N - 1]
    # p_N_1_2 = p(z_N_1_2)
    # f_N_1_2 = f(z_N_1_2)
    # f_N = f(z_N)
    # p_N = p(z_N)
    # K_N = - z_N_1_2 * kappa_N_1_2 / (R * R * h) + h * p_N_1_2 * z_N_1_2 / 8 - p_N * z_N * h / 4
    # M_N = z_N_1_2 * kappa_N_1_2 / (R * R * h) + h * p_N_1_2 * z_N_1_2 / 8
    # P_N = F_N * z_N / R - (f_N_1_2 * z_N_1_2 + f_N * z_N) * h / 4
    #
    # y = np.zeros(N)
    # y[N - 1] = (P_N - M_N * eta[N - 1]) / (M_N * xi[N - 1] + K_N)
    #
    # for n in range(N - 2, 0, -1):
    #     y[n] = xi[n + 1] * y[n + 1] + eta[n + 1]

    return z, y


def run(n_k):
    def graph(cur_z):
        T_z = np.array(list(map(T, z)))
        k = np.array(list(map(k_t, T_z)))

        fig, axs = plt.subplots(2, 3)
        for i in range(2):
            for j in range(3):
                axs[i][j].grid()

        # axs[0][0].plot(cur_z, f, color="red")  # , label="F(z)")
        # axs[0][0].set_title("F(z)")
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

    all_z, u = sweep(z, h, 0, 0)
    graph(all_z)

    # all_z, f, u = solve_2(k_t, u0, 0)
    # graph(all_z)

