import numpy as np
import matplotlib.pyplot as plt

T_TABLE = np.array([2000, 3000, 4000,
                    5000, 6000, 7000,
                    8000, 9000, 10000])

K_TABLE_1 = np.array([8.200E-03, 2.768E-02, 6.560E-02, 1.281E-01,
                    2.214E-01, 3.516E-01, 5.248E-01, 7.472E-01,
                    1.025E+00]) * 1000000

K_TABLE_2 = np.array([1.600e+00, 5.400e+00, 1.280e+01,
                      2.500e+01, 4.320e+01, 6.860e+01,
                      1.024e+02, 1.458e+02, 2.000e+02])

DI_EPS = 1e-4
R = 0.35
T_W = 2000  # 10000
T_0 = 10000
P = 4
c = 3e+10
min_z = 0
max_z = 1


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


def sweep(z_start, z_end, N, k):
    z, h = np.linspace(z_start, z_end, N, retstep=True)

    def lmbd(z):
        return c / (3 * k(T(z)))

    def p(z):
        return c * k(T(z))

    def f(z):
        return c * k(T(z)) * u_p(z)

    def V(z):
        return ((z + h / 2) ** 2 - (z - h / 2) ** 2) / 2

    def kappa(z):
        l1 = lmbd(z - h / 2)
        l2 = lmbd(z + h / 2)
        # return 2 * l1 * l2 / (l1 + l2)
        return (l1 + l2) / 2

    def start_MKP():
        #  F0 = 0

        z_1_2 = z_start + h / 2
        alpha = z_1_2 * kappa(z_1_2) / (R ** 2 * h)
        beta = h * p(z_1_2) * z_1_2 / 8

        M0 = -alpha - beta
        K0 = alpha - beta
        P0 = -f(z_1_2) * z_1_2 * h / 4

        return M0, K0, P0

    def end_MKP():
        #  FN = 0.393 * c * zN

        z_N_1_2 = z_end - h / 2
        alpha = z_N_1_2 * kappa(z_N_1_2) / (R ** 2 * h)
        beta = h * p(z_N_1_2) * z_N_1_2 / 8

        MN = alpha - beta
        KN = -alpha - beta - p(z_end) * z_end * h / 4 - 0.393 * z_end * c / R
        PN = -(f(z_N_1_2) * z_N_1_2 + f(z_end) * z_end) * h / 4

        return MN, KN, PN

    def A(z):
        return (z - h / 2) * kappa(z - h / 2) / (R ** 2 * h)

    def C(z):
        return (z + h / 2) * kappa(z + h / 2) / (R ** 2 * h)

    def B(z):
        return A(z) + C(z) + p(z) * V(z)

    def D(z):
        return f(z) * V(z)

    xi = np.zeros(N)
    eta = np.zeros(N)

    M0, K0, P0 = start_MKP()

    xi[0] = -K0 / M0
    eta[0] = P0 / M0

    for i in range(1, N):
        z_n = z[i]

        xi_n_1 = xi[i - 1]
        eta_n_1 = eta[i - 1]

        A_n = A(z_n)
        C_n = C(z_n)
        B_n = B(z_n)
        D_n = D(z_n)

        xi[i] = C_n / (B_n - A_n * xi_n_1)
        eta[i] = (A_n * eta_n_1 + D_n) / (B_n - A_n * xi_n_1)

    MN, KN, PN = end_MKP()
    y = np.zeros(N)
    y[N - 1] = (PN - MN * eta[N - 1]) / (MN * xi[N - 1] + KN)

    for i in range(N - 2, -1, -1):
        y[i] = xi[i + 1] * y[i + 1] + eta[i + 1]

    return z, y, h


def get_u_derivative(u, h):
    n = len(u)
    u_der = np.zeros(n)
    u_der[0] = (-3 * u[0] + 4 * u[1] - u[2]) / (2 * h)

    for i in range(1, n - 1):
        u_der[i] = (u[i + 1] - u[i - 1]) / (2 * h)

    u_der[n - 1] = (u[n - 3] - 4 * u[n - 2] + 3 * u[n - 1]) / (2 * h)

    return u_der


def get_f_from_u(u, z, h, k):
    n = len(u)
    f = np.zeros(n)
    der_u = get_u_derivative(u, h)

    for i in range(1, n):
        f[i] = -c * der_u[i] / (3 * R * k(T(z[i])))

    return f


def div_f(k, z, u):
    return c * k(T(z)) * (u_p(z) - u)


def run(k_n):
    n = 2000

    k = k1 if k_n == 1 else k2

    z, u, h = sweep(min_z, max_z, n, k)
    der_u = get_u_derivative(u, h)
    f = get_f_from_u(u, z, h, k)

    T_z = np.array(list(map(T, z)))
    _k = np.array(list(map(k, T_z)))

    fig, axs = plt.subplots(2, 3)
    for i in range(2):
        for j in range(3):
            axs[i][j].grid()

    axs[0][0].plot(z, f, color="red")  # , label="F(z)")
    axs[0][0].set_title("F(z)")
    axs[0][1].plot(z, u, color="red")  # , label="U(z)")
    axs[0][1].set_title("U(z)")
    axs[0][2].plot(z, u, color="red", label="U")  # , label="U(z)")
    axs[0][2].plot(z, u_p(z), color="blue", label="U_p")  # , label="U_p(z)")
    axs[0][2].set_title("U(z), U_p(z)")
    axs[0][2].legend(loc="best")
    axs[1][0].plot(z, div_f(k, z, u), color="red")  # , label="T(z)")
    axs[1][0].set_title("divF(z)")
    axs[1][1].plot(z, _k, color="red")  # , label="k(T(z))")
    axs[1][1].set_title(f"k{k_n}(z)")
    axs[1][2].plot(z, der_u, color="red")  # , label="k(T(z))")
    axs[1][2].set_title(f"u'(z)")

    fig.set_size_inches(12, 7)
    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    run(1)
