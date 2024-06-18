import numpy as np
import matplotlib.pyplot as plt

a1 = 0.0134
b1 = 1
c1 = 4.35e-4
m1 = 1
u0 = 300
a = 10
b = 10
x0 = 5
z0 = 5
alpha_1 = 0.05
alpha_2 = 0.1
alpha_3 = 0.5
alpha_4 = 1.0
f0 = 1
beta = 1
F0 = 30

def der_lambda(y):
    return a1 * c1 * m1 * (y ** (m1 - 1))

def _lambda(y):
    return a1 * (b1 + c1 * (y ** m1))

def _f(x, z):
    return f0 * np.exp(-beta * ((x - x0) ** 2 + (z - z0) ** 2))

def kappa(y1, y2):
    return (_lambda(y1) + _lambda(y2)) / 2


def x_left_bound_cond():
    xk0 = 1
    xm0 = 0
    xp0 = u0
    return xk0, xm0, xp0


def x_right_bound_cond():
    xkN = 0
    xmN = 1
    xpN = u0
    return xkN, xmN, xpN


def x_right_sweep(curr_z, y_prev, y_curr, x_0, xN, h_x, tau):
    eps = 1e-2
    flag = True
    approx_y_curr = y_curr
    delta_y = []

    while flag:
        k0, m0, p0 = x_left_bound_cond()
        kN, mN, pN = x_right_bound_cond()

        ksi = [0, -m0 / k0]
        eta = [0, p0 / k0]
        x = h_x
        n = 1

        for i in range(1, len(approx_y_curr) - 1):
            aN = tau / (h_x ** 2) * kappa(approx_y_curr[i - 1], approx_y_curr[i])
            dN = tau / (h_x ** 2) * kappa(approx_y_curr[i], approx_y_curr[i + 1])
            bN = tau / (h_x ** 2) + aN + dN
            fN = _f(x, curr_z) * tau / 2 + y_prev[i]

            aN_1 = tau / (2 * h_x ** 2) * der_lambda(approx_y_curr[i - 1]) * (approx_y_curr[i - 1] - approx_y_curr[i]) + aN
            bN_1 = tau / (2 * h_x ** 2) * der_lambda(approx_y_curr[i]) * (-approx_y_curr[i - 1] + 2 * approx_y_curr[i] - approx_y_curr[i + 1]) + bN
            dN_1 = tau / (2 * h_x ** 2) * der_lambda(approx_y_curr[i + 1]) * (-approx_y_curr[i] + approx_y_curr[i + 1]) + dN
            fN_1 = aN * approx_y_curr[i - 1] - bN * approx_y_curr[i] + dN * approx_y_curr[i + 1] + fN

            ksi.append(dN_1 / (bN_1 - aN_1 * ksi[n]))
            eta.append((aN_1 * eta[n] + fN_1) / (bN_1 - aN_1 * ksi[n]))

            n += 1
            x += h_x

        delta_y = [0] * (n + 1)
        delta_y[n] = (pN - kN * eta[n]) / (kN * ksi[n] + mN)

        for i in range(n - 1, -1, -1):
            delta_y[i] = ksi[i + 1] * delta_y[i + 1] + eta[i + 1]

        stop_approx_flag = True
        for i in range(len(delta_y)):
            # print("x sweep: ", abs(delta_y[i] / approx_y_curr[i]))
            if (abs(delta_y[i] / approx_y_curr[i]) >= eps):
                print(f"[i={i}] x sweep: ", abs(delta_y[i] / approx_y_curr[i]))
                stop_approx_flag = False
                break

        if stop_approx_flag:
            flag = False

        # !!! см. уравнение которое получено после линеаризации
        approx_y_curr = [approx_y_curr[i] + delta_y[i] for i in range(len(delta_y))]
        print('..........')

    print(approx_y_curr)
    exit(0)
    return approx_y_curr


def z_left_bound_cond():
    zk0 = 1
    zm0 = 0
    zp0 = u0
    return zk0, zm0, zp0


def z_right_bound_cond():
    zkN = 0
    zmN = 1
    zpN = u0
    return zkN, zmN, zpN


def z_right_sweep(curr_x, y_prev, y_curr, z_0, zM, h_z, tau):
    eps = 1e-4
    flag = True
    prev_res = y_curr
    y = []

    while flag:
        k0, m0, p0 = z_left_bound_cond()
        kN, mN, pN = z_right_bound_cond()

        ksi = [0, -m0 / k0]
        eta = [0, p0 / k0]
        z = h_z
        n = 1

        for i in range(1, len(prev_res) - 1):
            aN = tau / (h_z ** 2) * kappa(prev_res[i - 1], prev_res[i])
            dN = tau / (h_z ** 2) * kappa(prev_res[i], prev_res[i + 1])
            bN = tau / (h_z ** 2) + aN + dN
            fN = _f(z, curr_x) * tau / 2 + y_prev[i]

            aN_1 = tau / (2 * h_z ** 2) * der_lambda(prev_res[i - 1]) * (prev_res[i - 1] - prev_res[i]) + aN
            bN_1 = tau / (2 * h_z ** 2) * der_lambda(prev_res[i]) * (-prev_res[i - 1] + 2 * prev_res[i] - prev_res[i + 1]) + bN
            dN_1 = tau / (2 * h_z ** 2) * der_lambda(prev_res[i + 1]) * (-prev_res[i] + prev_res[i + 1]) + dN
            fN_1 = aN * prev_res[i - 1] - bN * prev_res[i] + dN * prev_res[i + 1] + fN

            ksi.append(dN_1 / (bN_1 - aN_1 * ksi[n]))
            eta.append((aN_1 * eta[n] + fN_1) / (bN_1 - aN_1 * ksi[n]))

            n += 1
            z += h_z

        y = [0] * (n + 1)
        y[n] = (pN - kN * eta[n]) / (kN * ksi[n] + mN)

        for i in range(n - 1, -1, -1):
            y[i] = ksi[i + 1] * y[i + 1] + eta[i + 1]

        count = 0
        for i in range(len(y)):
            if abs((y[i] - prev_res[i]) / y[i]) < eps:
                count += 1

        if count == len(y):
            flag = False

        prev_res = y

    return y


def is_fault_zero(y_delta, cur_y, eps):
    count = 0
    for i in range(len(y_delta)):
        if abs(y_delta[i] / cur_y[i]) < eps:
            count += 1

    if count == len(y_delta):
        return True

    return False


def get_approx_1(approx_0, x_0, xN, h_x, z_0, zM, h_z, tau):
    eps = 1e-4
    curr_z = z_0
    approx_1 = []
    i = 0

    while curr_z < zM - h_z:
        curr_approx_1 = []
        curr_approx_0 = []
        for approx_0_z in approx_0:
            curr_approx_1.append(approx_0_z[i])
            curr_approx_0.append(approx_0_z[i])

        fault_flag = False
        next_approx_1 = []

        while not fault_flag:
            curr_approx_delta_1 = x_right_sweep(curr_z, curr_approx_0, curr_approx_1, x_0, xN, h_x, tau)
            next_approx_1 = curr_approx_1 + curr_approx_delta_1
            fault_flag = is_fault_zero(curr_approx_delta_1, next_approx_1, eps)
            curr_approx_1 = next_approx_1

        curr_z += h_z
        i += 1
        approx_1.append(next_approx_1)

    return approx_1


def get_approx_2(approx_1, x_0, xN, h_x, z_0, zM, h_z, tau):
    eps = 1e-4
    curr_x = x_0
    approx_2 = []
    i = 0

    while curr_x < xN - h_x:
        curr_approx_2 = approx_1[i]
        curr_approx_1 = approx_1[i]

        fault_flag = False
        next_approx_2 = []

        while not fault_flag:
            curr_approx_delta_2 = z_right_sweep(curr_x, curr_approx_1, curr_approx_2, z_0, zM, h_z, tau)
            next_approx_2 = curr_approx_2 + curr_approx_delta_2
            fault_flag = is_fault_zero(curr_approx_delta_2, next_approx_2, eps)
            curr_approx_2 = next_approx_2

        curr_x += h_x
        i += 1
        approx_2.append(next_approx_2)

    return approx_2


def get_layer_solution(approx_0, x_0, xN, h_x, z_0, zM, h_z, tau):
    approx_1 = get_approx_1(approx_0, x_0, xN, h_x, z_0, zM, h_z, tau)
    approx_2 = get_approx_2(approx_1, x_0, xN, h_x, z_0, zM, h_z, tau)
    return approx_2


def is_left_side_zero(y_t_m, y_t_m_1, eps):
    count = 0
    for i in range(len(y_t_m)):
        if abs((y_t_m[i] - y_t_m_1[i]) / y_t_m_1[i]) < eps:
            count += 1

    if count == len(y_t_m):
        return True

    return False


def get_solution(x_0, xN, h_x, z_0, zM, h_z, t_0, tau):
    eps = 1e-4

    approx_0 = []
    tmp_x = x_0
    while tmp_x < xN - h_x:
        newapprox_0 = []
        tmp_z = z_0
        while tmp_z < zM - h_z:
            newapprox_0.append(u0)
            tmp_z += h_z
        approx_0.append(newapprox_0)
        tmp_x += h_x

    prev_res = approx_0.copy()
    res = []
    times = []
    cur_t = t_0
    left_side_zero = False

    while not left_side_zero:
        res = get_layer_solution(prev_res, x_0, xN, h_x, z_0, zM, h_z, tau)
        times.append(cur_t)

        count = 0
        for i in range(len(res)):
            if is_left_side_zero(prev_res[i], res[i], eps):
                count += 1

        if count == len(res):
            left_side_zero = True

        prev_res = res
        cur_t += tau

    zs = np.arange(z_0, zM + h_z, h_z)
    xs = np.arange(x_0, xN + h_x, h_x)

    return np.array(xs), np.array(zs), np.array(times), np.array(res)

x_0 = 0
z_0 = 0
xN = a
zM = b
t_0 = 0
N = 100
M = 100
T = 100

# hx = (xN - x_0) / N
# hz = (zM - z_0) / M
hx = 0.1
hz = 0.1
tau = 1

xs, zs, times, res = get_solution(x_0, xN, hx, z_0, zM, hz, t_0, tau)