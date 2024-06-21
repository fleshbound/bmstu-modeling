import numpy as np
import matplotlib.pyplot as plt

a1 = 0.0134
b1 = 1
c1 = 4.35e-4
m1 = 1
u0 = 300
a = 10
b = 10
x0 = 2
z0 = 2
f0 = 1000
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


def get_delta_y_by_x(curr_z, y_prev, y_curr, x_0, xN, h_x, tau):
    """Right sweep at curr_z (by x)"""
    # print("...call get_delta_y_by_x")

    k0, m0, p0 = x_left_bound_cond()
    kN, mN, pN = x_right_bound_cond()

    ksi = [0, -m0 / k0]
    eta = [0, p0 / k0]
    x = h_x
    n = 1

    for i in range(1, len(y_curr) - 1):
        aN = tau / (h_x ** 2) * kappa(y_curr[i - 1], y_curr[i])
        dN = tau / (h_x ** 2) * kappa(y_curr[i], y_curr[i + 1])
        bN = tau / (h_x ** 2) + aN + dN
        fN = _f(x, curr_z) / 2 + y_prev[i]

        aN_1 = tau / (2 * h_x ** 2) * der_lambda(y_curr[i - 1]) * (y_curr[i - 1] - y_curr[i]) + aN
        bN_1 = tau / (2 * h_x ** 2) * der_lambda(y_curr[i]) * (-y_curr[i - 1] + 2 * y_curr[i] - y_curr[i + 1]) + bN
        dN_1 = tau / (2 * h_x ** 2) * der_lambda(y_curr[i + 1]) * (-y_curr[i] + y_curr[i + 1]) + dN
        fN_1 = aN * y_curr[i - 1] - bN * y_curr[i] + dN * y_curr[i + 1] + fN

        ksi.append(dN_1 / (bN_1 - aN_1 * ksi[n]))
        eta.append((aN_1 * eta[n] + fN_1) / (bN_1 - aN_1 * ksi[n]))

        n += 1
        x += h_x

    delta_y = [0] * (n + 1)
    delta_y[n] = (pN - kN * eta[n]) / (kN * ksi[n] + mN)

    for i in range(n - 1, -1, -1):
        delta_y[i] = ksi[i + 1] * delta_y[i + 1] + eta[i + 1]

    # print("get_delta_y_by_x complete")

    return delta_y


def get_approx_y_by_x_fixed_z(curr_z, y_prev, start_approx_y_curr, x_0, xN, h_x, tau):
    """Simple iterations at curr_z (by x)"""
    # print(f"...call get_approx_y_by_x_fixed_z (curr_z = {curr_z})")

    eps = 1e-1
    run = True
    approx_y_curr = start_approx_y_curr
    iter_num = 0

    while run:
        approx_delta_y = get_delta_y_by_x(curr_z, y_prev, approx_y_curr, x_0, xN, h_x, tau)

        stop = True

        for i in range(len(approx_delta_y)):
            if abs(approx_delta_y[i] / approx_y_curr[i]) >= eps:
                stop = False
                break

        iter_num += 1
        run = not stop
        approx_y_curr = [approx_y_curr[i] + approx_delta_y[i] for i in range(len(approx_delta_y))]

    # print(f"get_approx_y_by_x_fixed_z complete ({iter_num} iterations)")

    return approx_y_curr


def get_approx_y_by_x(approx_0, x_0, xN, h_x, z_0, zM, h_z, tau):
    """Simple iterations for every z by x; \n
    returns matrix zs by xs"""
    print(f"...call get_approx_y_by_x")

    approx_1 = []
    curr_z = z_0
    i = 0

    while curr_z < zM - h_z:
        # массив массивов z-ов для каждого x
        curr_approx_0 = [approx_0[idx][i] for idx in range(len(approx_0))]
        curr_approx_1 = [approx_0[idx][i] for idx in range(len(approx_0))]

        next_approx_1 = get_approx_y_by_x_fixed_z(curr_z, curr_approx_0, curr_approx_1, x_0, xN, h_x, tau)
        approx_1.append(next_approx_1)
        curr_z += h_z
        i += 1

    print(f"get_approx_y_by_x complete")

    return np.array(approx_1)


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


def get_delta_y_by_z(curr_x, y_prev, y_curr, z_0, zM, h_z, tau):
    """Right sweep at curr_x (by z)"""
    # print("...call get_delta_y_by_z")

    k0, m0, p0 = z_left_bound_cond()
    kN, mN, pN = z_right_bound_cond()

    ksi = [0, -m0 / k0]
    eta = [0, p0 / k0]
    z = h_z
    n = 1

    for i in range(1, len(y_curr) - 1):
        aN = tau / (h_z ** 2) * kappa(y_curr[i - 1], y_curr[i])
        dN = tau / (h_z ** 2) * kappa(y_curr[i], y_curr[i + 1])
        bN = tau / (h_z ** 2) + aN + dN
        fN = _f(z, curr_x) / 2 + y_prev[i]

        aN_1 = tau / (2 * h_z ** 2) * der_lambda(y_curr[i - 1]) * (y_curr[i - 1] - y_curr[i]) + aN
        bN_1 = tau / (2 * h_z ** 2) * der_lambda(y_curr[i]) * (-y_curr[i - 1] + 2 * y_curr[i] - y_curr[i + 1]) + bN
        dN_1 = tau / (2 * h_z ** 2) * der_lambda(y_curr[i + 1]) * (-y_curr[i] + y_curr[i + 1]) + dN
        fN_1 = aN * y_curr[i - 1] - bN * y_curr[i] + dN * y_curr[i + 1] + fN

        ksi.append(dN_1 / (bN_1 - aN_1 * ksi[n]))
        eta.append((aN_1 * eta[n] + fN_1) / (bN_1 - aN_1 * ksi[n]))

        n += 1
        z += h_z

    delta_y = [0] * (n + 1)
    delta_y[n] = (pN - kN * eta[n]) / (kN * ksi[n] + mN)

    for i in range(n - 1, -1, -1):
        delta_y[i] = ksi[i + 1] * delta_y[i + 1] + eta[i + 1]

    # print("get_delta_y_by_z complete")

    return delta_y


def get_approx_y_by_z_fixed_x(curr_x, y_prev, start_approx_y_curr, z_0, zM, h_z, tau):
    """Simple iterations at curr_x (by z)"""
    # print(f"...call get_approx_y_by_z_fixed_x (curr_x = {curr_x})")

    eps = 1e-1
    run = True
    approx_y_curr = start_approx_y_curr
    iter_num = 0

    while run:
        approx_delta_y = get_delta_y_by_z(curr_x, y_prev, approx_y_curr, z_0, zM, h_z, tau)

        stop = True

        for i in range(len(approx_delta_y)):
            if abs(approx_delta_y[i] / approx_y_curr[i]) >= eps:
                stop = False
                break

        iter_num += 1
        run = not stop
        approx_y_curr = [approx_y_curr[i] + approx_delta_y[i] for i in range(len(approx_delta_y))]

    # print(f"get_approx_y_by_z_fixed_x complete ({iter_num} iterations)")

    return approx_y_curr


def get_approx_y_by_z(approx_1, x_0, xN, h_x, z_0, zM, h_z, tau):
    """Simple iterations for every x by z; \n
    returns matrix zs by xs"""
    print("...call get_approx_y_by_z")

    approx_2 = []
    curr_x = x_0
    i = 0

    while curr_x < xN - h_x:
        # approx_1 - массив массивов z-ов для каждого x
        curr_approx_1 = approx_1[i]
        curr_approx_2 = approx_1[i]

        next_approx_1 = get_approx_y_by_z_fixed_x(curr_x, curr_approx_1, curr_approx_2, z_0, zM, h_z, tau)
        approx_2.append(next_approx_1)
        curr_x += h_x
        i += 1

    print("get_approx_y_by_z complete")

    return np.array(approx_2)


def get_layer_solution(approx_0, x_0, xN, h_x, z_0, zM, h_z, tau):
    approx_1 = get_approx_y_by_x(approx_0, x_0, xN, h_x, z_0, zM, h_z, tau)
    approx_2 = get_approx_y_by_z(approx_1.T, x_0, xN, h_x, z_0, zM, h_z, tau)
    return approx_2


def is_left_side_zero(y_t_m, y_t_m_1):
    eps = 1e-1

    for i in range(len(y_t_m)):
        dt = np.mean(y_t_m - y_t_m_1) / tau
        mistake = abs(dt / tau)
        if mistake > eps:
            print(f"is_left_side_zero: !too big mistake: dt = {dt:.6f}, mistake = {mistake:.6f}")
            return False

    return True


def get_solution(x_0, xN, h_x, z_0, zM, h_z, t_0, tau):
    print("call get_solution")

    xs = np.arange(x_0, xN, h_x)
    zs = np.arange(z_0, zM, h_z)
    approx_0 = [[u0 for _ in zs] for _ in xs]
    prev_res = approx_0
    res = []
    times = []
    cur_t = t_0
    run = True
    iter_num = 0

    while run:
        res = get_layer_solution(prev_res, x_0, xN, h_x, z_0, zM, h_z, tau)
        times.append(cur_t)

        stop = True
        for i in range(len(res)):
            if not is_left_side_zero(prev_res[i], res[i]):
                stop = False
                break
        run = not stop

        prev_res = res.copy()
        cur_t += tau
        iter_num += 1

    print(f"get_solution ({iter_num} time iterations)")

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
tau = 0.1

xs1, zs1, times, res = get_solution(x_0, xN, hx, z_0, zM, hz, t_0, tau)
xs = np.linspace(x_0, xN, len(res))
zs = np.linspace(z_0, zM, len(res[0]))
X, Y = np.meshgrid(xs, zs)
farr = np.array([np.array(T_m) for T_m in res])

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.plot_surface(X, Y, farr, cmap="viridis")

ax.set_xlabel("x")
ax.set_ylabel("z")
ax.set_zlabel("U(x, z)")
plt.show()

plt.imshow(res, cmap='viridis')
plt.colorbar()
plt.show()