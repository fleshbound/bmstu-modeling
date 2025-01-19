import matplotlib.pyplot as plt


def draw_poisson(x, y_function, y_density, name):
    fig, axs = plt.subplots(nrows=2, figsize=(5, 10))

    fig.suptitle(name)
    axs[0].plot(x, y_function, 'v', color='red', ms=4)
    axs[1].vlines(x, 0, y_density, color='red', lw=4, alpha=0.5)

    axs[0].set_xlabel('x')
    axs[0].set_ylabel('F(x)')
    axs[1].set_xlabel('x')
    axs[1].set_ylabel('f(x)')

    axs[0].grid(True)
    axs[1].grid(True)
    plt.show()


def draw_universal(x, y_function, y_density, name):
    fig, axs = plt.subplots(nrows=2, figsize=(5, 10))

    fig.suptitle(name)
    axs[0].plot(x, y_function, color='red')
    axs[1].plot(x, y_density, color='red')

    axs[0].set_xlabel('x')
    axs[0].set_ylabel('F(x)')
    axs[1].set_xlabel('x')
    axs[1].set_ylabel('f(x)')

    axs[0].grid(True)
    axs[1].grid(True)
    plt.show()
