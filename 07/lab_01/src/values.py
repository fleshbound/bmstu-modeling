from scipy.stats import poisson

def ud_function(a, b, x):
    return (x - a) / (b - a) if a <= x < b else 0 if x < a else 1


def ud_density(a, b, x):
    return 1 / (b - a) if a <= x <= b else 0


def poisson_function(x, mu):
    poisson.cdf(x, mu)
    return poisson.cdf(x, mu)


def poisson_density(x, mu):
    return poisson.pmf(x, mu)
