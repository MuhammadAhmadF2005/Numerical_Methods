import math
from math import cos, sin, tan, sinh, cosh, tanh, asin, acos, atan, asinh, acosh, atanh
from math import exp, log, log10, log2, e, pi


def differentiate(f, h=1e-7):
    """
    Approximates the derivative f'(x) of function f using
    the central difference formula:
        f'(x) ≈ (f(x + h) - f(x - h)) / (2 * h)
    """
    return lambda x: (f(x + h) - f(x - h)) / (2 * h)


derivative = differentiate  # Alias for convenience


def newton_raphson(f, *args, df=None, tol=1e-5, max_iter=100, **kwargs):
    """
    Implements the Newton-Raphson Method from Burden & Faires.
    If df (derivative) is not provided, it is automatically computed using
    numerical differentiation.
    Stops when |p - p0| < tol or f(p) == 0.
    """
    p0 = kwargs.get("p0", None)

    if len(args) == 1:
        if callable(args[0]):
            df = args[0]
        else:
            p0 = args[0]
    elif len(args) >= 2:
        if callable(args[0]):
            df, p0 = args[0], args[1]
        else:
            p0, df = args[0], args[1]

    if p0 is None:
        raise ValueError("Initial approximation 'p0' must be provided.")

    if df is None:
        df = differentiate(f)

    # Header
    print(f"{'n':>3} | {'p_n':>15} | {'f(p_n)':>13} | {'|p_n - p_{n-1}|':>18}")
    print("-" * 60)

    for n in range(1, max_iter + 1):
        fp0 = f(p0)
        dfp0 = df(p0)

        if dfp0 == 0:
            raise ValueError(f"Derivative is zero at p = {p0}. Method fails.")

        # Newton-Raphson update
        p = p0 - fp0 / dfp0
        fp = f(p)
        error = abs(p - p0)

        print(f"{n:>3} | {p:>15.8f} | {fp:>13.6e} | {error:>18.6e}")

        # Check stopping criteria
        if fp == 0 or error < tol:
            print("-" * 60)
            print(f"Converged to root p ~ {p:.8f} in {n} iterations.")
            return p

        p0 = p

    print("-" * 60)
    print("Maximum iterations reached.")
    return p


# --- Examples from the book ---
if __name__ == "__main__":
    f = lambda x: sin(cos(exp(x)))
    newton_raphson(f, p0=0.5, tol=1e-5)

