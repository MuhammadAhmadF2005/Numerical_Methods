import math
from math import cos, sin, tan, sinh, cosh, tanh, asin, acos, atan, asinh, acosh, atanh
from math import exp, log, log10, log2, e, pi


def newton_raphson(f, df, p0, tol=1e-5, max_iter=100):
    """
    Implements the Newton-Raphson Method from Burden & Faires.
    Requires f and its derivative df.
    Stops when |p - p0| < tol or f(p) == 0.
    """
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
    # f(x) = cos(x) - x, f'(x) = -sin(x) - 1
    f = lambda x: cos(x) - x
    df = lambda x: -sin(x) - 1
    newton_raphson(f, df, p0=pi / 4, tol=1e-5)
