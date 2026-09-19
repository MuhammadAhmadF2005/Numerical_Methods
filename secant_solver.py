import math
from math import cos, sin, tan, sinh, cosh, tanh, asin, acos, atan, asinh, acosh, atanh
from math import exp, log, log10, log2, e, pi


def secant_method(f, p0, p1, tol=1e-5, max_iter=100):
    """
    Implements the Secant Method from Burden & Faires.
    Uses two initial approximations p0 and p1.
    Stops when |p - p1| < tol or f(p) == 0.
    """
    fp0 = f(p0)
    fp1 = f(p1)

    # Header
    print(f"{'n':>3} | {'p_n':>15} | {'f(p_n)':>13} | {'|p_n - p_{n-1}|':>18}")
    print("-" * 60)

    for n in range(2, max_iter + 2):
        if fp1 - fp0 == 0:
            raise ValueError("Division by zero: f(p1) - f(p0) = 0. Method fails.")

        # Secant update
        p = p1 - fp1 * (p1 - p0) / (fp1 - fp0)
        fp = f(p)
        error = abs(p - p1)

        print(f"{n:>3} | {p:>15.8f} | {fp:>13.6e} | {error:>18.6e}")

        # Check stopping criteria
        if fp == 0 or error < tol:
            print("-" * 60)
            print(f"Converged to root p ~ {p:.8f} in {n} iterations.")
            return p

        # Shift for next iteration
        p0, fp0 = p1, fp1
        p1, fp1 = p, fp

    print("-" * 60)
    print("Maximum iterations reached.")
    return p


# --- Examples from the book ---
if __name__ == "__main__":
    # f(x) = cos(x) - x
    f = lambda x: cos(x) - x
    secant_method(f, p0=0.5, p1=pi / 4, tol=1e-5)
