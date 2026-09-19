import math
from math import cos, sin, tan, sinh, cosh, tanh, asin, acos, atan, asinh, acosh, atanh
from math import exp, log, log10, log2, e, pi


def fixed_point_iteration(g, p0, tol=1e-5, max_iter=100):
    """
    Implements the Fixed-Point Iteration Method from Burden & Faires.
    Finds p such that g(p) = p (i.e., a fixed point of g).
    Stops when |p - p0| < tol.
    """
    # Header
    print(f"{'n':>3} | {'p_n':>15} | {'|p_n - p_{n-1}|':>18}")
    print("-" * 44)

    for n in range(1, max_iter + 1):
        # Fixed-point update
        p = g(p0)
        error = abs(p - p0)

        print(f"{n:>3} | {p:>15.8f} | {error:>18.6e}")

        # Check stopping criteria
        if error < tol:
            print("-" * 44)
            print(f"Converged to fixed point p ~ {p:.8f} in {n} iterations.")
            return p

        p0 = p

    print("-" * 44)
    print("Maximum iterations reached.")
    return p


# --- Examples from the book ---
if __name__ == "__main__":
    # Solve x = (2 - exp(x) + x^2) / 3, derived from f(x) = exp(x) - x^2 + 3x - 2 = 0
    g = lambda x: (2 - exp(x) + x**2) / 3
    fixed_point_iteration(g, p0=0.0, tol=1e-5)
