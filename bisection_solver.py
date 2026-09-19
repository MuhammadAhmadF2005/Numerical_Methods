from math import exp
import math

def bisection_method(f, a, b, tol=1e-5, max_iter=100):
    """
    Implements the standard Bisection Method from Burden & Faires.
    Stops when (b - a) / 2 < tol or f(p) == 0.
    """
    fa = f(a)
    fb = f(b)

    if fa * fb > 0:
        raise ValueError(f"Function must have opposite signs at endpoints: f({a})={fa}, f({b})={fb}")

    # Header
    print(f"{'n':>3} | {'a_n':>11} | {'b_n':>11} | {'p_n':>11} | {'f(p_n)':>13} | {'(b-a)/2':>11}")
    print("-" * 75)

    for n in range(1, max_iter + 1):
        p = a + (b - a) / 2.0
        fp = f(p)
        error = (b - a) / 2.0

        print(f"{n:3d} | {a:11.6f} | {b:11.6f} | {p:11.6f} | {fp:13.6e} | {error:11.6e}")

        # Check stopping criteria
        if fp == 0 or error < tol:
            print("-" * 75)
            print(f"Converged to root p ~ {p:.6f} in {n} iterations.")
            return p

        # Interval reduction
        if fa * fp > 0:
            a = p
            fa = fp
        else:
            b = p
            fb = fp

    print("-" * 75)
    print("Maximum iterations reached.")
    return p


# --- Examples from the book ---
if __name__ == "__main__":
    f_5a = lambda x: exp(x)-(x**2)+(3*x)-2
    bisection_method(f_5a, 0.0, 1.0, tol=1e-5)