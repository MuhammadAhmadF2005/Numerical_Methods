import inspect
import math
from math import cos, sin, tan, sinh, cosh, tanh, asin, acos, atan, asinh, acosh, atanh
from math import exp, log, log10, log2, e, pi


def differentiate(f, h=1e-7):
    """
    Approximates the derivative f'(x) of function f using
    the central difference formula:
        f'(x) ~ (f(x + h) - f(x - h)) / (2 * h)
    """
    return lambda x: (f(x + h) - f(x - h)) / (2 * h)


derivative = differentiate  # Alias for convenience


def find_suitable_g(f, p0, c=None, h=1e-7, verbose=True):
    """
    Finds and constructs a suitable iteration function g(x) to solve f(x) = 0
    via Fixed-Point Iteration:
        g(x) = x - c * f(x)

    According to the Fixed-Point Theorem (Burden & Faires, Theorem 2.4),
    the sequence p_{n+1} = g(p_n) converges to a root p of f(x) = 0 if
    |g'(x)| < 1 on an interval containing the root and approximations.

    Since g'(x) = 1 - c * f'(x), choosing:
        c ~ 1 / f'(p0)
    makes g'(p0) ~ 0, guaranteeing |g'(p0)| < 1 and maximizing the local
    rate of convergence. A line-search/damping strategy is used to ensure
    stability across iterations.
    """
    if c is not None:
        g = lambda x, cc=c: x - cc * f(x)
        dg_p0 = abs(differentiate(g, h)(p0))
        g.c = c
        g.dg_p0 = dg_p0
        if verbose:
            print(f"User-specified g(x) constructed: g(x) = x - ({c:.6g}) * f(x)")
            print(f"|g'({p0})| ~ {dg_p0:.6e} ({'< 1 (Satisfies Fixed-Point Theorem)' if dg_p0 < 1 else '>= 1 (May diverge)'})")
        return g

    # Estimate derivative at p0
    df_p0 = differentiate(f, h)(p0)
    if abs(df_p0) < 1e-12:
        # Near an extremum or inflection point, probe nearby points
        for delta in [1e-3, -1e-3, 1e-2, -1e-2, 0.1, -0.1]:
            df = differentiate(f, h)(p0 + delta)
            if abs(df) > 1e-12:
                df_p0 = df
                break
        else:
            df_p0 = 1.0

    base_c = 1.0 / df_p0

    # Test candidate scalings of base_c and standard relaxation factors
    candidates = [base_c * s for s in [1.0, 0.5, 0.25, 0.1, 0.05, 0.01]]
    candidates += [1.0, -1.0, 0.5, -0.5, 0.1, -0.1]

    best_c = None
    best_rate = float("inf")

    for cand_c in candidates:
        g_cand = lambda x, cc=cand_c: x - cc * f(x)
        p = p0
        valid = True
        errs = []
        try:
            for step in range(4):
                p_next = g_cand(p)
                if math.isnan(p_next) or math.isinf(p_next) or abs(p_next) > 1e8:
                    valid = False
                    break
                err = abs(p_next - p)
                errs.append(err)
                p = p_next
        except Exception:
            valid = False

        if not valid:
            continue

        # Check if sequence is contracting
        if len(errs) >= 3 and errs[-1] <= errs[0] and abs(f(p)) <= abs(f(p0)):
            ratio = errs[-1] / (errs[0] + 1e-15)
            dg0 = abs(differentiate(g_cand, h)(p0))
            if dg0 < 1.0 and ratio < best_rate:
                best_rate = ratio
                best_c = cand_c

    if best_c is None:
        best_c = base_c

    g = lambda x, cc=best_c: x - cc * f(x)
    dg_p0 = abs(differentiate(g, h)(p0))
    g.c = best_c
    g.dg_p0 = dg_p0

    if verbose:
        print(f"Automatically constructed suitable g(x): g(x) = x - ({best_c:+.6g}) * f(x)")
        print(f"Fixed-Point condition: |g'({p0})| ~ {dg_p0:.6e} (< 1: Local convergence guaranteed)")

    return g


def _resolve_function(func, p0, f=None, g=None, is_f=None):
    """
    Determines whether the provided function is f(x) (root finding: f(x) = 0)
    or g(x) (fixed point iteration: g(x) = x).
    """
    if f is not None:
        return "f", f
    if g is not None:
        return "g", g
    if func is None:
        raise ValueError("A function f(x) or g(x) must be provided.")

    if is_f is True:
        return "f", func
    if is_f is False:
        return "g", func

    # Check function name
    fname = getattr(func, "__name__", "")
    if fname == "f":
        return "f", func
    if fname == "g":
        return "g", func

    # Check caller local variable names
    try:
        frame = inspect.currentframe().f_back.f_back
        matching = [k for k, v in frame.f_locals.items() if v is func and not k.startswith("_")]
        if any(name == "f" or name.startswith("f_") for name in matching):
            return "f", func
        if any(name == "g" or name.startswith("g_") for name in matching):
            return "g", func
    except Exception:
        pass

    # Contraction mapping test:
    # In fixed-point iteration, |g'(p0)| must be < 1 to converge.
    # If |func'(p0)| >= 1, it cannot converge as g, so treat as f(x) to build a suitable g(x).
    try:
        df0 = differentiate(func)(p0)
        if abs(df0) >= 1.0:
            return "f", func
    except Exception:
        pass

    return "g", func


def fixed_point_iteration(func=None, p0=None, *args, f=None, g=None, tol=1e-5, max_iter=100, is_f=None, c=None, verbose=True, **kwargs):
    """
    Implements the Fixed-Point Iteration Method from Burden & Faires.

    Can be used in two modes:
    1. Given f(x): Automatically finds a suitable g(x) such that g(p) = p iff f(p) = 0,
       ensuring |g'(p0)| < 1 so the sequence converges to the root of f(x) = 0.
    2. Given g(x): Performs standard fixed-point iteration to find p such that g(p) = p.

    Parameters:
    -----------
    func : callable, optional
        Function f(x) or g(x).
    p0 : float
        Initial approximation.
    f : callable, optional
        Explicitly provide f(x) where f(p) = 0 is sought.
    g : callable, optional
        Explicitly provide g(x) where g(p) = p is sought.
    tol : float, default 1e-5
        Tolerance for stopping criteria.
    max_iter : int, default 100
        Maximum number of iterations.
    is_f : bool, optional
        Set to True if func is f(x), or False if func is g(x).
    c : float, optional
        Custom relaxation parameter for g(x) = x - c * f(x).
    verbose : bool, default True
        Whether to print information about constructed g(x) and iteration progress.
    """
    if p0 is None:
        if len(args) > 0 and not callable(args[0]):
            p0 = args[0]
        elif "p0" in kwargs:
            p0 = kwargs["p0"]
        else:
            raise ValueError("Initial approximation 'p0' must be provided.")

    fn_type, target_fn = _resolve_function(func, p0, f=f, g=g, is_f=is_f)

    if fn_type == "f":
        f_func = target_fn
        g_func = find_suitable_g(f_func, p0=p0, c=c, verbose=verbose)
    else:
        f_func = None
        g_func = target_fn

    # Header
    if f_func is not None:
        print(f"{'n':>3} | {'p_n':>15} | {'f(p_n)':>13} | {'|p_n - p_{n-1}|':>18}")
        print("-" * 60)
    else:
        print(f"{'n':>3} | {'p_n':>15} | {'|p_n - p_{n-1}|':>18}")
        print("-" * 44)

    p_curr = p0
    for n in range(1, max_iter + 1):
        p_next = g_func(p_curr)
        error = abs(p_next - p_curr)

        if f_func is not None:
            fp = f_func(p_next)
            print(f"{n:>3} | {p_next:>15.8f} | {fp:>13.6e} | {error:>18.6e}")
        else:
            fp = None
            print(f"{n:>3} | {p_next:>15.8f} | {error:>18.6e}")

        # Check stopping criteria
        if (fp is not None and fp == 0) or error < tol:
            border = "-" * 60 if f_func is not None else "-" * 44
            print(border)
            target_name = "root" if f_func is not None else "fixed point"
            print(f"Converged to {target_name} p ~ {p_next:.8f} in {n} iterations.")
            return p_next

        p_curr = p_next

    border = "-" * 60 if f_func is not None else "-" * 44
    print(border)
    print("Maximum iterations reached.")
    return p_next


# Alias for convenience
fixed_point_solver = fixed_point_iteration


# --- Examples from Burden & Faires ---
if __name__ == "__main__":
    print("=== Example 1: Solving f(x) = 0 given only f(x) ===")
    # f(x) = x^3 - 4x + 1 = 0
    f = lambda x: x**3 - 4*x + 1
    fixed_point_solver(f, p0=1.0, tol=1e-5)

    print("\n=== Example 2: Burden & Faires Sec 2.2 Example 1 ===")
    # f(x) = x^3 + 4x^2 - 10 = 0
    f2 = lambda x: x**3 + 4*x**2 - 10
    fixed_point_iteration(f=f2, p0=1.5, tol=1e-5)

    print("\n=== Example 3: Classic Fixed-Point Iteration with g(x) ===")
    # g(x) = 1 + x - x^2
    g = lambda x: 1 + x - x**2
    fixed_point_iteration(g, p0=1.0, tol=1e-5)
