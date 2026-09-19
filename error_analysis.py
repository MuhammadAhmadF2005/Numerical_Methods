import math
from math import exp, log, log10, log2, e, pi


def absolute_error(true_val, approx_val):
    """Absolute error: |true - approx|"""
    return abs(true_val - approx_val)


def relative_error(true_val, approx_val):
    """Relative error: |true - approx| / |true|"""
    if true_val == 0:
        raise ValueError("True value is zero; relative error is undefined.")
    return abs(true_val - approx_val) / abs(true_val)


def percentage_relative_error(true_val, approx_val):
    """Percentage relative error: (|true - approx| / |true|) * 100"""
    return relative_error(true_val, approx_val) * 100


def significant_digits(approx_val, prev_approx_val):
    """
    Estimates significant digits using the stopping criterion from Burden & Faires:
        |p_n - p_{n-1}| / |p_n| < 0.5 * 10^(2-d)
    Returns approximate number of significant digits d.
    """
    if approx_val == 0:
        return float('inf')
    eps = abs(approx_val - prev_approx_val) / abs(approx_val)
    if eps == 0:
        return float('inf')
    return 2 - log10(2 * eps)


def truncation_error_analysis(true_val, approx_series, labels=None):
    """
    Prints a table of absolute and relative errors for a sequence of approximations.
    Useful for Taylor series truncation error analysis.
    """
    if labels is None:
        labels = [f"Term {i+1}" for i in range(len(approx_series))]

    print(f"{'Approximation':<20} | {'Value':>18} | {'Abs. Error':>14} | {'Rel. Error (%)':>16}")
    print("-" * 76)

    for label, approx in zip(labels, approx_series):
        abs_err = absolute_error(true_val, approx)
        pct_err = percentage_relative_error(true_val, approx)
        print(f"{label:<20} | {approx:>18.10f} | {abs_err:>14.6e} | {pct_err:>15.6f}%")

    print("-" * 76)
    print(f"{'True value':<20} | {true_val:>18.10f}")


# --- Examples from the book ---
if __name__ == "__main__":
    # Taylor series approximation of e^1 using partial sums
    import math

    true_value = math.e  # e = 2.718281828...
    terms = []
    partial_sum = 0.0
    labels = []

    print("Taylor Series Approximation of e^1 = sum(1/n!) for n = 0, 1, 2, ...\n")

    for k in range(9):
        partial_sum += 1.0 / math.factorial(k)
        terms.append(partial_sum)
        labels.append(f"n = {k}  (k+1 terms)")

    truncation_error_analysis(true_value, terms, labels)
