# Numerical Methods

Standard numerical analysis algorithms implemented in Python, based on Burden & Faires (*Numerical Analysis*).

## Implemented Methods

### Bisection Method (`bisection_solver.py`)

Finds a root $p \in (a, b)$ of a continuous function $f(x)$ where $f(a) \cdot f(b) < 0$.

- **Iteration Formula**:
  $$p_n = a_n + \frac{b_n - a_n}{2}$$
- **Stopping Criteria**:
  $$\frac{b_n - a_n}{2} < \text{TOL} \quad \text{or} \quad f(p_n) = 0$$

#### Usage

```python
from math import sin, pi
from bisection_solver import bisection_method

f = lambda x: x + 1 - 2 * sin(pi * x)
root = bisection_method(f, a=0.5, b=1.0, tol=1e-5)
```

#### Output

```text
  n |         a_n |         b_n |         p_n |        f(p_n) |     (b-a)/2
---------------------------------------------------------------------------
  1 |    0.500000 |    1.000000 |    0.750000 |  3.357864e-01 | 2.500000e-01
  2 |    0.500000 |    0.750000 |    0.625000 | -2.227591e-01 | 1.250000e-01
...
 16 |    0.681961 |    0.681976 |    0.681969 | -2.692392e-05 | 7.629395e-06
---------------------------------------------------------------------------
Converged to root p ~ 0.681969 in 16 iterations.
```

## Execution

```bash
python bisection_solver.py
```

## References

- Burden, R. L., & Faires, J. D. *Numerical Analysis*. Cengage Learning.
