# Numerical Methods in Python

A structured repository implementing standard numerical analysis algorithms in Python, following the formulations from **Burden & Faires (*Numerical Analysis*)**.

---

## 📌 Features

- **Standard Implementations**: Clear, textbook-aligned numerical methods with robust stopping criteria.
- **Detailed Iteration Logs**: Formatted console output showing iteration counts ($n$), endpoints ($a_n, b_n$), approximations ($p_n$), function evaluations ($f(p_n)$), and error bounds.
- **Minimal Dependencies**: Built purely with standard Python mathematical libraries.

---

## 🗂️ Repository Structure

```text
.
├── README.md
└── bisection_solver.py     # Implementation of the Bisection Method
```

---

## 🚀 Implemented Methods

### 1. Bisection Method (`bisection_solver.py`)

The **Bisection Method** is a bracketed root-finding method based on the Intermediate Value Theorem. If a continuous function $f(x)$ satisfies $f(a) \cdot f(b) < 0$, a root must exist in $(a, b)$.

- **Formula**:
  $$p_n = a_n + \frac{b_n - a_n}{2}$$
- **Stopping Criteria**:
  $$\frac{b_n - a_n}{2} < \text{TOL} \quad \text{or} \quad f(p_n) = 0$$

#### Usage Example

```python
from bisection_solver import bisection_method
import math

# Define your function
f = lambda x: math.exp(x) - (x**2) + (3 * x) - 2

# Solve on interval [0, 1] with tolerance 1e-5
root = bisection_method(f, a=0.0, b=1.0, tol=1e-5)
print(f"Calculated root: {root}")
```

#### Sample Output

```text
  n |         a_n |         b_n |         p_n |        f(p_n) |     (b-a)/2
---------------------------------------------------------------------------
  1 |    0.000000 |    1.000000 |    0.500000 |  8.987213e-01 | 5.000000e-01
  2 |    0.000000 |    0.500000 |    0.250000 | -2.847458e-02 | 2.500000e-01
...
 17 |    0.257523 |    0.257538 |    0.257530 | -2.759847e-07 | 7.629395e-06
---------------------------------------------------------------------------
Converged to root p ~ 0.257530 in 17 iterations.
```

---

## 🛠️ Getting Started

### Prerequisites

- Python 3.8+ (no external dependencies required)

### Running the Solvers

Run directly using Python:

```bash
python bisection_solver.py
```

---

## 🗺️ Roadmap / Upcoming Methods

- [ ] Fixed-Point Iteration
- [ ] Newton-Raphson Method
- [ ] Secant Method & Method of False Position
- [ ] Gaussian Elimination & LU Decomposition
- [ ] Lagrange & Newton Divided Difference Interpolation
- [ ] Numerical Integration (Trapezoidal, Simpson's Rules)
- [ ] Runge-Kutta Methods for ODEs

---

## 📚 References

- Burden, R. L., & Faires, J. D. *Numerical Analysis* (9th / 10th Edition). Cengage Learning.
