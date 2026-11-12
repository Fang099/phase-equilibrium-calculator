"""
numerics.py -- small root finders in plain Python and NumPy.

equilibrium_solver.py uses SciPy's brentq and fsolve when SciPy is installed.
These are drop-in replacements for the two calls it makes, used where SciPy
is not available: in particular in the browser, where the website runs these
same modules with only NumPy loaded. Set the environment variable
PEC_NO_SCIPY=1 to force their use (the test suite passes either way).
"""

import numpy as np


def brentq(f, a, b, xtol=2e-12, rtol=8.881784197001252e-16, maxiter=100, args=()):
    """Root of f in [a, b] by Brent's method. f(a) and f(b) must differ in sign."""
    fa, fb = f(a, *args), f(b, *args)
    if fa == 0:
        return a
    if fb == 0:
        return b
    if fa * fb > 0:
        raise ValueError("f(a) and f(b) must have different signs")
    c, fc = a, fa
    d = e = b - a
    for _ in range(maxiter):
        if fb * fc > 0:
            c, fc = a, fa
            d = e = b - a
        if abs(fc) < abs(fb):
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb
        tol1 = 2 * rtol * abs(b) + 0.5 * xtol
        xm = 0.5 * (c - b)
        if abs(xm) <= tol1 or fb == 0:
            return b
        if abs(e) >= tol1 and abs(fa) > abs(fb):
            s = fb / fa
            if a == c:                       # secant step
                p = 2 * xm * s
                q = 1 - s
            else:                            # inverse quadratic interpolation
                q, r = fa / fc, fb / fc
                p = s * (2 * xm * q * (q - r) - (b - a) * (r - 1))
                q = (q - 1) * (r - 1) * (s - 1)
            if p > 0:
                q = -q
            p = abs(p)
            if 2 * p < min(3 * xm * q - abs(tol1 * q), abs(e * q)):
                e, d = d, p / q
            else:
                d = e = xm
        else:
            d = e = xm
        a, fa = b, fb
        b += d if abs(d) > tol1 else (tol1 if xm > 0 else -tol1)
        fb = f(b, *args)
    raise RuntimeError("brentq did not converge")


def fsolve(func, x0, full_output=False, xtol=1.49012e-08, maxiter=80):
    """Solve func(x) = 0 by Newton's method with a finite-difference Jacobian
    and step halving. With full_output=True returns (x, {'fvec': f}, ier, msg)
    like scipy.optimize.fsolve; ier == 1 means converged."""
    x = np.array(x0, dtype=float).ravel()
    f = np.asarray(func(x), dtype=float).ravel()
    n = len(x)
    ier, msg = 5, "no convergence"
    for _ in range(maxiter):
        fn = np.linalg.norm(f)
        if fn < 1e-13:
            ier, msg = 1, "converged"
            break
        J = np.empty((len(f), n))
        for j in range(n):
            h = 1e-7 * max(1.0, abs(x[j]))
            xp = x.copy()
            xp[j] += h
            J[:, j] = (np.asarray(func(xp), dtype=float).ravel() - f) / h
        try:
            dx = np.linalg.solve(J, -f)
        except np.linalg.LinAlgError:
            dx = np.linalg.lstsq(J, -f, rcond=None)[0]
        lam, accepted = 1.0, False
        while lam > 1e-6:
            xn = x + lam * dx
            f_new = np.asarray(func(xn), dtype=float).ravel()
            if np.all(np.isfinite(f_new)) and np.linalg.norm(f_new) < fn:
                accepted = True
                break
            lam *= 0.5
        if not accepted:
            msg = "line search failed"
            break
        step = np.linalg.norm(xn - x)
        x, f = xn, f_new
        if step <= xtol * (xtol + np.linalg.norm(x)):
            ier, msg = 1, "converged"
            break
    if full_output:
        return x, {"fvec": f}, ier, msg
    return x
