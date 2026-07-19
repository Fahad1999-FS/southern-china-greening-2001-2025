"""
trend_stats.py
Non-parametric and autocorrelation-aware trend statistics used throughout the study.

Implements:
  * Mann-Kendall test (S, Var(S), Z, tau, p)  -- Mann (1945), Kendall (1975)
  * Theil-Sen slope estimator                 -- Theil (1950), Sen (1968)
  * GLS with AR(1) errors (Prais-Winsten)     -- Hamed & Rao (1998) rationale

Dependencies: numpy only (Student-t p-values via an incomplete beta function,
so scipy is not required).

Author: Kaleem Mehmood
License: MIT
"""

from __future__ import annotations
import math
import numpy as np

__all__ = ["mann_kendall", "theil_sen", "gls_ar1", "zscore"]


# --------------------------------------------------------------------------
# Student-t survival function without scipy
# --------------------------------------------------------------------------
def _betacf(a: float, b: float, x: float, maxit: int = 300, eps: float = 3e-16) -> float:
    fpmin = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, maxit + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        c = 1.0 + aa / c
        d = fpmin if abs(d) < fpmin else d
        c = fpmin if abs(c) < fpmin else c
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        c = 1.0 + aa / c
        d = fpmin if abs(d) < fpmin else d
        c = fpmin if abs(c) < fpmin else c
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < eps:
            break
    return h


def _betai(a: float, b: float, x: float) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lb = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
          + a * math.log(x) + b * math.log(1.0 - x))
    bt = math.exp(lb)
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def _t_two_sided(t: float, dof: int) -> float:
    """Two-sided p-value for a Student-t statistic."""
    return _betai(dof / 2.0, 0.5, dof / (dof + t * t))


def _norm_two_sided(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0))


# --------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------
def zscore(x) -> np.ndarray:
    """Standardise to zero mean and unit variance (sample sd)."""
    x = np.asarray(x, dtype=float)
    return (x - x.mean()) / x.std(ddof=1)


def mann_kendall(y) -> dict:
    """
    Mann-Kendall trend test with tie correction.

    Returns dict with S, var_S, Z, tau, p (two-sided).
    """
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 3:
        raise ValueError("Mann-Kendall requires at least 3 observations")

    s = 0
    for i in range(n - 1):
        s += int(np.sign(y[i + 1:] - y[i]).sum())

    # tie-corrected variance
    _, counts = np.unique(y, return_counts=True)
    tie_term = sum(c * (c - 1) * (2 * c + 5) for c in counts if c > 1)
    var_s = (n * (n - 1) * (2 * n + 5) - tie_term) / 18.0

    if s > 0:
        z = (s - 1) / math.sqrt(var_s)
    elif s < 0:
        z = (s + 1) / math.sqrt(var_s)
    else:
        z = 0.0

    return {"S": s, "var_S": var_s, "Z": z,
            "tau": s / (0.5 * n * (n - 1)), "p": _norm_two_sided(z)}


def theil_sen(y, x=None) -> float:
    """Theil-Sen slope: median of all pairwise slopes."""
    y = np.asarray(y, dtype=float)
    x = np.arange(len(y), dtype=float) if x is None else np.asarray(x, dtype=float)
    slopes = [(y[j] - y[i]) / (x[j] - x[i])
              for i in range(len(y) - 1) for j in range(i + 1, len(y))
              if x[j] != x[i]]
    return float(np.median(slopes))


def gls_ar1(y, X, max_iter: int = 400, tol: float = 1e-11) -> dict:
    """
    Generalised least squares with AR(1) errors, fitted by Prais-Winsten
    iteration.

    Parameters
    ----------
    y : (n,) response vector
    X : (n, k) design matrix; include an intercept column explicitly.

    Returns
    -------
    dict with beta, se, p, phi, aic, dof
    """
    y = np.asarray(y, dtype=float)
    X = np.asarray(X, dtype=float)
    n, k = X.shape

    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    phi = 0.0

    for _ in range(max_iter):
        resid = y - X @ beta
        phi_new = float(np.sum(resid[1:] * resid[:-1]) / np.sum(resid[:-1] ** 2))
        phi_new = max(-0.99, min(0.99, phi_new))
        yt, Xt = _pw_transform(y, X, phi_new)
        beta_new = np.linalg.lstsq(Xt, yt, rcond=None)[0]
        converged = (np.max(np.abs(beta_new - beta)) < tol
                     and abs(phi_new - phi) < tol)
        beta, phi = beta_new, phi_new
        if converged:
            break

    yt, Xt = _pw_transform(y, X, phi)
    resid = yt - Xt @ beta
    dof = n - k
    sigma2 = float(resid @ resid) / dof
    cov = sigma2 * np.linalg.inv(Xt.T @ Xt)
    se = np.sqrt(np.diag(cov))
    p = np.array([_t_two_sided(b / s, dof) for b, s in zip(beta, se)])

    rss = float(resid @ resid)
    loglik = -0.5 * n * (math.log(2 * math.pi * rss / n) + 1.0)
    aic = 2 * (k + 1) - 2 * loglik

    return {"beta": beta, "se": se, "p": p, "phi": phi, "aic": aic, "dof": dof}


def _pw_transform(y, X, phi):
    """Prais-Winsten transformation (retains the first observation)."""
    n = len(y)
    s = math.sqrt(1.0 - phi * phi)
    yt = np.empty(n)
    Xt = np.empty_like(X)
    yt[0] = s * y[0]
    Xt[0] = s * X[0]
    yt[1:] = y[1:] - phi * y[:-1]
    Xt[1:] = X[1:] - phi * X[:-1]
    return yt, Xt


if __name__ == "__main__":
    rng = np.random.default_rng(42)
    t = np.arange(25, dtype=float)
    y = 0.03 * t + rng.normal(0, 0.05, 25)
    mk = mann_kendall(y)
    print(f"Theil-Sen slope : {theil_sen(y, t):.4f}")
    print(f"Mann-Kendall    : tau={mk['tau']:.3f}  p={mk['p']:.3g}")
    Xd = np.column_stack([np.ones(25), zscore(t)])
    fit = gls_ar1(zscore(y), Xd)
    print(f"GLS-AR(1) slope : {fit['beta'][1]:.3f} "
          f"(se {fit['se'][1]:.3f}, p {fit['p'][1]:.3g}, phi {fit['phi']:.3f})")
