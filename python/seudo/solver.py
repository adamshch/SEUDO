"""FISTA solver for the nonnegative weighted-L1 problem TFOCS's
solver_L1RLS (with the 'nonneg' option) solves in the MATLAB code:

    minimize_{x >= 0}  0.5 * ||A(x) - b||^2 + lambda' * x

Since x is constrained nonnegative, the L1 penalty lambda' * |x| is linear
on the feasible set, so the proximal operator of the penalty + indicator is
just a shifted ReLU -- see extras/prox_l1_semipos.m in the MATLAB code
(the commented-out `prox_f`: x = max(0, x - t*q)).
"""

import numpy as np


def fista_nonneg_weighted_l1(A, At, b, lam, x0, tol=0.01, max_iter=1000, l0=1.0):
    """Minimize 0.5*||A(x) - b||^2 + lam^T x over x >= 0.

    A, At: callables implementing the forward and adjoint linear operator.
    lam: per-coordinate L1 weight (elementwise, same shape as x), >= 0.
    x0: initial point.
    tol: relative objective-change stopping tolerance.

    Returns the solution vector x.
    """
    x = np.array(x0, dtype=float, copy=True)
    y = x.copy()
    t = 1.0
    L = l0

    Ax = A(x)
    f_prev = 0.5 * np.dot(Ax - b, Ax - b) + np.dot(lam, x)

    for _ in range(max_iter):
        Ay = A(y)
        resid = Ay - b
        grad = At(resid)
        fy = 0.5 * np.dot(resid, resid)

        while True:
            x_new = np.maximum(0.0, y - (grad + lam) / L)
            diff = x_new - y
            Ax_new = A(x_new)
            lhs = 0.5 * np.dot(Ax_new - b, Ax_new - b)
            rhs = fy + np.dot(grad, diff) + 0.5 * L * np.dot(diff, diff)
            if lhs <= rhs + 1e-12 or L > 1e10:
                break
            L *= 2.0

        t_new = 0.5 * (1.0 + np.sqrt(1.0 + 4.0 * t * t))
        y = x_new + ((t - 1.0) / t_new) * (x_new - x)

        f_new = lhs + np.dot(lam, x_new)
        rel_change = abs(f_new - f_prev) / max(1.0, abs(f_prev))

        x, t, f_prev = x_new, t_new, f_new

        if rel_change < tol:
            break

    return x
