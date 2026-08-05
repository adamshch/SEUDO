import numpy as np

from seudo.solver import fista_nonneg_weighted_l1


def test_recovers_sparse_nonneg_signal():
    rng = np.random.default_rng(0)
    n, m = 40, 15
    A_mat = rng.normal(size=(n, m))

    x_true = np.zeros(m)
    x_true[[2, 7, 11]] = [3.0, 1.5, 2.0]
    b = A_mat @ x_true

    A = lambda z: A_mat @ z
    At = lambda v: A_mat.T @ v

    lam = np.full(m, 0.01)
    x0 = np.zeros(m)
    x_hat = fista_nonneg_weighted_l1(A, At, b, lam, x0, tol=1e-6, max_iter=2000)

    assert np.all(x_hat >= -1e-8)
    assert np.linalg.norm(x_hat - x_true) / np.linalg.norm(x_true) < 0.05


def test_nonnegativity_enforced_on_negative_target():
    rng = np.random.default_rng(1)
    n, m = 20, 5
    A_mat = rng.normal(size=(n, m))
    b = -np.abs(rng.normal(size=n))  # a target that pulls toward negative x

    A = lambda z: A_mat @ z
    At = lambda v: A_mat.T @ v
    lam = np.zeros(m)
    x0 = np.zeros(m)

    x_hat = fista_nonneg_weighted_l1(A, At, b, lam, x0, tol=1e-8, max_iter=2000)
    assert np.all(x_hat >= -1e-8)
