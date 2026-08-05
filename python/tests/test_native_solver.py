"""Tests for the optional native (C++) SEUDO accelerator (seudo/_native).

All native-dependent tests skip gracefully if the extension hasn't been
built (see seudo/_native/build_native.sh) -- it's an optional accelerator,
not a hard dependency of the package.
"""

import numpy as np
import pytest

from seudo._native import NATIVE_AVAILABLE, seudo_solve
from seudo.estimate import estimate_time_courses_with_seudo

from test_estimate_synthetic import COMMON_KWARGS, make_synthetic_movie

requires_native = pytest.mark.skipif(not NATIVE_AVAILABLE, reason='native extension not built')


@requires_native
def test_native_solve_recovers_known_coefficient():
    # a single ROI, no noise/contamination -- solver should recover its
    # exact scale and leave the blob weights at 0
    wd = ht = 5
    roi = np.zeros((ht, wd))
    roi[2, 2] = 1.0
    roi[1, 2] = 0.5
    roi[2, 1] = 0.5
    image = roi * 2.0

    rois = roi.reshape(-1, 1)
    blob = np.array([[1.0]])
    n_weights = 1 + wd * ht
    lam = np.zeros(n_weights)
    lam[1:] = 0.5

    weights, n_steps, log = seudo_solve(
        image, blob, 1.0, rois, np.array([]), lam,
        0.001, 1000, 2, 0, False, 1,
    )
    assert weights[0] == pytest.approx(2.0, abs=0.01)
    assert np.max(weights[1:]) == pytest.approx(0.0, abs=1e-9)
    assert n_steps > 0


@requires_native
def test_native_and_python_solvers_agree_on_contamination_rejection():
    # same behavioral test as test_estimate_synthetic.py, run through the
    # native path -- should show the same qualitative (and roughly
    # quantitative) contamination-rejection behavior as the Python solver
    movie, profiles, cell_activity, contaminant_activity = make_synthetic_movie()

    result_py = estimate_time_courses_with_seudo(movie, profiles, use_native=False, **COMMON_KWARGS)
    result_native = estimate_time_courses_with_seudo(movie, profiles, use_native=True, **COMMON_KWARGS)

    tc_py = result_py['tc'][:, 0]
    tc_native = result_native['tc'][:, 0]

    assert not np.any(np.isnan(tc_native))
    # not bit-identical (different FISTA implementations/stopping paths),
    # but should be close on this well-separated, low-noise problem
    np.testing.assert_allclose(tc_native, tc_py, atol=0.05)

    contaminant_only = slice(15, 20)
    assert np.mean(np.abs(tc_native[contaminant_only] - cell_activity[contaminant_only])) < 1.0
    cell_only = slice(5, 10)
    assert np.mean(np.abs(tc_native[cell_only] - cell_activity[cell_only])) < 1.0


@requires_native
def test_n_jobs_gives_bit_identical_results_to_sequential():
    # pure parallelism (frame-level thread pool) must never change the
    # result, only wall-clock time -- see estimate.py's Phase 1/Phase 2 split
    movie, profiles, _, _ = make_synthetic_movie()

    result_seq = estimate_time_courses_with_seudo(movie, profiles, use_native=True, n_jobs=1, **COMMON_KWARGS)
    result_par = estimate_time_courses_with_seudo(movie, profiles, use_native=True, n_jobs=4, **COMMON_KWARGS)

    np.testing.assert_array_equal(result_seq['tc'], result_par['tc'])
    np.testing.assert_array_equal(result_seq['tc_lsq'], result_par['tc_lsq'])


def test_n_jobs_with_python_solver_also_bit_identical():
    # threading with the pure-Python solver doesn't help much (GIL-bound),
    # but it must still be correct
    movie, profiles, _, _ = make_synthetic_movie()

    result_seq = estimate_time_courses_with_seudo(movie, profiles, use_native=False, n_jobs=1, **COMMON_KWARGS)
    result_par = estimate_time_courses_with_seudo(movie, profiles, use_native=False, n_jobs=4, **COMMON_KWARGS)

    np.testing.assert_array_equal(result_seq['tc'], result_par['tc'])


def test_use_native_true_raises_clearly_when_unavailable(monkeypatch):
    import seudo.estimate as estimate_mod
    monkeypatch.setattr(estimate_mod, '_NATIVE_AVAILABLE', False)

    movie, profiles, _, _ = make_synthetic_movie()
    with pytest.raises(RuntimeError, match='not available'):
        estimate_time_courses_with_seudo(movie, profiles, use_native=True, **COMMON_KWARGS)


def test_use_native_auto_resolves_without_error():
    movie, profiles, _, _ = make_synthetic_movie()
    # should not raise regardless of whether the extension is built
    result = estimate_time_courses_with_seudo(movie, profiles, use_native='auto', **COMMON_KWARGS)
    assert result['params']['use_native'] == NATIVE_AVAILABLE
