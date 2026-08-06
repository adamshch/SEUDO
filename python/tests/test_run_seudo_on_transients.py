import numpy as np
import pytest

from seudo.core import Seudo
from seudo.run_seudo_on_transients import frame_blocks_for_cell, run_seudo_restricted_to_transients


def gaussian_patch(shape, center, sigma):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    g = np.exp(-(((yy - center[0]) ** 2 + (xx - center[1]) ** 2) / (2 * sigma ** 2)))
    g[g < 0.05 * g.max()] = 0.0
    return g


def make_two_cell_seudo():
    y, x, f = 20, 40, 600
    prof0 = gaussian_patch((y, x), center=(10, 10), sigma=2.0)
    prof1 = gaussian_patch((y, x), center=(10, 30), sigma=2.0)

    act0 = np.zeros(f)
    act0[100:110] = 3.0
    act1 = np.zeros(f)
    act1[300:310] = 3.0
    act1[450:460] = 3.0

    rng = np.random.default_rng(7)
    movie = np.zeros((y, x, f))
    for ff in range(f):
        movie[:, :, ff] = prof0 * act0[ff] + prof1 * act1[ff] + 0.01 * rng.normal(size=(y, x))

    profiles = np.stack([prof0, prof1], axis=2)
    tc = np.stack([act0, act1], axis=1)

    se = Seudo(movie, profiles, time_courses=tc)
    se.compute_transient_info('default', min_duration=3, blur_radius=1)
    return se


@pytest.fixture
def se():
    return make_two_cell_seudo()


COMMON = dict(p=1e-4, sigma2=0.02, lambda_blob=1.0, blob_radius=1.2, pad_space=6, solver_max_iter=200)


def test_frame_blocks_for_cell_matches_transient_info(se):
    ti = se.tc_default['transient_info']
    blocks0 = frame_blocks_for_cell(ti, 0)
    blocks1 = frame_blocks_for_cell(ti, 1)

    assert len(blocks0) == 1
    assert len(blocks1) == 2
    assert blocks0 == [(int(s), int(e)) for s, e in ti[0]['times']]


def test_run_one_cell_only_computes_that_cells_frames(se):
    result = run_seudo_restricted_to_transients(se, 'default', which_cells=[0], verbose=False, **COMMON)

    assert result['tc'].shape == (se.mov_f, 1)
    covered = ~np.isnan(result['tc'][:, 0])
    ti0 = se.tc_default['transient_info'][0]
    expected_covered = np.zeros(se.mov_f, dtype=bool)
    for s, e in ti0['times']:
        expected_covered[s:e + 1] = True

    assert np.array_equal(covered, expected_covered)
    # most of the movie was never touched -> still NaN
    assert covered.sum() < se.mov_f * 0.5

    assert se.tc_seudo[-1] is result


def test_run_all_cells_gives_each_its_own_columns(se):
    result = run_seudo_restricted_to_transients(se, 'default', verbose=False, **COMMON)

    assert result['tc'].shape == (se.mov_f, 2)
    ti = se.tc_default['transient_info']

    for col, cell_idx in enumerate(result['params']['which_cells']):
        covered = ~np.isnan(result['tc'][:, col])
        expected = np.zeros(se.mov_f, dtype=bool)
        for s, e in ti[cell_idx]['times']:
            expected[s:e + 1] = True
        assert np.array_equal(covered, expected)


def test_cell_with_no_transients_is_skipped_gracefully():
    y, x, f = 20, 20, 100
    prof = gaussian_patch((y, x), (10, 10), 2.0)
    movie = np.zeros((y, x, f))
    tc = np.zeros((f, 1))
    se = Seudo(movie, prof[:, :, None], time_courses=tc)
    se.compute_transient_info('default', transient_frames=np.zeros((f, 1), dtype=bool))

    result = run_seudo_restricted_to_transients(se, 'default', verbose=False, **COMMON)
    assert np.all(np.isnan(result['tc']))


def test_progress_callback_invoked_per_cell(se):
    calls = []
    run_seudo_restricted_to_transients(
        se, 'default', verbose=False,
        progress_callback=lambda done, total, cell_idx: calls.append((done, total, cell_idx)),
        **COMMON,
    )
    assert calls == [(1, 2, 0), (2, 2, 1)]


def test_frame_progress_callback_invoked_per_frame(se):
    # cell 0 has a single, multi-frame transient block -- progress_callback
    # alone would report nothing until the *whole* cell finishes, which is
    # what made the GUI's single-cell progress bar look stuck at 0% then
    # jump straight to 100%. frame_progress_callback should report smooth,
    # increasing within-cell progress instead.
    calls = []
    run_seudo_restricted_to_transients(
        se, 'default', which_cells=[0], verbose=False,
        frame_progress_callback=lambda col, cell_idx, done, total: calls.append((col, cell_idx, done, total)),
        **COMMON,
    )
    assert len(calls) > 1
    assert all(col == 0 and cell_idx == 0 for col, cell_idx, done, total in calls)
    assert [c[2] for c in calls] == list(range(1, len(calls) + 1))
    total_frames = calls[0][3]
    assert all(c[3] == total_frames for c in calls)
    assert calls[-1][2] == total_frames
