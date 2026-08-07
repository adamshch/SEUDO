"""End-to-end sanity test for estimate_time_courses_with_seudo on synthetic
data, since we can't run the real MATLAB code headless in this environment
to generate reference numbers (see conversation / commit notes).

This does not check numeric parity with MATLAB -- it checks the behavior
SEUDO is meant to provide: given a movie containing one "found" cell (whose
profile is supplied) and one spatially-overlapping "unmodeled" contaminant
(whose profile is withheld), SEUDO should recover the found cell's true
activity more accurately than naive least squares during frames where the
contaminant is active.
"""

import numpy as np
import pytest

from seudo._native import NATIVE_AVAILABLE
from seudo.estimate import estimate_time_courses_with_seudo

requires_native = pytest.mark.skipif(not NATIVE_AVAILABLE, reason='native extension not built')


def gaussian_patch(shape, center, sigma, amplitude=1.0):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    return amplitude * np.exp(-(((yy - center[0]) ** 2 + (xx - center[1]) ** 2) / (2 * sigma ** 2)))


def make_synthetic_movie(seed=0):
    rng = np.random.default_rng(seed)
    y, x, f = 30, 30, 40

    cell_profile = gaussian_patch((y, x), center=(15, 15), sigma=2.0)
    cell_profile[cell_profile < 0.05 * cell_profile.max()] = 0.0

    contaminant_profile = gaussian_patch((y, x), center=(19, 19), sigma=2.0)
    contaminant_profile[contaminant_profile < 0.05 * contaminant_profile.max()] = 0.0

    cell_activity = np.zeros(f)
    contaminant_activity = np.zeros(f)

    # cell active alone in frames 5-9, contaminant active alone in frames 15-19,
    # both active together in frames 25-29
    cell_activity[5:10] = 4.0
    contaminant_activity[15:20] = 12.0
    cell_activity[25:30] = 4.0
    contaminant_activity[25:30] = 12.0

    movie = np.zeros((y, x, f))
    for ff in range(f):
        frame = (cell_profile * cell_activity[ff]
                 + contaminant_profile * contaminant_activity[ff]
                 + 0.02 * rng.normal(size=(y, x)))
        movie[:, :, ff] = frame

    profiles = cell_profile[:, :, np.newaxis]  # only the "found" cell is known

    return movie, profiles, cell_activity, contaminant_activity


def test_seudo_rejects_contamination_better_than_lsq():
    movie, profiles, cell_activity, contaminant_activity = make_synthetic_movie()

    result = estimate_time_courses_with_seudo(
        movie, profiles,
        p=1e-4, sigma2=0.02, lambda_blob=1.0, blob_radius=1.2,
        pad_space=6, verbose=False, solver_max_iter=300,
    )

    tc_seudo = result['tc'][:, 0]
    tc_lsq = result['tc_lsq'][:, 0]

    assert tc_seudo.shape == (40,)
    assert not np.any(np.isnan(tc_seudo))

    contaminant_only = slice(15, 20)  # contaminant active, cell truly silent

    lsq_error = np.abs(tc_lsq[contaminant_only] - cell_activity[contaminant_only])
    seudo_error = np.abs(tc_seudo[contaminant_only] - cell_activity[contaminant_only])

    # naive LSQ should be fooled into reporting spurious activity here
    assert np.mean(lsq_error) > 1.0
    # SEUDO should substantially reduce that false signal
    assert np.mean(seudo_error) < np.mean(lsq_error)

    cell_only = slice(5, 10)  # cell truly active, contaminant silent
    assert np.mean(np.abs(tc_seudo[cell_only] - cell_activity[cell_only])) < 1.0


COMMON_KWARGS = dict(p=1e-4, sigma2=0.02, lambda_blob=1.0, blob_radius=1.2,
                      pad_space=6, verbose=False, solver_max_iter=300)


def test_frame_blocks_restricts_to_specified_frames_and_matches_full_run():
    # port of @seudo/parallelSEUDO.m's keepFrames restriction
    movie, profiles, _, _ = make_synthetic_movie()

    full = estimate_time_courses_with_seudo(movie, profiles, **COMMON_KWARGS)
    blocks = [(5, 9), (25, 29)]
    restricted = estimate_time_courses_with_seudo(movie, profiles, frame_blocks=blocks, **COMMON_KWARGS)

    tc_full = full['tc'][:, 0]
    tc_restricted = restricted['tc'][:, 0]

    covered = np.zeros(40, dtype=bool)
    for s, e in blocks:
        covered[s:e + 1] = True

    assert np.all(np.isnan(tc_restricted[~covered]))
    assert not np.any(np.isnan(tc_restricted[covered]))
    # ds_time=1 (the default) means each frame is computed independently,
    # so restricting to a subset of frames must not change their values
    np.testing.assert_allclose(tc_restricted[covered], tc_full[covered])


def test_frame_blocks_resets_moving_average_between_blocks():
    movie, profiles, _, _ = make_synthetic_movie()
    kwargs = dict(COMMON_KWARGS, ds_time=3)

    multi_block = estimate_time_courses_with_seudo(movie, profiles, frame_blocks=[(0, 4), (25, 29)], **kwargs)
    single_block = estimate_time_courses_with_seudo(movie, profiles, frame_blocks=[(25, 29)], **kwargs)

    # the second block's results must be identical whether or not an earlier,
    # unrelated block was processed first -- the ds_time moving average must
    # reset at each block's start, not carry over stale frames across the gap
    np.testing.assert_allclose(multi_block['tc'][25:30, 0], single_block['tc'][25:30, 0])


def test_overlapping_frame_blocks_are_merged():
    movie, profiles, _, _ = make_synthetic_movie()
    merged = estimate_time_courses_with_seudo(movie, profiles, frame_blocks=[(5, 20), (10, 29)], **COMMON_KWARGS)
    assert merged['params']['frame_blocks'] == [(5, 29)]


class _FakeLazyMovie:
    """Minimal stand-in for matlab_io.Hdf5Movie -- only .shape and .get_frame."""

    def __init__(self, arr):
        self._arr = arr
        self.shape = arr.shape

    def get_frame(self, frame_index):
        return self._arr[:, :, frame_index].astype(float)


def test_lazy_movie_source_matches_in_memory_array():
    movie, profiles, _, _ = make_synthetic_movie()

    from_array = estimate_time_courses_with_seudo(movie, profiles, **COMMON_KWARGS)
    from_lazy = estimate_time_courses_with_seudo(_FakeLazyMovie(movie), profiles, **COMMON_KWARGS)

    np.testing.assert_allclose(from_array['tc'], from_lazy['tc'])


def test_blob_spacing_defaults_to_one_and_is_recorded_in_params():
    movie, profiles, _, _ = make_synthetic_movie()
    result = estimate_time_courses_with_seudo(movie, profiles, **COMMON_KWARGS)
    assert result['params']['blob_spacing'] == 1.0


@requires_native
def test_blob_spacing_changes_native_result():
    # blob_spacing only affects the native solver's internal grid -- default
    # (1, one blob per pixel) vs coarser (4) should give genuinely different
    # results, confirming the parameter is actually threaded through, not
    # silently ignored
    movie, profiles, _, _ = make_synthetic_movie()
    fine = estimate_time_courses_with_seudo(movie, profiles, use_native=True, blob_spacing=1, **COMMON_KWARGS)
    coarse = estimate_time_courses_with_seudo(movie, profiles, use_native=True, blob_spacing=4, **COMMON_KWARGS)
    assert not np.allclose(fine['tc'], coarse['tc'])


def test_blob_spacing_is_a_noop_for_python_fallback():
    # the pure-Python solver has no coarser-grid concept -- always one blob
    # per pixel regardless of blob_spacing
    movie, profiles, _, _ = make_synthetic_movie()
    fine = estimate_time_courses_with_seudo(movie, profiles, use_native=False, blob_spacing=1, **COMMON_KWARGS)
    coarse = estimate_time_courses_with_seudo(movie, profiles, use_native=False, blob_spacing=4, **COMMON_KWARGS)
    np.testing.assert_array_equal(fine['tc'], coarse['tc'])


@requires_native
def test_blob_spacing_still_rejects_contamination_at_coarser_grid():
    # the realSEUDO paper's core claim: a coarser blob grid shouldn't
    # substantially degrade false-transient rejection quality
    movie, profiles, cell_activity, _ = make_synthetic_movie()
    result = estimate_time_courses_with_seudo(
        movie, profiles, use_native=True, blob_spacing=4, **COMMON_KWARGS)
    tc_seudo = result['tc'][:, 0]

    contaminant_only = slice(15, 20)
    seudo_error = np.abs(tc_seudo[contaminant_only] - cell_activity[contaminant_only])
    assert np.mean(seudo_error) < 1.0  # still suppresses the false signal reasonably well

    cell_only = slice(5, 10)
    assert np.mean(np.abs(tc_seudo[cell_only] - cell_activity[cell_only])) < 1.0
