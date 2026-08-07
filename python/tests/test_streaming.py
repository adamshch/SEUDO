"""Tests for the streaming/online SEUDO pipeline (seudo/streaming.py):
realSEUDOfit() fits one frame at a time against a growing cell-profile set,
detecting and promoting new cells as they appear.
"""

import numpy as np
import pytest

from seudo._native import NATIVE_AVAILABLE
from seudo.estimate import estimate_time_courses_with_seudo
from seudo.geometry import compute_roi_coms
from seudo.streaming import DetectionParams, FitParams, StreamingState, TilingConfig, realSEUDOfit

requires_native = pytest.mark.skipif(not NATIVE_AVAILABLE, reason='native extension not built')

COMMON_FIT = dict(p=1e-4, sigma2=0.02, lambda_blob=1.0, blob_radius=1.2, pad_space=6,
                   solver_tol=0.01, solver_max_iter=500, use_native=False)


def gaussian_patch(shape, center, sigma, amplitude=1.0):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    return amplitude * np.exp(-(((yy - center[0]) ** 2 + (xx - center[1]) ** 2) / (2 * sigma ** 2)))


def test_parity_with_offline_solver_frame_by_frame():
    # one known cell, one spatially-overlapping unmodeled contaminant --
    # same scenario shape as test_estimate_synthetic.py's, exercising the
    # same LSQ-vs-sparse decision branching, not just a trivial single-profile fit.
    rng = np.random.default_rng(0)
    y, x, f = 30, 30, 40
    cell_profile = gaussian_patch((y, x), center=(15, 15), sigma=2.0)
    cell_profile[cell_profile < 0.05 * cell_profile.max()] = 0.0
    contam_profile = gaussian_patch((y, x), center=(19, 19), sigma=2.0)
    contam_profile[contam_profile < 0.05 * contam_profile.max()] = 0.0

    cell_activity = np.zeros(f)
    contam_activity = np.zeros(f)
    cell_activity[5:10] = 4.0
    contam_activity[15:20] = 12.0
    cell_activity[25:30] = 4.0
    contam_activity[25:30] = 12.0

    movie = np.zeros((y, x, f))
    for ff in range(f):
        movie[:, :, ff] = (cell_profile * cell_activity[ff] + contam_profile * contam_activity[ff]
                            + 0.02 * rng.normal(size=(y, x)))

    profiles = cell_profile[:, :, np.newaxis]

    offline = estimate_time_courses_with_seudo(movie, profiles, ds_time=1, **COMMON_FIT)
    offline_tc = offline['tc'][:, 0]

    state = StreamingState((y, x), initial_profiles=profiles, fit=FitParams(**COMMON_FIT))
    streaming_tc = np.array([realSEUDOfit(movie[:, :, ff], state).activity[0] for ff in range(f)])

    np.testing.assert_allclose(streaming_tc, offline_tc, atol=1e-8)


def make_two_cell_movie(seed=1, onset=50, y=30, x=60, f=100):
    yy, xx = np.mgrid[0:y, 0:x]
    prof0 = gaussian_patch((y, x), center=(15, 15), sigma=2.0)
    prof0[prof0 < 0.05 * prof0.max()] = 0.0
    prof1 = gaussian_patch((y, x), center=(15, 45), sigma=2.0)
    prof1[prof1 < 0.05 * prof1.max()] = 0.0

    rng = np.random.default_rng(seed)
    act0 = np.full(f, 2.0)
    act1 = np.zeros(f)
    act1[onset:] = 3.0

    movie = np.zeros((y, x, f))
    for ff in range(f):
        movie[:, :, ff] = prof0 * act0[ff] + prof1 * act1[ff] + 0.01 * rng.normal(size=(y, x))

    return movie, prof0, prof1, act0, act1


def default_streaming_fit():
    return FitParams(sigma2=0.002, lambda_blob=10.0, blob_radius=3.0, pad_space=5, use_native=False)


def test_bounded_latency_new_cell_detection():
    onset = 50
    movie, prof0, prof1, act0, act1 = make_two_cell_movie(onset=onset)
    y, x, f = movie.shape

    state = StreamingState((y, x), initial_profiles=prof0[:, :, np.newaxis], fit=default_streaming_fit())

    promoted_at = None
    for ff in range(f):
        result = realSEUDOfit(movie[:, :, ff], state)
        if result.new_cells:
            assert promoted_at is None, 'should only promote once'
            promoted_at = ff
            assert result.new_cells == [1]
            assert 1 in result.activity

    assert promoted_at is not None, 'the new cell was never detected'
    assert onset <= promoted_at <= onset + 10, (
        f'promotion latency too high: onset={onset}, promoted at {promoted_at}')

    # activity tracks the true (steady) amplitude after promotion
    late_activity = [realSEUDOfit(movie[:, :, ff], state).activity[1] for ff in range(80, 90)]
    assert np.allclose(late_activity, 3.0, atol=0.5)


def test_contamination_never_promoted():
    y, x, f = 30, 60, 100
    prof0 = gaussian_patch((y, x), center=(15, 15), sigma=2.0)
    prof0[prof0 < 0.05 * prof0.max()] = 0.0
    prof1 = gaussian_patch((y, x), center=(15, 45), sigma=2.0)
    prof1[prof1 < 0.05 * prof1.max()] = 0.0

    rng = np.random.default_rng(2)
    act0 = np.full(f, 2.0)
    act1 = np.zeros(f)
    # flickers on for 2 frames, off for 3, repeating -- never sustains
    # PromotionParams.consecutive_frames_required (default 5) consecutive matches
    for i in range(0, f, 5):
        act1[i:i + 2] = 3.0

    movie = np.zeros((y, x, f))
    for ff in range(f):
        movie[:, :, ff] = prof0 * act0[ff] + prof1 * act1[ff] + 0.01 * rng.normal(size=(y, x))

    state = StreamingState((y, x), initial_profiles=prof0[:, :, np.newaxis], fit=default_streaming_fit())
    for ff in range(f):
        result = realSEUDOfit(movie[:, :, ff], state)
        assert result.new_cells == []

    assert state.profiles.shape[2] == 1


def test_promoted_profile_recovers_true_shape():
    onset = 50
    movie, prof0, prof1, act0, act1 = make_two_cell_movie(onset=onset)
    y, x, f = movie.shape

    state = StreamingState((y, x), initial_profiles=prof0[:, :, np.newaxis], fit=default_streaming_fit())
    for ff in range(f):
        realSEUDOfit(movie[:, :, ff], state)

    assert state.profiles.shape[2] == 2
    true_coms, _ = compute_roi_coms(prof1)
    built_coms, _ = compute_roi_coms(state.profiles[:, :, 1])
    np.testing.assert_allclose(built_coms[0], true_coms[0], atol=1.0)


def test_tiling_boundary_no_double_detect_or_miss():
    y, x, f = 40, 40, 60
    prof = gaussian_patch((y, x), center=(20, 20), sigma=2.0)  # exactly on a 20x20 tile-grid corner
    prof[prof < 0.05 * prof.max()] = 0.0

    rng = np.random.default_rng(3)
    act = np.zeros(f)
    act[10:] = 3.0
    movie = np.zeros((y, x, f))
    for ff in range(f):
        movie[:, :, ff] = prof * act[ff] + 0.01 * rng.normal(size=(y, x))

    tiling = TilingConfig(tile_shape=(20, 20), overlap=15)
    state = StreamingState((y, x), tiling=tiling, fit=default_streaming_fit())
    for ff in range(f):
        realSEUDOfit(movie[:, :, ff], state)

    assert state.profiles.shape[2] == 1, 'should promote exactly once, not zero or two'


def test_tiling_vs_no_tiling_agreement():
    onset = 50
    movie, prof0, prof1, act0, act1 = make_two_cell_movie(onset=onset)
    y, x, f = movie.shape

    fit = default_streaming_fit()

    state_untiled = StreamingState((y, x), initial_profiles=prof0[:, :, np.newaxis], fit=fit)
    untiled_activity0 = []
    promoted_untiled = None
    for ff in range(f):
        r = realSEUDOfit(movie[:, :, ff], state_untiled)
        untiled_activity0.append(r.activity[0])
        if r.new_cells:
            promoted_untiled = ff

    state_tiled = StreamingState(
        (y, x), initial_profiles=prof0[:, :, np.newaxis], fit=fit,
        tiling=TilingConfig(tile_shape=(15, 20), overlap=15))
    tiled_activity0 = []
    promoted_tiled = None
    for ff in range(f):
        r = realSEUDOfit(movie[:, :, ff], state_tiled)
        tiled_activity0.append(r.activity[0])
        if r.new_cells:
            promoted_tiled = ff

    assert promoted_untiled == promoted_tiled
    assert state_untiled.profiles.shape[2] == state_tiled.profiles.shape[2] == 2
    np.testing.assert_allclose(untiled_activity0, tiled_activity0, atol=1e-8)


def test_parallel_tile_detection_matches_sequential():
    # a cell centered exactly on a tile-grid corner exercises multiple
    # tiles' detection genuinely (see test_tiling_boundary_no_double_detect_or_miss)
    y, x, f = 40, 40, 60
    prof = gaussian_patch((y, x), center=(20, 20), sigma=2.0)
    prof[prof < 0.05 * prof.max()] = 0.0

    rng = np.random.default_rng(3)
    act = np.zeros(f)
    act[10:] = 3.0
    movie = np.zeros((y, x, f))
    for ff in range(f):
        movie[:, :, ff] = prof * act[ff] + 0.01 * rng.normal(size=(y, x))

    tiling = TilingConfig(tile_shape=(20, 20), overlap=15)

    def run(n_jobs):
        fit = FitParams(**{**vars(default_streaming_fit()), 'n_jobs': n_jobs})
        state = StreamingState((y, x), tiling=tiling, fit=fit)
        promoted_at = None
        for ff in range(f):
            r = realSEUDOfit(movie[:, :, ff], state)
            if r.new_cells:
                promoted_at = ff
        state.close()
        return promoted_at, state.profiles.shape[2]

    seq = run(1)
    par = run(4)
    assert seq == par


def test_blanking_skips_expensive_detection_for_empty_tiles(monkeypatch):
    # a tile with no signal at all should never reach ndimage.label --
    # verifies the "blanking" pre-filter is actually exercised, not just
    # present in the source with no effect
    import seudo.streaming as streaming_mod

    y, x, f = 40, 40, 15
    movie = np.zeros((y, x, f))  # pure background, no cell anywhere

    calls = []
    orig_label = streaming_mod.ndimage.label

    def spy_label(*args, **kwargs):
        calls.append(1)
        return orig_label(*args, **kwargs)

    monkeypatch.setattr(streaming_mod.ndimage, 'label', spy_label)

    state = StreamingState((y, x), fit=default_streaming_fit())
    for ff in range(f):
        realSEUDOfit(movie[:, :, ff], state)

    assert calls == [], 'ndimage.label should never be reached for an all-background frame'


def test_single_active_tile_never_reaches_the_executor(monkeypatch):
    # the sparse-activity regression this split was built to fix: with only
    # one tile ever non-blank, detection should never touch the thread pool
    # at all -- not "submit it and let it be fast," genuinely never call
    # submit() for detection, since submission/sync overhead on an
    # already-cheap tile is a net loss (see commit history/memory notes)
    y, x, f = 120, 120, 10
    # cell centered well inside one tile's core (>overlap px from every edge
    # of its own tile), so its signal can never leak into a neighboring
    # tile's halo -- see build_tiles: tile_shape=60/overlap=15 gives tile 0
    # core (0,59,0,59), and (30,30) is 30px from the nearest edge
    prof = gaussian_patch((y, x), center=(30, 30), sigma=2.0)
    prof[prof < 0.05 * prof.max()] = 0.0

    rng = np.random.default_rng(6)
    movie = np.zeros((y, x, f))
    for ff in range(f):
        movie[:, :, ff] = prof * 3.0 + 0.01 * rng.normal(size=(y, x))

    fit = FitParams(**{**vars(default_streaming_fit()), 'n_jobs': 4})
    tiling = TilingConfig(tile_shape=(60, 60), overlap=15)  # 2x2 = 4 tiles, only tile 0 ever near the cell
    state = StreamingState((y, x), tiling=tiling, fit=fit)  # no known cells -- isolates detection submissions

    submit_calls = []
    orig_submit = state._executor.submit

    def spy_submit(*args, **kwargs):
        submit_calls.append(args[0])  # the submitted function
        return orig_submit(*args, **kwargs)

    monkeypatch.setattr(state._executor, 'submit', spy_submit)

    for ff in range(f):
        realSEUDOfit(movie[:, :, ff], state)
    state.close()

    detect_submissions = [c for c in submit_calls if c.__name__ == '_detect_in_tile']
    assert detect_submissions == [], 'a single non-blank tile should never be submitted to the executor'


def make_multi_cell_movie(seed=4, n_cells=6, y=40, x=80, f=30):
    rng = np.random.default_rng(seed)
    xs = np.linspace(10, x - 10, n_cells)
    profiles = np.zeros((y, x, n_cells))
    for i, cx in enumerate(xs):
        p = gaussian_patch((y, x), center=(y // 2, cx), sigma=2.0)
        p[p < 0.05 * p.max()] = 0.0
        profiles[:, :, i] = p

    activity = rng.uniform(0.5, 4.0, size=(f, n_cells))
    movie = np.zeros((y, x, f))
    for ff in range(f):
        frame = np.zeros((y, x))
        for i in range(n_cells):
            frame += profiles[:, :, i] * activity[ff, i]
        movie[:, :, ff] = frame + 0.01 * rng.normal(size=(y, x))

    return movie, profiles


def test_n_jobs_gives_bit_identical_results_to_sequential():
    movie, profiles = make_multi_cell_movie()
    y, x, f = movie.shape

    fit_seq = default_streaming_fit()
    fit_par = FitParams(**{**vars(fit_seq), 'n_jobs': 4})

    state_seq = StreamingState((y, x), initial_profiles=profiles, fit=fit_seq)
    state_par = StreamingState((y, x), initial_profiles=profiles, fit=fit_par)
    try:
        for ff in range(f):
            r_seq = realSEUDOfit(movie[:, :, ff], state_seq)
            r_par = realSEUDOfit(movie[:, :, ff], state_par)
            assert r_seq.activity.keys() == r_par.activity.keys()
            for cell_id in r_seq.activity:
                assert r_seq.activity[cell_id] == r_par.activity[cell_id]
    finally:
        state_par.close()


def test_executor_lifecycle():
    y, x = 20, 20
    state_seq = StreamingState((y, x), fit=default_streaming_fit())
    assert state_seq._executor is None  # n_jobs=1 -- no pool created

    fit_par = FitParams(**{**vars(default_streaming_fit()), 'n_jobs': 4})
    with StreamingState((y, x), fit=fit_par) as state_par:
        assert state_par._executor is not None
    assert state_par._executor is None  # closed on context-manager exit


def test_state_contract_sanity():
    onset = 20
    movie, prof0, prof1, act0, act1 = make_two_cell_movie(onset=onset, f=40)
    y, x, f = movie.shape

    state = StreamingState((y, x), initial_profiles=prof0[:, :, np.newaxis], fit=default_streaming_fit())
    state_id = id(state)

    n_before = state.profiles.shape[2]
    for ff in range(f):
        result = realSEUDOfit(movie[:, :, ff], state)
        assert id(state) == state_id, 'state should mutate in place, not be replaced'

        if result.new_cells:
            assert state.profiles.shape[2] == n_before + len(result.new_cells)
            for cell_id in result.new_cells:
                assert cell_id in result.activity
            n_before = state.profiles.shape[2]
        else:
            assert state.profiles.shape[2] == n_before


def make_contaminated_frames(seed=5, y=30, x=30, n_frames=20):
    """A single known cell plus a withheld, overlapping 'contaminant' --
    unlike make_multi_cell_movie (clean, no interference, so LSQ always
    wins the cost comparison and the blob dictionary is never engaged
    regardless of spacing), this scenario actually exercises the blob
    branch, which is where blob_spacing has any effect at all."""
    yy, xx = np.mgrid[0:y, 0:x]
    prof0 = gaussian_patch((y, x), center=(15, 15), sigma=2.0)
    prof0[prof0 < 0.05 * prof0.max()] = 0.0
    contaminant = gaussian_patch((y, x), center=(19, 19), sigma=2.0)
    contaminant[contaminant < 0.05 * contaminant.max()] = 0.0

    rng = np.random.default_rng(seed)
    frames = [prof0 * 4.0 + contaminant * 12.0 + 0.02 * rng.normal(size=(y, x)) for _ in range(n_frames)]
    return frames, prof0[:, :, np.newaxis], (y, x)


@requires_native
def test_blob_spacing_changes_native_result():
    frames, profiles, (y, x) = make_contaminated_frames()

    def run(blob_spacing):
        fit = FitParams(sigma2=0.02, lambda_blob=1.0, blob_radius=1.2, pad_space=6,
                         use_native=True, blob_spacing=blob_spacing, solver_max_iter=300)
        state = StreamingState((y, x), initial_profiles=profiles, fit=fit)
        return [realSEUDOfit(frame, state).activity[0] for frame in frames]

    fine = run(1)
    coarse = run(4)
    assert not np.allclose(fine, coarse)


def test_blob_spacing_is_a_noop_for_python_fallback():
    frames, profiles, (y, x) = make_contaminated_frames()

    def run(blob_spacing):
        fit = FitParams(sigma2=0.02, lambda_blob=1.0, blob_radius=1.2, pad_space=6,
                         use_native=False, blob_spacing=blob_spacing, solver_max_iter=300)
        state = StreamingState((y, x), initial_profiles=profiles, fit=fit)
        return [realSEUDOfit(frame, state).activity[0] for frame in frames]

    fine = run(1)
    coarse = run(4)
    np.testing.assert_array_equal(fine, coarse)


def test_noise_grid_shape_default_matches_old_single_scalar_behavior():
    # noise_grid_shape=(1,1) (the default) must reproduce the pre-existing
    # single-global-scalar behavior exactly -- this is the backward
    # compatibility guarantee for every earlier test in this file
    onset = 20
    movie, prof0, prof1, act0, act1 = make_two_cell_movie(onset=onset, f=40)
    y, x, f = movie.shape

    fit = default_streaming_fit()
    default_detection = DetectionParams()
    assert default_detection.noise_grid_shape == (1, 1)

    state = StreamingState((y, x), initial_profiles=prof0[:, :, np.newaxis], fit=fit, detection=default_detection)
    activity = [realSEUDOfit(movie[:, :, ff], state).activity[0] for ff in range(f)]

    state_ref = StreamingState((y, x), initial_profiles=prof0[:, :, np.newaxis], fit=fit)
    activity_ref = [realSEUDOfit(movie[:, :, ff], state_ref).activity[0] for ff in range(f)]

    np.testing.assert_array_equal(activity, activity_ref)


def test_spatially_varying_noise_detects_modest_cell_faster_in_quiet_region():
    # a modest-amplitude cell in a quiet region, next to an unrelated noisy
    # region -- a single global noise scalar gets dragged up by the noisy
    # region, delaying detection of the modest cell; per-section estimation
    # should detect it much faster since the quiet region's own noise floor
    # is correctly low from the start
    y, x, f = 60, 60, 30
    yy, xx = np.mgrid[0:y, 0:x]
    prof = gaussian_patch((y, x), center=(30, 15), sigma=2.0)  # in the quiet left half
    prof[prof < 0.05 * prof.max()] = 0.0

    rng = np.random.default_rng(2)
    act = np.zeros(f)
    act[10:] = 1.2  # modest amplitude

    frames = []
    for ff in range(f):
        frame = prof * act[ff]
        noise = np.zeros((y, x))
        noise[:, :30] = 0.01 * rng.normal(size=(y, 30))       # quiet left half
        noise[:, 30:] = 0.15 * rng.normal(size=(y, x - 30))   # noisy right half
        frames.append(frame + noise)

    def run(noise_grid_shape):
        fit = default_streaming_fit()
        detection = DetectionParams(noise_grid_shape=noise_grid_shape)
        state = StreamingState((y, x), fit=fit, detection=detection)
        promoted_at = None
        for ff in range(f):
            r = realSEUDOfit(frames[ff], state)
            if r.new_cells:
                promoted_at = ff
        return promoted_at, state.profiles.shape[2]

    global_promoted, global_n = run((1, 1))
    split_promoted, split_n = run((1, 2))  # left/right halves get separate noise estimates

    assert global_n == split_n == 1, 'both should eventually find the cell'
    assert split_promoted < global_promoted, (
        'spatially-varying noise estimation should detect the modest cell faster '
        'than a global scalar dragged up by the unrelated noisy region')


def test_subtract_tracked_contributions_removes_tracked_signal():
    # a single active (not-yet-promoted) candidate track whose accumulated
    # history already captures a clean cell shape -- stage 1 of two-stage
    # detection should recognize and remove most of its contribution from a
    # fresh residual, not leave it untouched (a no-op would defeat the whole
    # point: telling an already-tracked candidate's own signal apart from
    # something genuinely new -- see _subtract_tracked_contributions).
    from seudo.streaming import CandidateTrack, _subtract_tracked_contributions

    y, x = 30, 30
    prof = gaussian_patch((y, x), center=(15, 15), sigma=2.0)
    prof[prof < 0.05 * prof.max()] = 0.0

    state = StreamingState((y, x), fit=default_streaming_fit())
    bbox = (10, 20, 10, 20)
    crop = prof[bbox[0]:bbox[1] + 1, bbox[2]:bbox[3] + 1]
    mask = crop > 0
    history_crop = crop * 3.0
    track = CandidateTrack(
        centroid=(15.0, 15.0), consecutive_frames=3, gap=0,
        history=[(0, bbox, mask, history_crop)] * 3,
    )
    state.candidate_tracks[0] = track

    rng = np.random.default_rng(0)
    residual = prof * 3.0 + 0.01 * rng.normal(size=(y, x))
    residual2, exclude_mask = _subtract_tracked_contributions(state, residual)

    footprint = prof > 0
    assert np.abs(residual2[footprint]).sum() < 0.5 * np.abs(residual[footprint]).sum(), (
        'stage 1 should explain away most of an already-tracked candidate\'s own signal')
    assert exclude_mask[15, 15], 'the tracked footprint should be excluded from new-candidate scanning'


def test_tracked_candidate_never_spawns_a_duplicate_track():
    # while a real cell is being tracked (pre-promotion), the two-stage
    # subtraction should keep explaining its own signal away every frame --
    # if it didn't, leftover residual at the same location could spuriously
    # trigger a second, duplicate track for the same physical source
    onset = 30
    movie, prof0, prof1, act0, act1 = make_two_cell_movie(onset=onset, f=onset + 15)
    y, x, f = movie.shape

    state = StreamingState((y, x), initial_profiles=prof0[:, :, np.newaxis], fit=default_streaming_fit())
    for ff in range(f):
        realSEUDOfit(movie[:, :, ff], state)

    assert state._next_track_id == 1, (
        f'expected exactly one track ever created for the single new cell, got {state._next_track_id}')


def test_ds_time_default_is_a_noop():
    # ds_time=1 (the default) must apply no averaging at all -- sigma2_ds
    # equals the raw sigma2 exactly, and a single call's avg_frame equals
    # the input frame exactly (np.mean of one array introduces no
    # floating-point change), so behavior is bit-identical to every test
    # above that never mentions ds_time at all.
    y, x = 20, 20
    fit = FitParams(sigma2=0.02)
    state = StreamingState((y, x), fit=fit)
    assert state.sigma2_ds == fit.sigma2

    frame = np.random.default_rng(0).normal(size=(y, x))
    state._frame_buffer.append(frame)
    avg_frame = np.mean(state._frame_buffer, axis=0)
    np.testing.assert_array_equal(avg_frame, frame)


def test_ds_time_smooths_a_stepped_activity_trace():
    # a single known cell with a noise-free step in true activity -- with
    # ds_time=3, the FIRST frame after the step should read back as a
    # partial (averaged-with-the-past) value, not the new true value
    # immediately, and should reach the new true value only once the
    # trailing window has fully slid past the step -- direct evidence the
    # averaging math is doing the right causal thing, not just "some
    # smoothing happens somewhere."
    y, x, f = 30, 30, 20
    prof = gaussian_patch((y, x), center=(15, 15), sigma=2.0)
    prof[prof < 0.05 * prof.max()] = 0.0

    act = np.full(f, 2.0)
    step_at = 10
    act[step_at:] = 8.0

    movie = np.stack([prof * act[ff] for ff in range(f)], axis=2)

    fit = FitParams(**{**vars(default_streaming_fit()), 'ds_time': 3})
    state = StreamingState((y, x), initial_profiles=prof[:, :, np.newaxis], fit=fit)
    activity = [realSEUDOfit(movie[:, :, ff], state).activity[0] for ff in range(f)]

    # before the step: steady at the true value
    assert np.allclose(activity[step_at - 1], 2.0, atol=1e-6)
    # right at the step: averages 1 new frame at 8.0 with 2 trailing frames at 2.0
    assert np.allclose(activity[step_at], (8.0 + 2.0 + 2.0) / 3.0, atol=1e-6)
    # one frame later: 2 frames at 8.0, 1 trailing frame at 2.0
    assert np.allclose(activity[step_at + 1], (8.0 + 8.0 + 2.0) / 3.0, atol=1e-6)
    # once the window has fully slid past the step: back to steady, at the new value
    assert np.allclose(activity[step_at + 2], 8.0, atol=1e-6)
    assert np.allclose(activity[-1], 8.0, atol=1e-6)


def test_ds_time_improves_detection_in_high_noise():
    # a modest-amplitude cell under noise high enough to make single-frame
    # detection unreliable -- ds_time=3's trailing average should reduce
    # the effective per-pixel noise and detect it, while ds_time=1 (no
    # averaging) fails to promote at all within the same frame budget
    y, x, f = 30, 30, 100
    prof = gaussian_patch((y, x), center=(15, 15), sigma=2.0)
    prof[prof < 0.05 * prof.max()] = 0.0

    rng = np.random.default_rng(7)
    act = np.zeros(f)
    act[10:] = 1.0

    frames = [prof * act[ff] + 0.3 * rng.normal(size=(y, x)) for ff in range(f)]

    def run(ds_time):
        fit = FitParams(**{**vars(default_streaming_fit()), 'ds_time': ds_time})
        state = StreamingState((y, x), fit=fit)
        promoted_at = None
        for ff in range(f):
            r = realSEUDOfit(frames[ff], state)
            if r.new_cells:
                promoted_at = ff
        return promoted_at, state.profiles.shape[2]

    no_avg_promoted, no_avg_n = run(1)
    smoothed_promoted, smoothed_n = run(3)

    assert smoothed_n == 1, 'ds_time=3 should reliably detect the modest cell under high noise'
    assert no_avg_n == 0 or (no_avg_promoted is not None and smoothed_promoted < no_avg_promoted), (
        'ds_time=3 should detect the cell either faster than ds_time=1, or where ds_time=1 fails entirely')
