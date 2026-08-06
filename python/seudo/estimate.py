"""Port of @seudo/estimateTimeCoursesWithSEUDO.m -- 'static' dyn_mod path only.

This is a proof-of-concept port of the core SEUDO solver. Each frame is fit
independently as

    minimize_{x >= 0}  0.5 * ||A x - frame||^2 + lambda' * x

where A's columns are the (padded, windowed) known cell profiles plus a
bank of Gaussian "blob" basis functions covering every pixel in the window
-- one blob per pixel, meant to absorb activity from unmodeled ("not
found") sources. The frame's cell activation is then taken from whichever
of (a) plain least squares against the known profiles, or (b) this sparse
fit, has lower cost -- matching the MATLAB logic exactly (see
estimateTimeCoursesWithSEUDO.m lines ~582-600).

NOT yet ported from the MATLAB version:
    - 'rwl1df' / 'bpdndf' dynamic modes (temporal regularization across frames)
    - tif / function-handle movie sources (in-memory arrays and lazy
      get_frame()-style sources -- see matlab_io.Hdf5Movie -- are both fine)
    - saveBlobTimeCourse / saveOtherCellTimeCourses output options

frame_blocks (port of @seudo/parallelSEUDO.m's keepFrames restriction):
    an optional list of (start, end) 0-indexed INCLUSIVE frame ranges. When
    given, only frames within these ranges are processed -- everywhere else
    in the output arrays stays NaN. This is what makes running SEUDO
    tractable on a real, tens-of-thousands-of-frame movie: restrict it to
    just the frames spanned by a cell's detected transients (its
    transient_info times) rather than the whole movie. The moving-average
    (ds_time) window is reset at the start of each block, rather than
    blending across the (possibly large) gap between two unrelated blocks.
"""

from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy.signal import convolve2d

from .blob import make_seudo_blob
from .geometry import compute_roi_coms
from .solver import fista_nonneg_weighted_l1
from ._native import seudo_solve as _native_seudo_solve, NATIVE_AVAILABLE as _NATIVE_AVAILABLE


def _get_frame(movie, frame_index, zero_level):
    if hasattr(movie, 'get_frame'):
        return movie.get_frame(frame_index) - zero_level
    return movie[:, :, frame_index].astype(float) - zero_level


def _setup_cell_window(
    profiles, this_cell, y0, y1, x0, x1, min_pix_for_inclusion,
    lambda_prof, lambda_blob, sigma2_ds, p, one_blob,
):
    """Build the ROI/blob-basis/lambda setup for one cell's window -- shared
    by the per-frame loop below and seudo_residual.py's per-transient
    single-image fit, so both use identical math for "what SEUDO regresses
    against" in a given window. y0/y1/x0/x1 are the window bounds (already
    padded/clamped by the caller)."""
    n_y = y1 - y0 + 1
    n_x = x1 - x0 + 1
    n_blobs = n_y * n_x

    pix_per_profile = np.sum(profiles[y0:y1 + 1, x0:x1 + 1, :] > 0, axis=(0, 1))
    include = pix_per_profile > min_pix_for_inclusion
    include[this_cell] = True

    n_cells_window = int(np.sum(include))
    window_profiles = profiles[y0:y1 + 1, x0:x1 + 1, include]
    rois = window_profiles.reshape(n_y * n_x, n_cells_window)

    prof_norms = np.sqrt(np.sum(rois ** 2, axis=0))

    lambdas0 = np.concatenate([
        np.full(n_cells_window, lambda_prof),
        np.full(n_blobs, lambda_blob),
    ])
    k1 = 2 * sigma2_ds * lambdas0
    pos = lambdas0 > 0
    k2 = 2 * sigma2_ds * (np.log(1.0 / p - 1.0) - np.sum(np.log(lambdas0[pos])))

    norm_factors = np.concatenate([prof_norms, np.ones(n_blobs)])
    lambdas = k1 / norm_factors

    cell_index_within = int(np.where(np.flatnonzero(include) == this_cell)[0][0])

    rois_scaled = rois / prof_norms[np.newaxis, :]

    def A(z):
        z_cells = z[:n_cells_window]
        z_blob = z[n_cells_window:].reshape(n_y, n_x)
        out = rois_scaled @ z_cells
        out = out + convolve2d(z_blob, one_blob, mode='same').ravel()
        return out

    def At(v):
        top = rois_scaled.T @ v
        bottom = convolve2d(v.reshape(n_y, n_x), one_blob, mode='same').ravel()
        return np.concatenate([top, bottom])

    return dict(
        n_y=n_y, n_x=n_x, rois=rois, rois_scaled=rois_scaled, norm_factors=norm_factors,
        lambdas=lambdas, k1=k1, k2=k2, cell_index_within=cell_index_within,
        n_cells_window=n_cells_window, operators=(A, At), include=include,
    )


def _solve_one_frame_cell(
    this_frame, rois, rois_scaled, lambdas, norm_factors, k1, k2,
    n_y, n_x, one_blob, operators_cc,
    use_native, native_l_mode, native_nthreads, solver_tol, solver_max_iter,
):
    """Solve one (frame, cell) pair -- independent of every other frame and
    cell, which is what makes this safe to run concurrently (see n_jobs on
    estimate_time_courses_with_seudo). No shared mutable state is touched:
    everything here is either a local or a read-only array/closure."""
    n_cells_window = rois.shape[1]
    tc_lsq_frame = np.linalg.solve(rois.T @ rois, rois.T @ this_frame)

    if use_native:
        fit_weights, _n_steps, _log = _native_seudo_solve(
            this_frame.reshape(n_y, n_x), one_blob, 1.0,
            rois_scaled, np.zeros(0), lambdas,
            solver_tol, solver_max_iter, native_l_mode, 0, False, native_nthreads,
        )
    else:
        A_cc, At_cc = operators_cc
        x0 = np.zeros(n_cells_window + n_y * n_x)
        fit_weights = fista_nonneg_weighted_l1(
            A_cc, At_cc, this_frame, lambdas, x0,
            tol=solver_tol, max_iter=solver_max_iter,
        )
    fit_weights = fit_weights / norm_factors
    fit_x = fit_weights[:n_cells_window]
    fit_bob = fit_weights[n_cells_window:]

    lsq_cost = np.sum((rois @ tc_lsq_frame - this_frame) ** 2)
    blob_contrib = convolve2d(fit_bob.reshape(n_y, n_x), one_blob, mode='same').ravel()
    bob_cost = (np.sum((this_frame - rois @ fit_x - blob_contrib) ** 2)
                + np.sum(np.abs(k1 * fit_weights)) - k2)

    if lsq_cost < bob_cost:
        fit_fancy = tc_lsq_frame
        fit_x = np.zeros_like(fit_x)
    else:
        fit_fancy = fit_x

    return tc_lsq_frame, fit_fancy, fit_x, lsq_cost, bob_cost


def _merge_frame_blocks(blocks):
    blocks = sorted((int(s), int(e)) for s, e in blocks)
    merged = []
    for s, e in blocks:
        if merged and s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def estimate_time_courses_with_seudo(
    movie,
    profiles,
    which_cells=None,
    p=1e-5,
    sigma2=0.01,
    lambda_blob=20.0,
    blob_radius=1.2,
    ds_time=1,
    lambda_prof=0.0,
    min_pix_for_inclusion=1,
    pad_space=10,
    use_com=False,
    zero_level=0.0,
    solver_tol=0.01,
    solver_max_iter=1000,
    verbose=True,
    frame_blocks=None,
    use_native='auto',
    native_l_mode=2,
    native_nthreads=1,
    n_jobs=1,
):
    """movie: (Y, X, F) array, or an object with `.shape` and `.get_frame(i)`
    for lazy access (see matlab_io.Hdf5Movie) -- movie data is only ever
    touched frame-by-frame, so a lazy source combined with `frame_blocks`
    avoids ever materializing the full movie.
    profiles: (Y, X, nCellsTotal) array.
    which_cells: 0-indexed sequence of cell indices to analyze (default: all).
    frame_blocks: see module docstring.

    use_native: 'auto' (use the compiled C++ accelerator -- see seudo/_native --
        if it's been built, else fall back to the pure-Python solver silently),
        True (require it, raising if not built), or False (always use the
        pure-Python solver). Both solve the identical problem (same nonneg
        weighted-L1 formulation, same lambda/normalization scheme); the native
        one is just compiled and multithreaded.
    native_l_mode: only used when the native solver runs -- 0-3, see
        _native/seudo_native_binding.cpp; 2 (dynamic L + fast brake) is a
        reasonable default.
    native_nthreads: threads the native solver uses *within* a single frame's
        solve (independent of any outer parallelism).
    n_jobs: number of (frame, cell) solves to run concurrently in a thread
        pool. Each is fully independent (see _solve_one_frame_cell), so this
        changes wall-clock time only, never the result. Most effective with
        use_native=True: that solver releases the GIL for the whole compute()
        call, so real multi-core parallelism happens even though Python
        threads (not processes) are used; with the pure-Python solver the
        benefit is smaller (NumPy releases the GIL for individual BLAS/ufunc
        calls, but not for the whole per-frame loop). n_jobs=1 (default) runs
        sequentially, identical to earlier versions of this function.

    Returns a dict with keys 'tc', 'tc_lsq', 'params', 'extras', mirroring
    se.tcSeudo(idx) in the MATLAB version.
    """
    if use_native == 'auto':
        use_native = _NATIVE_AVAILABLE
    elif use_native and not _NATIVE_AVAILABLE:
        raise RuntimeError(
            'use_native=True but the compiled native SEUDO extension is not available; '
            'build it with seudo/_native/build_native.sh, or pass use_native=False'
        )
    profiles = np.asarray(profiles, dtype=float)
    mov_y, mov_x, n_frames = movie.shape
    n_cells_total = profiles.shape[2]

    if which_cells is None:
        which_cells = np.arange(n_cells_total)
    else:
        which_cells = np.asarray(which_cells)
    n_cells_analyzed = len(which_cells)

    sigma2_ds = sigma2 / ds_time
    min_pix_for_inclusion = max(min_pix_for_inclusion, 1)

    one_blob = make_seudo_blob(blob_radius)

    # ---- per-cell precomputation (mirrors the setup loop before the frame loop) ----

    which_pixels = np.zeros((mov_y * mov_x, n_cells_analyzed), dtype=bool)
    which_profiles = np.zeros((n_cells_total, n_cells_analyzed), dtype=bool)

    n_y_list, n_x_list = [], []
    rois_list, rois_scaled_list, norm_factors_list, lambdas_list = [], [], [], []
    k1_list, k2_list = [], []
    cell_index_within_list = []
    operators = []

    for cc, this_cell in enumerate(which_cells):
        prof = profiles[:, :, this_cell]
        coms, outer_bounds = compute_roi_coms(prof)
        if use_com:
            cy, cx = int(round(coms[0, 1])), int(round(coms[0, 0]))
            y0 = y1 = cy
            x0 = x1 = cx
        else:
            y0, y1, x0, x1 = outer_bounds[0]

        y0 = max(0, y0 - pad_space)
        y1 = min(mov_y - 1, y1 + pad_space)
        x0 = max(0, x0 - pad_space)
        x1 = min(mov_x - 1, x1 + pad_space)

        pix_mask = np.zeros((mov_y, mov_x), dtype=bool)
        pix_mask[y0:y1 + 1, x0:x1 + 1] = True
        which_pixels[:, cc] = pix_mask.ravel()

        setup = _setup_cell_window(
            profiles, this_cell, y0, y1, x0, x1, min_pix_for_inclusion,
            lambda_prof, lambda_blob, sigma2_ds, p, one_blob,
        )
        n_y_list.append(setup['n_y'])
        n_x_list.append(setup['n_x'])
        rois_list.append(setup['rois'])
        k1_list.append(setup['k1'])
        k2_list.append(setup['k2'])
        norm_factors_list.append(setup['norm_factors'])
        lambdas_list.append(setup['lambdas'])
        cell_index_within_list.append(setup['cell_index_within'])
        rois_scaled_list.append(setup['rois_scaled'])
        operators.append(setup['operators'])
        which_profiles[:, cc] = setup['include']

    # ---- main per-frame loop ----

    tc_seudo = np.full((n_frames, n_cells_analyzed), np.nan)
    tc_lsq = np.full((n_frames, n_cells_analyzed), np.nan)
    lsq_costs = np.full((n_frames, n_cells_analyzed), np.nan)
    bob_costs = np.full((n_frames, n_cells_analyzed), np.nan)
    tc_cells_with_blobs = np.full((n_frames, n_cells_analyzed), np.nan)

    if frame_blocks is None:
        blocks = [(0, n_frames - 1)]
    else:
        blocks = _merge_frame_blocks(frame_blocks)

    n_block_frames = sum(e - s + 1 for s, e in blocks)
    frames_done = 0

    executor = ThreadPoolExecutor(max_workers=n_jobs) if n_jobs and n_jobs > 1 else None
    try:
        for block_start, block_end in blocks:
            # Phase 1 (sequential): the ds_time moving average has a genuine
            # order dependency (each frame's average depends on the raw
            # pixels of the preceding ds_time-1 frames), but is cheap --
            # just array reads/writes, no optimization. Precompute all of
            # them up front so Phase 2's actual solves are fully independent.
            sliding_window = np.full((mov_y * mov_x, ds_time), np.nan)
            frame_full_by_ff = {}
            for ff in range(block_start, block_end + 1):
                local_ff = ff - block_start
                frame_flat = _get_frame(movie, ff, zero_level).reshape(-1)
                if local_ff < ds_time:
                    sliding_window[:, local_ff] = frame_flat
                else:
                    sliding_window[:, :-1] = sliding_window[:, 1:]
                    sliding_window[:, -1] = frame_flat
                frame_full_by_ff[ff] = np.nanmean(sliding_window, axis=1)

            # Phase 2: solve every (frame, cell) pair in this block. Each is
            # fully independent (own read-only per-cell setup + its own
            # precomputed frame average), so safe to run concurrently --
            # results are written into the shared output arrays only after
            # they come back, sequentially, in the main thread.
            def _task(ff, cc):
                this_frame = frame_full_by_ff[ff][which_pixels[:, cc]]
                return ff, cc, _solve_one_frame_cell(
                    this_frame, rois_list[cc], rois_scaled_list[cc], lambdas_list[cc],
                    norm_factors_list[cc], k1_list[cc], k2_list[cc],
                    n_y_list[cc], n_x_list[cc], one_blob, operators[cc],
                    use_native, native_l_mode, native_nthreads, solver_tol, solver_max_iter,
                )

            work_items = [(ff, cc) for ff in range(block_start, block_end + 1)
                          for cc in range(n_cells_analyzed)]

            if executor is not None:
                results_iter = executor.map(lambda item: _task(*item), work_items)
            else:
                results_iter = (_task(*item) for item in work_items)

            for ff, cc, (tc_lsq_frame, fit_fancy, fit_x, lsq_cost, bob_cost) in results_iter:
                idx_within = cell_index_within_list[cc]
                tc_seudo[ff, cc] = fit_fancy[idx_within]
                tc_cells_with_blobs[ff, cc] = fit_x[idx_within]
                tc_lsq[ff, cc] = tc_lsq_frame[idx_within]
                lsq_costs[ff, cc] = lsq_cost
                bob_costs[ff, cc] = bob_cost

                if cc == n_cells_analyzed - 1:
                    frames_done += 1
                    if verbose and (frames_done % max(1, n_block_frames // 10) == 0):
                        print(f'  frame {frames_done} of {n_block_frames}')
    finally:
        if executor is not None:
            executor.shutdown()

    params = dict(
        p=p, sigma2=sigma2, lambda_blob=lambda_blob, blob_radius=blob_radius,
        ds_time=ds_time, lambda_prof=lambda_prof,
        min_pix_for_inclusion=min_pix_for_inclusion, pad_space=pad_space,
        use_com=use_com, which_cells=which_cells, frame_blocks=blocks,
        use_native=use_native, n_jobs=n_jobs,
    )
    extras = dict(
        lsq_costs=lsq_costs, bob_costs=bob_costs,
        tc_cells_with_blobs=tc_cells_with_blobs, one_blob=one_blob,
        which_pixels=which_pixels, which_profiles=which_profiles,
    )

    return dict(tc=tc_seudo, tc_lsq=tc_lsq, params=params, extras=extras)
