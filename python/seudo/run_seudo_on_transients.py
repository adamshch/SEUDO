"""Port of the core idea in @seudo/parallelSEUDO.m: run SEUDO for one or
more cells, each restricted to just the frames spanned by its own detected
transients (transient_info times), instead of the whole movie.

Unlike parallelSEUDO, cells are run sequentially rather than via parfor --
a simplification, not a fidelity gap: each cell's computation is already
independent of the others, so sequential vs. parallel execution doesn't
change the result, only wall-clock time.
"""

import time

import numpy as np

from .estimate import estimate_time_courses_with_seudo


def frame_blocks_for_cell(transient_info, cell_idx):
    """(start, end) 0-indexed inclusive frame ranges spanning cell_idx's
    detected transients (already tPre/tPost padded by compute_transient_info)."""
    times = transient_info[cell_idx]['times']
    return [(int(s), int(e)) for s, e in times]


def run_seudo_restricted_to_transients(
    se,
    which_struct='default',
    which_cells=None,
    verbose=True,
    progress_callback=None,
    frame_progress_callback=None,
    **seudo_params,
):
    """Run estimate_time_courses_with_seudo for each of which_cells (default:
    all cells), restricted to that cell's own detected-transient frames.

    progress_callback, if given, is called as progress_callback(done, total,
    cell_idx) after each cell finishes -- lets a caller (e.g. a GUI) report
    progress without depending on stdout.

    frame_progress_callback, if given, is called as
    frame_progress_callback(col, cell_idx, frames_done, total_frames) as each
    *frame* within the current cell finishes -- finer-grained than
    progress_callback, which only fires once a whole cell is done. Useful for
    a GUI progress bar that shouldn't sit at 0% for the entire duration of a
    single (possibly slow) cell's run.

    Returns a combined result dict (tc/tc_lsq of shape (movF, len(which_cells)))
    and appends it to se.tc_seudo, matching estimate_time_courses_with_seudo's
    own convention.
    """
    tc_struct = se._resolve_tc_struct(which_struct)
    transient_info = tc_struct['transient_info']

    if which_cells is None:
        which_cells = list(range(se.n_cells))
    else:
        which_cells = list(which_cells)

    tc = np.full((se.mov_f, len(which_cells)), np.nan)
    tc_lsq = np.full((se.mov_f, len(which_cells)), np.nan)
    per_cell_results = {}

    for col, cell_idx in enumerate(which_cells):
        blocks = frame_blocks_for_cell(transient_info, cell_idx)
        if not blocks:
            if verbose:
                print(f'cell {cell_idx}: no transients, skipping')
        else:
            if verbose:
                n_frames_total = sum(e - s + 1 for s, e in blocks)
                print(f'cell {cell_idx}: running SEUDO on {len(blocks)} '
                      f'transient(s), {n_frames_total} frames...')

            t0 = time.time()
            result = estimate_time_courses_with_seudo(
                se.movie, se.profiles, which_cells=[cell_idx], zero_level=se.zero_level,
                frame_blocks=blocks, verbose=False,
                progress_callback=(
                    None if frame_progress_callback is None
                    else lambda done, total, col=col, cell_idx=cell_idx: (
                        frame_progress_callback(col, cell_idx, done, total))
                ),
                **seudo_params,
            )
            elapsed = time.time() - t0
            tc[:, col] = result['tc'][:, 0]
            tc_lsq[:, col] = result['tc_lsq'][:, 0]
            per_cell_results[cell_idx] = result

            if verbose:
                print(f'cell {cell_idx}: done in {elapsed:0.2f}s')

        if progress_callback is not None:
            progress_callback(col + 1, len(which_cells), cell_idx)

    combined = dict(
        tc=tc,
        tc_lsq=tc_lsq,
        params=dict(which_cells=which_cells, which_struct=which_struct, **seudo_params),
        extras=dict(per_cell_results=per_cell_results),
    )
    se.tc_seudo.append(combined)
    return combined
