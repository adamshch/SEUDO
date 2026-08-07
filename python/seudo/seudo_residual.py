"""The "SEUDO residual" auto-classification criterion: for each transient,
fit its spatial shape (a single image -- see transients.compute_transient_info's
'shapes' field, not a full movie) two ways, restricted to the cell's own
SEUDO window (same profile/blob-basis setup as estimate.py, just solved
once per transient instead of once per movie frame):

    coef_lsq   -- plain least squares against known (overlapping) profiles
    coef_seudo -- the nonneg weighted-L1 sparse fit (native or Python),
                  which can divert amplitude into the unmodeled "blob"
                  basis instead of this cell's profile

    fraction = (coef_lsq - coef_seudo) / coef_lsq

A real transient, well explained by this cell's own profile, should barely
move between the two fits (fraction near 0). Contamination -- absorbed by
the blob rather than attributed to this cell -- pulls coef_seudo well below
coef_lsq (fraction near 1).

Reuses estimate.py's per-cell window setup and single-fit solve directly
(_setup_cell_window, _solve_one_frame_cell) so this always agrees with what
"Run SEUDO" itself would compute, rather than re-deriving the same math.

Note: auto_classify.py's 'seudo_residual' criterion classifies *high*
fraction as true (per explicit user direction) -- the opposite of the
"low fraction = real" reading above, and the opposite convention from the
'res_ratio' criterion. This module only computes the fraction; the
classification direction lives in auto_classify.py.
"""

import numpy as np

from .blob import make_seudo_blob
from .estimate import _setup_cell_window, _solve_one_frame_cell
from ._native import NATIVE_AVAILABLE as _NATIVE_AVAILABLE

SEUDO_RESIDUAL_DEFAULTS = dict(
    sigma2=0.0020, lambda_blob=10.0, blob_radius=3.0, lambda_prof=0.0, p=1e-5,
    min_pix_for_inclusion=1, solver_tol=0.01, solver_max_iter=1000,
    use_native='auto', native_l_mode=2, native_nthreads=1, blob_spacing=1.0,
)


def compute_seudo_residual_fractions(se, cell_id, ti, **seudo_kwargs):
    """(coef_lsq - coef_seudo) / coef_lsq for each of cell_id's transients in
    ti (a transient_info entry, see compute_transient_info). Returns an
    array of length ti['times'].shape[0]; nan where coef_lsq is ~0 (no
    meaningful fraction to report)."""
    params = dict(SEUDO_RESIDUAL_DEFAULTS)
    params.update(seudo_kwargs)

    use_native = params['use_native']
    if use_native == 'auto':
        use_native = _NATIVE_AVAILABLE
    elif use_native and not _NATIVE_AVAILABLE:
        raise RuntimeError(
            'use_native=True but the compiled native SEUDO extension is not available; '
            'build it with seudo/_native/build_native.sh, or pass use_native=False'
        )

    n_trans = ti['times'].shape[0]
    if n_trans == 0 or ti['shapes'] is None:
        return np.array([])

    y0, y1, x0, x1 = ti['window']
    one_blob = make_seudo_blob(params['blob_radius'])
    setup = _setup_cell_window(
        se.profiles, cell_id, y0, y1, x0, x1, params['min_pix_for_inclusion'],
        params['lambda_prof'], params['lambda_blob'], params['sigma2'], params['p'], one_blob,
    )
    idx = setup['cell_index_within']

    fractions = np.full(n_trans, np.nan)
    for tt in range(n_trans):
        this_frame = ti['shapes'][:, :, tt].reshape(-1)

        tc_lsq_frame, _fit_fancy, fit_x, _lsq_cost, _bob_cost = _solve_one_frame_cell(
            this_frame, setup['rois'], setup['rois_scaled'], setup['lambdas'],
            setup['norm_factors'], setup['k1'], setup['k2'],
            setup['n_y'], setup['n_x'], one_blob, setup['operators'],
            use_native, params['native_l_mode'], params['native_nthreads'],
            params['solver_tol'], params['solver_max_iter'], params['blob_spacing'],
        )

        coef_lsq = tc_lsq_frame[idx]
        coef_seudo = fit_x[idx]
        fractions[tt] = np.nan if coef_lsq == 0 else (coef_lsq - coef_seudo) / coef_lsq

    return fractions
