"""Port of @seudo/autoClassifyTransients.m -- 'corr' method only (the
default, and what demo.m and the classification GUI actually use: the
GUI's contamination-severity scatter plot and 'residual ratio' sort/title
options are exactly this function's 'resRatios' output).

The 'lbq' method (2D Ljung-Box autocorrelation test, the paper's more
elaborate classifier) is NOT yet ported.

For each transient, this fits a linear combination of [this cell's
profile, every other cell's profile] to the transient's (blurred) spatial
shape, then measures how much residual is left after including this
cell's own profile in the fit vs. leaving it out ('resRatios' -- high
means the shape isn't well explained by this cell, i.e. likely
contamination) and the plain correlation between the (unblurred) shape and
this cell's profile, which is thresholded to produce an automatic
true/false classification.
"""

import numpy as np
from scipy.ndimage import gaussian_filter

from .constants import VAL_FALSE, VAL_TRUE
from .stats import correlation_vector_matrix


def _blur_stack(images, sigma):
    """images: (Y, X, T). Blur each frame spatially, independently."""
    if images.size == 0 or sigma == 0:
        return images
    return gaussian_filter(images, sigma=(sigma, sigma, 0), mode='nearest')


def _crop_to_window(arr, win_sub_y, win_sub_x):
    return arr[win_sub_y, :][:, win_sub_x]


def auto_classify_transients(
    se,
    which_struct='default',
    which_cells=None,
    overwrite=False,
    save_results=True,
    profile_field_name='shapes',
    ignore_zeros=False,
    minimum_window=True,
    blur_radius=1.0,
    blur_profiles_for_fitting=False,
    corr_thresh=0.4,
    weighted_trans=False,
    res_ratio_mean=True,
):
    tc_struct = se._resolve_tc_struct(which_struct)
    if tc_struct.get('transient_info') is None:
        se.compute_transient_info(which_struct)
        tc_struct = se._resolve_tc_struct(which_struct)

    transient_info = tc_struct['transient_info']
    n_cells = se.n_cells

    if which_cells is None:
        which_cells = list(range(n_cells))

    if save_results:
        has_class = {
            cc: not np.all(np.isnan(transient_info[cc]['classification']))
            for cc in which_cells
        }
        if any(has_class.values()) and not overwrite:
            n_bad = sum(has_class.values())
            raise ValueError(
                f'classification exists for {n_bad} of {len(which_cells)} cells, '
                f'to overwrite call with overwrite=True'
            )

    params = dict(
        which_cells=which_cells, overwrite=overwrite, save_results=save_results,
        profile_field_name=profile_field_name, ignore_zeros=ignore_zeros,
        minimum_window=minimum_window, blur_radius=blur_radius,
        blur_profiles_for_fitting=blur_profiles_for_fitting, method='corr',
        corr_thresh=corr_thresh, weighted_trans=weighted_trans, res_ratio_mean=res_ratio_mean,
    )

    results = []

    for cell_id in which_cells:
        ti = transient_info[cell_id]
        n_trans = ti['times'].shape[0]

        if n_trans == 0:
            metric = dict(cfrac=np.nan, rfrac=np.nan, corrs=np.array([]),
                          classification=np.array([]), params=params)
            results.append(metric)
            if save_results:
                ti['classification'] = np.array([])
                ti['auto_class'] = metric
            continue

        tc = tc_struct['tc'][:, cell_id]
        evts = ti['times']
        tc_amps = np.array([np.max(tc[s:e + 1]) for s, e in evts], dtype=float)

        y0, y1, x0, x1 = ti['window']
        x_prof = se.profiles[y0:y1 + 1, x0:x1 + 1, cell_id]
        other_idx = [c for c in range(n_cells) if c != cell_id]
        x_other = se.profiles[y0:y1 + 1, x0:x1 + 1, other_idx]

        if minimum_window:
            win_sub_y = np.any(x_prof > 0, axis=1)
            win_sub_x = np.any(x_prof > 0, axis=0)
        else:
            win_sub_y = np.ones(x_prof.shape[0], dtype=bool)
            win_sub_x = np.ones(x_prof.shape[1], dtype=bool)

        x_prof = _crop_to_window(x_prof, win_sub_y, win_sub_x)
        f_fits = ti[profile_field_name][win_sub_y, :, :][:, win_sub_x, :]
        x_other = x_other[win_sub_y, :, :][:, win_sub_x, :]

        w_y, w_x = x_prof.shape

        if blur_profiles_for_fitting:
            x_prof_for_fit = gaussian_filter(x_prof, sigma=blur_radius, mode='nearest')
            x_other_for_fit = _blur_stack(x_other, blur_radius)
        else:
            x_prof_for_fit = x_prof
            x_other_for_fit = x_other

        x_prof_for_fit = x_prof_for_fit.reshape(-1)
        x_other_for_fit = x_other_for_fit.reshape(w_y * w_x, -1)
        x_other_for_fit = x_other_for_fit[:, x_other_for_fit.sum(axis=0) != 0]

        x_prof_flat = x_prof.reshape(-1).astype(float)
        f_fits_nofilt = f_fits.reshape(w_y * w_x, n_trans)
        x_other_flat = x_other.reshape(w_y * w_x, -1)
        x_other_flat = x_other_flat[:, x_other_flat.sum(axis=0) != 0]  # noqa: F841 (parity w/ MATLAB, unused downstream)

        f_fits_blurred = _blur_stack(f_fits, blur_radius).reshape(w_y * w_x, n_trans)

        design = np.column_stack([x_prof_for_fit, x_other_for_fit])
        fit_vals, *_ = np.linalg.lstsq(design, f_fits_blurred, rcond=None)

        res_evt = f_fits_blurred - x_other_for_fit @ fit_vals[1:, :]
        res = f_fits_blurred - design @ fit_vals
        res_img = res.reshape(w_y, w_x, n_trans)

        res_norms = np.sum(res ** 2, axis=0)
        res_no_cell = np.sum(res_evt ** 2, axis=0)
        with np.errstate(divide='ignore', invalid='ignore'):
            res_ratios = res_norms / res_no_cell

        metric = dict(
            resNorms=res_norms, resNoCell=res_no_cell, resRatios=res_ratios,
            residuals=res_img, residualsWithoutSource=res_evt.reshape(w_y, w_x, n_trans),
        )

        x_prof_corr = x_prof_flat.copy()
        f_fits_corr = f_fits_nofilt.astype(float).copy()
        if ignore_zeros:
            zero_pix = x_prof_flat == 0
            x_prof_corr[zero_pix] = np.nan
            f_fits_corr[zero_pix, :] = np.nan

        corrs = correlation_vector_matrix(x_prof_corr, f_fits_corr)
        metric['corrs'] = corrs

        new_classification = np.full(n_trans, VAL_FALSE)
        new_classification[corrs >= corr_thresh] = VAL_TRUE
        is_false = new_classification == VAL_FALSE
        is_true = ~is_false

        res_den2 = res_no_cell[is_true].sum()
        res_den = res_no_cell[is_false].sum()
        weights = (tc_amps / tc_amps.sum()) if weighted_trans else np.full(n_trans, 1.0 / n_trans)

        cfrac = rfrac = rfrac2 = rfrac3 = 0.0
        for kk in range(n_trans):
            tmp_res = res_img[:, :, kk].copy()
            if ignore_zeros:
                tmp_res.reshape(-1)[x_prof_flat == 0] = np.nan
            res_sq_sum = np.nansum(tmp_res ** 2)

            if is_false[kk]:
                cfrac += weights[kk]
                rd = res_den
            else:
                rd = res_den2

            if res_ratio_mean:
                contrib = res_sq_sum / rd
                rfrac3_contrib = res_sq_sum / res_no_cell.sum()
            else:
                denom = np.nansum(res_evt[:, kk] ** 2)
                contrib = res_sq_sum / denom
                rfrac3_contrib = contrib

            if is_false[kk]:
                rfrac += contrib
            else:
                rfrac2 += contrib
            rfrac3 += rfrac3_contrib

        if not res_ratio_mean:
            if is_false.sum() > 0:
                rfrac /= is_false.sum()
            if is_true.sum() > 0:
                rfrac2 /= is_true.sum()
            rfrac3 /= n_trans

        metric['cfrac'] = cfrac
        metric['rfrac'] = rfrac
        metric['rfrac2'] = rfrac2
        metric['rfrac3'] = rfrac3
        metric['fullVal'] = np.sum(tc_amps * res_ratios) / np.sum(tc_amps)
        metric['classification'] = new_classification
        metric['params'] = params

        results.append(metric)

        if save_results:
            ti['classification'] = new_classification
            ti['auto_class'] = metric

    return results
