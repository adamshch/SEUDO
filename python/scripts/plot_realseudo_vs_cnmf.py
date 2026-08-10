"""Comparison plots between realSEUDOfit's discovered cells
(code/demoDataRealSEUDO.mat, produced by run_realseudo_full_movie.py) and the
CNMF ground truth already bundled in code/demoData1.mat.

Three figures, saved under scripts/output/:
- realseudo_vs_cnmf_contours.png: side-by-side spatial comparison, ROI
  boundaries drawn as colored CONTOUR OUTLINES over a shared grayscale
  backdrop (the standard CNMF/CaImAn convention) rather than filled/blended
  regions -- individual ROIs stay legible even where several overlap.
  Corresponding (matched) neurons share the same color across both panels.
- realseudo_vs_cnmf_traces.png: for every matched cell, the CNMF ground-truth
  time course and realSEUDO's discovered time course plotted together
  (discovered rescaled to the CNMF trace's own peak for display -- the two
  are on different, unrelated calibration scales -- with the applied factor
  shown in the row label), stacked in one figure with a vertical offset per
  cell, matching each cell's contour-plot color.
- realseudo_rois_grid.png: every discovered ROI's own profile, cropped to
  its own footprint and rendered alone, one panel per cell -- addresses ROIs
  that are hard to make out in the combined spatial view.
"""

import math
import os
import sys

import colorsys
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, '.')
sys.path.insert(0, 'scripts')
from seudo.geometry import compute_roi_coms  # noqa: E402
from seudo.matlab_io import load_profiles, load_time_courses  # noqa: E402
from compare_realseudo_to_cnmf import match_rois  # noqa: E402

CNMF_PATH = '../code/demoData1.mat'
DISCOVERED_PATH = '../code/demoDataRealSEUDO.mat'
OUT_DIR = 'scripts/output'
MATCH_MAX_DIST = 10.0
CONTOUR_LEVEL = 0.2  # fraction of each cell's own peak
BACKGROUND_GAIN = 1.4


def distinct_colors(n):
    """n perceptually well-separated colors via golden-angle hue stepping --
    avoids nearby indices landing on visually similar hues, unlike a plain
    linspace over [0, 1)."""
    golden = 0.6180339887498949
    hue = 0.0
    colors = []
    for _ in range(n):
        colors.append(colorsys.hsv_to_rgb(hue, 0.85, 0.95))
        hue = (hue + golden) % 1.0
    return colors


def peak_normalize(profiles):
    peaks = profiles.max(axis=(0, 1), keepdims=True)
    peaks[peaks == 0] = 1.0
    return profiles / peaks


def grayscale_backdrop(normed_profiles, gain=BACKGROUND_GAIN):
    """A dim grayscale "where the cells generally are" reference image
    (sum of every cell's own peak-normalized profile), shared by both
    contour panels so a viewer has one consistent visual anchor."""
    bg = np.clip(normed_profiles.sum(axis=2) * gain, 0.0, 1.0)
    return np.stack([bg, bg, bg], axis=-1)


def draw_contours(ax, normed_profiles, colors, level=CONTOUR_LEVEL):
    n_cells = normed_profiles.shape[2]
    for i in range(n_cells):
        prof = normed_profiles[:, :, i]
        if prof.max() <= 0:
            continue
        ax.contour(prof, levels=[level], colors=[colors[i]], linewidths=1.3)


def build_color_assignment(n_cnmf, n_disc, matched_cnmf, matched_disc):
    unmatched_cnmf = [i for i in range(n_cnmf) if i not in matched_cnmf]
    unmatched_disc = [j for j in range(n_disc) if j not in matched_disc]
    palette = distinct_colors(len(matched_cnmf) + len(unmatched_cnmf) + len(unmatched_disc))

    cnmf_colors = [None] * n_cnmf
    disc_colors = [None] * n_disc
    pos = 0
    for cnmf_idx, (disc_idx, _dist) in matched_cnmf.items():
        cnmf_colors[cnmf_idx] = palette[pos]
        disc_colors[disc_idx] = palette[pos]
        pos += 1
    for cnmf_idx in unmatched_cnmf:
        cnmf_colors[cnmf_idx] = palette[pos]
        pos += 1
    for disc_idx in unmatched_disc:
        disc_colors[disc_idx] = palette[pos]
        pos += 1
    return cnmf_colors, disc_colors, unmatched_cnmf, unmatched_disc


def plot_contour_comparison(cnmf_profiles, disc_profiles, cnmf_colors, disc_colors,
                             n_matched, out_path):
    cnmf_normed = peak_normalize(cnmf_profiles)
    disc_normed = peak_normalize(disc_profiles)
    backdrop = grayscale_backdrop(cnmf_normed)

    fig, axes = plt.subplots(1, 2, figsize=(14, 7), facecolor='white')
    axes[0].imshow(backdrop)
    draw_contours(axes[0], disc_normed, disc_colors)
    axes[0].set_title(f'realSEUDOfit discovered ({disc_profiles.shape[2]} cells, {n_matched} matched to CNMF)')
    axes[0].axis('off')

    axes[1].imshow(backdrop)
    draw_contours(axes[1], cnmf_normed, cnmf_colors)
    axes[1].set_title(f'CNMF ground truth ({cnmf_profiles.shape[2]} cells, {n_matched} matched)')
    axes[1].axis('off')

    fig.suptitle('ROI boundaries at 20% of each cell\'s own peak; corresponding neurons share a color', y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f'saved contour comparison to {out_path}')


def plot_trace_comparison(T, discovered_tc, matched_cnmf, disc_colors, first_detected_frame, out_path):
    pairs = sorted(matched_cnmf.items())
    if not pairs:
        print('no matched cells -- skipping trace comparison plot')
        return

    n = len(pairs)
    fig, ax = plt.subplots(figsize=(16, max(4, 0.55 * n)), facecolor='white')

    offset_unit = 0.0
    for cnmf_idx, (disc_idx, _dist) in pairs:
        t_true = T[:, cnmf_idx]
        finite = t_true[np.isfinite(t_true)]
        if finite.size:
            offset_unit = max(offset_unit, float(np.nanmax(np.abs(finite))))
    offset_unit = offset_unit * 1.4 or 1.0

    yticks, ylabels = [], []
    for row, (cnmf_idx, (disc_idx, _dist)) in enumerate(pairs):
        y0 = row * offset_unit
        t_true = T[:, cnmf_idx]
        t_disc = discovered_tc[:, disc_idx]
        valid = np.isfinite(t_disc)

        true_peak = float(np.nanmax(np.abs(t_true))) if np.any(np.isfinite(t_true)) else 0.0
        disc_peak = float(np.nanmax(np.abs(t_disc[valid]))) if valid.any() else 0.0
        scale = (true_peak / disc_peak) if disc_peak > 0 else 1.0

        ax.plot(np.arange(len(t_true)), y0 + t_true, color='0.55', linewidth=0.8)
        start = first_detected_frame.get(disc_idx) or 0
        frames = np.arange(len(t_disc))
        ax.plot(frames[valid], y0 + scale * t_disc[valid], color=disc_colors[disc_idx], linewidth=1.0)
        ax.axvline(start, color=disc_colors[disc_idx], linewidth=0.6, linestyle=':', alpha=0.6)

        yticks.append(y0)
        ylabels.append(f'cnmf {cnmf_idx} / disc {disc_idx}  (x{scale:.3g})')

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=7)
    ax.set_xlabel('frame')
    ax.set_title('CNMF ground truth (gray) vs. realSEUDOfit discovered (colored, rescaled to CNMF peak) -- '
                  'dotted line marks the promotion frame')
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f'saved trace comparison to {out_path}')


def plot_roi_grid(disc_profiles, disc_colors, matched_disc_to_cnmf, out_path, n_cols=6, pad=8):
    n_cells = disc_profiles.shape[2]
    if n_cells == 0:
        print('no discovered cells -- skipping ROI grid')
        return

    _coms, bounds = compute_roi_coms(disc_profiles)
    mov_y, mov_x = disc_profiles.shape[:2]
    n_rows = math.ceil(n_cells / n_cols)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(2.1 * n_cols, 2.1 * n_rows), facecolor='white')
    axes = np.atleast_2d(axes)

    for i in range(n_rows * n_cols):
        ax = axes[i // n_cols, i % n_cols]
        ax.axis('off')
        if i >= n_cells:
            continue

        y0, y1, x0, x1 = bounds[i]
        if np.any(np.isnan([y0, y1, x0, x1])):
            continue
        y0 = max(0, int(y0) - pad)
        y1 = min(mov_y - 1, int(y1) + pad)
        x0 = max(0, int(x0) - pad)
        x1 = min(mov_x - 1, int(x1) + pad)
        crop = disc_profiles[y0:y1 + 1, x0:x1 + 1, i]

        ax.imshow(crop, cmap='hot')
        label = f'cell {i}'
        if i in matched_disc_to_cnmf:
            label += f' (= cnmf {matched_disc_to_cnmf[i]})'
        ax.set_title(label, fontsize=8, color=disc_colors[i])

    fig.suptitle('Every discovered ROI, cropped to its own footprint', y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    print(f'saved per-ROI grid to {out_path}')


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    with h5py.File(CNMF_PATH, 'r') as f:
        cnmf_profiles = load_profiles(f)
        T = load_time_courses(f, 'T')
    with h5py.File(DISCOVERED_PATH, 'r') as f:
        disc_profiles = load_profiles(f)
        discovered_tc = load_time_courses(f, 'T')
        first_detected_frame = {i: (v if v >= 0 else None)
                                 for i, v in enumerate(np.asarray(f['firstDetectedFrame']).ravel())}

    n_cnmf = cnmf_profiles.shape[2]
    n_disc = disc_profiles.shape[2]
    matched_cnmf, matched_disc = match_rois(cnmf_profiles, disc_profiles, MATCH_MAX_DIST)
    matched_disc_to_cnmf = {disc_idx: cnmf_idx for cnmf_idx, (disc_idx, _dist) in matched_cnmf.items()}

    cnmf_colors, disc_colors, _unmatched_cnmf, unmatched_disc = build_color_assignment(
        n_cnmf, n_disc, matched_cnmf, matched_disc)

    plot_contour_comparison(cnmf_profiles, disc_profiles, cnmf_colors, disc_colors,
                             len(matched_cnmf), f'{OUT_DIR}/realseudo_vs_cnmf_contours.png')
    plot_trace_comparison(T, discovered_tc, matched_cnmf, disc_colors, first_detected_frame,
                           f'{OUT_DIR}/realseudo_vs_cnmf_traces.png')
    plot_roi_grid(disc_profiles, disc_colors, matched_disc_to_cnmf, f'{OUT_DIR}/realseudo_rois_grid.png')

    print(f'\n{len(matched_cnmf)}/{n_cnmf} CNMF cells matched, {len(unmatched_disc)}/{n_disc} discovered cells unmatched')


if __name__ == '__main__':
    main()
