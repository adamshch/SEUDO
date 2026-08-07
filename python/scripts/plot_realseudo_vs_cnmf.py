"""Field-of-view comparison between realSEUDOfit's discovered cells
(code/demoDataRealSEUDO.mat, produced by run_realseudo_full_movie.py) and the
CNMF ground truth already bundled in code/demoData1.mat -- two color-coded
maps side by side, one neuron = one color, with corresponding (matched)
neurons sharing the same color across both subplots, so agreement/disagreement
is visible directly rather than having to cross-reference a table.

Each pixel is colored by whichever cell has the highest peak-normalized
profile value there (a clean "winner take all" segmentation, not alpha
blending -- overlapping-cell pixels would otherwise muddy the color and
defeat the point of a per-cell color key).
"""

import sys

import colorsys
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, '.')
sys.path.insert(0, 'scripts')
from seudo.matlab_io import load_profiles  # noqa: E402
from compare_realseudo_to_cnmf import match_rois  # noqa: E402

CNMF_PATH = '../code/demoData1.mat'
DISCOVERED_PATH = '../code/demoDataRealSEUDO.mat'
OUT_PATH = 'scripts/output/realseudo_vs_cnmf_fov.png'
MATCH_MAX_DIST = 10.0
BACKGROUND = (0.08, 0.08, 0.1)


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


def segmentation_map(profiles, colors, threshold=0.05):
    """profiles: (Y, X, nCells). colors: list of nCells RGB tuples.
    Returns an (Y, X, 3) image where each pixel takes the color of whichever
    cell peaks highest there (peak-normalized per cell first, so a dim cell
    isn't always shadowed by a brighter one), or BACKGROUND if no cell's
    normalized value clears `threshold` at that pixel."""
    mov_y, mov_x, n_cells = profiles.shape
    peaks = profiles.max(axis=(0, 1), keepdims=True)
    peaks[peaks == 0] = 1.0
    normed = profiles / peaks

    best_val = normed.max(axis=2)
    best_idx = normed.argmax(axis=2)

    img = np.full((mov_y, mov_x, 3), BACKGROUND, dtype=float)
    color_arr = np.array(colors)
    has_cell = best_val > threshold
    img[has_cell] = color_arr[best_idx[has_cell]]
    return img


def main():
    with h5py.File(CNMF_PATH, 'r') as f:
        cnmf_profiles = load_profiles(f)
    with h5py.File(DISCOVERED_PATH, 'r') as f:
        disc_profiles = load_profiles(f)

    n_cnmf = cnmf_profiles.shape[2]
    n_disc = disc_profiles.shape[2]
    matched_cnmf, matched_disc = match_rois(cnmf_profiles, disc_profiles, MATCH_MAX_DIST)

    n_matched = len(matched_cnmf)
    unmatched_cnmf = [i for i in range(n_cnmf) if i not in matched_cnmf]
    unmatched_disc = [j for j in range(n_disc) if j not in matched_disc]

    palette = distinct_colors(n_matched + len(unmatched_cnmf) + len(unmatched_disc))
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

    cnmf_img = segmentation_map(cnmf_profiles, cnmf_colors)
    disc_img = segmentation_map(disc_profiles, disc_colors)

    fig, axes = plt.subplots(1, 2, figsize=(14, 7), facecolor='white')
    axes[0].imshow(disc_img)
    axes[0].set_title(f'realSEUDOfit discovered ({n_disc} cells, {n_matched} matched to CNMF)')
    axes[0].axis('off')

    axes[1].imshow(cnmf_img)
    axes[1].set_title(f'CNMF ground truth ({n_cnmf} cells, {n_matched} matched)')
    axes[1].axis('off')

    fig.suptitle('Corresponding neurons share the same color between panels; '
                  'unmatched neurons get their own unique color', y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT_PATH, dpi=150)
    print(f'saved comparison figure to {OUT_PATH}')
    print(f'{n_matched}/{n_cnmf} CNMF cells matched, {len(unmatched_disc)}/{n_disc} discovered cells unmatched')


if __name__ == '__main__':
    main()
