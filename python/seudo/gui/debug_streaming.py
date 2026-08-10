"""Three-panel streaming-pipeline debugger, for visually auditing what
realSEUDOfit actually saw and did on a real run -- built after repeated
real-data findings this session (illumination artifacts, frame-to-frame
unstable Xtemp shapes, under/over-segmentation) that were hard to pin down
from aggregate discovered/matched/unmatched numbers alone.

Panel 1: the actual denoised input frames (avg_frame -- post lookahead-
average + spatial denoise), scrubbable with a slider -- "what did the
detector actually see."
Panel 2: every Xtemp candidate track ever created, in order of creation,
each with its own last-known profile snapshot, creation-frame timestamp,
and actual fate (promoted / merged into an existing cell / dropped via
gap -- captured by hooking _promote_candidate directly during the capture
run, not inferred after the fact).
Panel 3: every promoted Xstab cell, in order of promotion, with its
promotion-frame timestamp.

Data comes from scripts/debug_capture_run.py's HDF5 output (a real
realSEUDOfit run's internals aren't otherwise exposed/retained -- Xtemp
tracks are ephemeral and avg_frame is a local variable, so a dedicated
capture pass is needed before this GUI has anything to show).

Usage:
    python scripts/debug_capture_run.py   # once, writes scripts/output/debug_capture.h5
    python scripts/open_debug_streaming_gui.py
"""

import h5py
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5 import QtCore, QtWidgets

CROP_PAD = 8

# must match scripts/debug_capture_run.py's OUTCOME_* constants
OUTCOME_LABELS = {-1: 'outcome unavailable', 0: 'dropped', 1: 'promoted', 2: 'merged'}


def _crop_bounds(profile, mov_shape, pad=CROP_PAD):
    """Bounding box of a profile's own nonzero footprint, padded, clamped
    to the frame -- falls back to the full frame if the profile is all
    zero (e.g. a degenerate capture)."""
    mov_y, mov_x = mov_shape
    ys, xs = np.nonzero(profile > 0)
    if ys.size == 0:
        return 0, mov_y - 1, 0, mov_x - 1
    y0, y1 = max(0, int(ys.min()) - pad), min(mov_y - 1, int(ys.max()) + pad)
    x0, x1 = max(0, int(xs.min()) - pad), min(mov_x - 1, int(xs.max()) + pad)
    return y0, y1, x0, x1


class DebugStreamingWindow(QtWidgets.QMainWindow):
    def __init__(self, h5_path):
        super().__init__()
        self.setWindowTitle(f'realSEUDOfit debug viewer -- {h5_path}')
        self._load(h5_path)

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QHBoxLayout(central)

        layout.addWidget(self._build_frame_panel())
        layout.addWidget(self._build_xtemp_panel())
        layout.addWidget(self._build_xstab_panel())

        self.resize(1500, 550)

    def _load(self, h5_path):
        with h5py.File(h5_path, 'r') as f:
            self.avg_frames = f['avg_frames'][:]
            self.avg_frame_indices = f['avg_frame_indices'][:]
            self.xtemp_profiles = f['xtemp_profiles'][:]
            self.xtemp_track_ids = f['xtemp_track_ids'][:]
            self.xtemp_creation_frames = f['xtemp_creation_frames'][:]
            self.xtemp_last_seen_frames = f['xtemp_last_seen_frames'][:]
            # older capture files predate these fields -- fall back to
            # sentinels rather than fail to load
            n_tracks = self.xtemp_creation_frames.shape
            self.xtemp_match_counts = f['xtemp_match_counts'][:] if 'xtemp_match_counts' in f else \
                np.full(n_tracks, -1, dtype=np.int64)
            self.xtemp_outcome_codes = f['xtemp_outcome_codes'][:] if 'xtemp_outcome_codes' in f else \
                np.full(n_tracks, -1, dtype=np.int64)
            self.xtemp_outcome_cell_ids = f['xtemp_outcome_cell_ids'][:] if 'xtemp_outcome_cell_ids' in f else \
                np.full(n_tracks, -1, dtype=np.int64)
            self.xtemp_outcome_frames = f['xtemp_outcome_frames'][:] if 'xtemp_outcome_frames' in f else \
                np.full(n_tracks, -1, dtype=np.int64)
            self.xstab_profiles = f['xstab_profiles'][:]  # (Y, X, n_cells)
            self.xstab_promotion_frames = f['xstab_promotion_frames'][:]

        self.mov_shape = self.xstab_profiles.shape[:2] if self.xstab_profiles.shape[2] else self.avg_frames.shape[1:]
        if self.avg_frames.shape[0]:
            lo, hi = np.percentile(self.avg_frames, [1, 99.5])
            self._frame_vmin, self._frame_vmax = float(lo), float(hi)
        else:
            self._frame_vmin, self._frame_vmax = 0.0, 1.0

    # ---- panel 1: denoised input frames --------------------------------

    def _build_frame_panel(self):
        n = self.avg_frames.shape[0]
        panel, canvas, ax, slider, label = self._build_generic_panel(
            'Panel 1: denoised input frames (post lookahead + spatial denoise)', n)
        self._frame_canvas, self._frame_ax, self._frame_slider, self._frame_label = canvas, ax, slider, label
        slider.valueChanged.connect(self._on_frame_slider)
        if n:
            self._draw_frame(0)
        else:
            label.setText('no frames captured')
        return panel

    def _draw_frame(self, i):
        self._frame_ax.clear()
        self._frame_ax.imshow(self.avg_frames[i], cmap='hot', vmin=self._frame_vmin, vmax=self._frame_vmax)
        self._frame_ax.axis('off')
        self._frame_canvas.draw_idle()
        self._frame_label.setText(f'frame {i + 1}/{self.avg_frames.shape[0]}  '
                                   f'(reports frame_index {self.avg_frame_indices[i]})')

    def _on_frame_slider(self, value):
        self._draw_frame(value)

    # ---- panel 2: Xtemp tracks ------------------------------------------

    def _build_xtemp_panel(self):
        n = self.xtemp_profiles.shape[0]
        panel, canvas, ax, slider, label = self._build_generic_panel(
            'Panel 2: Xtemp candidate tracks, in order of creation', n)
        self._xtemp_canvas, self._xtemp_ax, self._xtemp_slider, self._xtemp_label = canvas, ax, slider, label
        slider.valueChanged.connect(self._on_xtemp_slider)
        if n:
            self._draw_xtemp(0)
        else:
            label.setText('no Xtemp tracks captured')
        return panel

    def _draw_xtemp(self, i):
        profile = self.xtemp_profiles[i]
        y0, y1, x0, x1 = _crop_bounds(profile, self.mov_shape)
        self._xtemp_ax.clear()
        self._xtemp_ax.imshow(profile[y0:y1 + 1, x0:x1 + 1], cmap='hot')
        self._xtemp_ax.axis('off')
        self._xtemp_canvas.draw_idle()
        match_count = self.xtemp_match_counts[i]
        match_str = f'{match_count} match(es)' if match_count >= 0 else 'match count unavailable'
        outcome_code = int(self.xtemp_outcome_codes[i])
        outcome_str = OUTCOME_LABELS.get(outcome_code, 'outcome unavailable')
        if outcome_code in (1, 2):
            outcome_str += f' -> cell {self.xtemp_outcome_cell_ids[i]} at frame {self.xtemp_outcome_frames[i]}'
        self._xtemp_label.setText(
            f'track {i + 1}/{self.xtemp_profiles.shape[0]}  (id={self.xtemp_track_ids[i]})\n'
            f'added at frame {self.xtemp_creation_frames[i]}, last matched at frame {self.xtemp_last_seen_frames[i]}\n'
            f'{match_str}  --  {outcome_str}')

    def _on_xtemp_slider(self, value):
        self._draw_xtemp(value)

    # ---- panel 3: Xstab cells --------------------------------------------

    def _build_xstab_panel(self):
        n = self.xstab_profiles.shape[2]
        panel, canvas, ax, slider, label = self._build_generic_panel(
            'Panel 3: Xstab cells, in order of promotion', n)
        self._xstab_canvas, self._xstab_ax, self._xstab_slider, self._xstab_label = canvas, ax, slider, label
        slider.valueChanged.connect(self._on_xstab_slider)
        if n:
            self._draw_xstab(0)
        else:
            label.setText('no Xstab cells captured')
        return panel

    def _draw_xstab(self, i):
        profile = self.xstab_profiles[:, :, i]
        y0, y1, x0, x1 = _crop_bounds(profile, self.mov_shape)
        self._xstab_ax.clear()
        self._xstab_ax.imshow(profile[y0:y1 + 1, x0:x1 + 1], cmap='hot')
        self._xstab_ax.axis('off')
        self._xstab_canvas.draw_idle()
        self._xstab_label.setText(
            f'cell {i + 1}/{self.xstab_profiles.shape[2]}  (cell_id={i})\n'
            f'promoted at frame {self.xstab_promotion_frames[i]}')

    def _on_xstab_slider(self, value):
        self._draw_xstab(value)

    # ---- shared panel scaffolding -----------------------------------------

    def _build_generic_panel(self, title, n_items):
        panel = QtWidgets.QWidget()
        vbox = QtWidgets.QVBoxLayout(panel)

        title_label = QtWidgets.QLabel(title)
        title_label.setWordWrap(True)
        vbox.addWidget(title_label)

        fig = Figure(figsize=(4, 4))
        canvas = FigureCanvas(fig)
        canvas.setMinimumSize(350, 350)
        ax = fig.add_subplot(111)
        vbox.addWidget(canvas)

        info_label = QtWidgets.QLabel('')
        vbox.addWidget(info_label)

        slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        slider.setMinimum(0)
        slider.setMaximum(max(0, n_items - 1))
        slider.setValue(0)
        slider.setEnabled(n_items > 1)
        vbox.addWidget(slider)

        return panel, canvas, ax, slider, info_label
