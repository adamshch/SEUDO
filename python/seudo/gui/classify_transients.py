"""Python/PyQt5 port of seudoClassifyTransients.m -- core manual-labeling
loop, including the autoClassifyTransients-derived contamination-severity
scatter plot, residual-ratio sort/title options, the click-on-field-of-view
cell picker, a Gaussian blur control for the thumbnails (std dev in
pixels), continuous scrolling through a cell's transients (a QScrollArea
replacing MATLAB's page-by-page slider), running SEUDO restricted to
detected-transient frames, and horizontal-only zoom/pan on the time course
(see conversation/memory for what's still deferred: clustering/amplitude
sort orders, exhaustive-click mode, the right-click detail-viewer sub-GUI,
alternate colormaps).

Layout: a full-height left sidebar (classification counts, single-cell
profile, clickable field-of-view picture, contamination-severity scatter,
sort/title/blur controls, save/load, run-SEUDO panel) next to a right
column (time course with zoom/pan, scrollable thumbnail grid), with the
cell navigation bar along the bottom.

For responsiveness, classification/artifact changes recolor just the
affected transient's thumbnail border / time-course segment / scatter
point in place (see _update_transient_colors) instead of rebuilding every
plot from scratch -- a full refresh() is only needed when the cell changes
or the underlying data changes (new transients loaded, SEUDO re-run).

Usage:
    se.compute_transient_info('default', ...)   # must be run first
    app = QtWidgets.QApplication([])
    win = ClassifyTransientsWindow(se, 'default')
    win.show()
    app.exec_()
"""

import os

import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5 import QtCore, QtWidgets
from scipy.ndimage import binary_dilation, gaussian_filter

from ..auto_classify import auto_classify_transients
from ..classification_io import save_classification, load_classification
from ..constants import classification_color, cycle_classification, pick_color
from .._native import NATIVE_AVAILABLE
from ..run_seudo_on_transients import run_seudo_restricted_to_transients

SORT_ORDERS = ['time', 'correlation', 'residual ratio']
TITLE_MODES = ['transient id', 'correlation', 'residual ratio']

FOV_HIGHLIGHT_COLOR = (0.0, 1.0, 1.0)  # cellHighlight in the MATLAB code
SEUDO_OVERLAY_COLOR = (0.6, 0.1, 0.8)


def _scale_seudo_for_display(tc, seudo_tc, valid_mask):
    """Display-only rescaling factor for the SEUDO overlay: SEUDO's raw
    output can sit on a very different amplitude scale than the plotted
    default/LSQ time course (different profile normalization, sigma2/lambda
    choices, ...), which makes the overlay hard to compare by eye even when
    it's doing the right thing. Match peak-absolute-height over the frames
    SEUDO actually computed, rather than rescaling the underlying data --
    the stored se.tc_seudo values are left untouched."""
    if not np.any(valid_mask):
        return 1.0
    seudo_amp = np.max(np.abs(seudo_tc[valid_mask]))
    if seudo_amp <= 1e-12:
        return 1.0
    default_amp = np.max(np.abs(tc[valid_mask]))
    return default_amp / seudo_amp


def _contiguous_runs(mask):
    """mask: 1D bool array. Returns [(start, end), ...] inclusive index
    ranges of each contiguous True run, so plotted line segments don't
    connect across gaps (e.g. between two unrelated transient blocks)."""
    runs = []
    start = None
    for i, v in enumerate(mask):
        if v and start is None:
            start = i
        elif not v and start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, len(mask) - 1))
    return runs


class _SeudoRunWorker(QtCore.QThread):
    """Runs run_seudo_restricted_to_transients on a background thread so a
    multi-cell run doesn't freeze the GUI's event loop."""

    progress = QtCore.pyqtSignal(int, int, int)  # done, total, cell_idx
    finished_ok = QtCore.pyqtSignal(dict)
    failed = QtCore.pyqtSignal(str)

    def __init__(self, se, which_struct, which_cells, seudo_params, parent=None):
        super().__init__(parent)
        self.se = se
        self.which_struct = which_struct
        self.which_cells = which_cells
        self.seudo_params = seudo_params

    def run(self):
        try:
            result = run_seudo_restricted_to_transients(
                self.se, self.which_struct, which_cells=self.which_cells, verbose=False,
                progress_callback=lambda done, total, cc: self.progress.emit(done, total, cc),
                **self.seudo_params,
            )
            self.finished_ok.emit(result)
        except Exception as exc:  # noqa: BLE001 -- surface any failure to the GUI, don't crash the thread silently
            self.failed.emit(str(exc))


class ClassifyTransientsWindow(QtWidgets.QMainWindow):
    def __init__(self, se, which_struct='default', n_trans_x=5, parent=None):
        super().__init__(parent)

        self.se = se
        self.which_struct = which_struct
        self.tc_struct = se._resolve_tc_struct(which_struct)
        if self.tc_struct.get('transient_info') is None:
            raise ValueError(
                'compute_transient_info() must be run on this struct before '
                'opening the classification window'
            )

        self.n_trans_x = n_trans_x
        self.this_cell = 0
        self.sort_order = 'time'
        self.title_mode = 'transient id'
        self.blur_sigma = 0.0
        self._ax_to_trans = {}
        self._trans_to_ax = {}
        self._tc_transient_lines = {}
        self._tc_transient_labels = {}
        self._scatter_points = {}
        self._fov_backdrop_rgb = None
        self._last_click_yx = None
        self._last_click_count = 0

        # horizontal-only zoom/pan state for the time-course plot (port of
        # seudo.plotTransients' `set(h,'motion','horizontal','enable','on')`)
        self._tc_xlim = None  # (xmin, xmax) currently shown; None = full range
        self._tc_pan_start = None  # (start_xdata, start_xlim) while dragging

        self.setWindowTitle(
            f"{getattr(se, 'name', 'untitled')}: transient classification for "
            f"{which_struct}, [{se.mov_y} x {se.mov_x}] pixels, {se.mov_f} frames, "
            f"{se.n_cells} cells"
        )

        self._recompute_auto_class()
        self.pixel_lookup = self._build_pixel_lookup()
        self._build_ui()
        self._go_to_cell(0)

    def _recompute_auto_class(self):
        # port of handles.transStats = se.autoClassifyTransients('default','saveResults',0)
        # in seudoClassifyTransients_OpeningFcn -- computed once (not per cell-switch), and
        # crucially with save_results=False so it never touches the manual classification.
        self.auto_class = auto_classify_transients(self.se, self.which_struct, save_results=False)

    def _build_pixel_lookup(self):
        # port of handles.pLookup{y}{x} built from
        # pMap = convn(se.profiles>0, ones(3,3,1), 'same') in the OpeningFcn --
        # for each pixel, which cells have profile support within 1px of it.
        mask = self.se.profiles > 0
        dilated = binary_dilation(mask, structure=np.ones((3, 3, 1), dtype=bool))
        lookup = {}
        ys, xs = np.nonzero(np.any(dilated, axis=2))
        for y, x in zip(ys, xs):
            lookup[(int(y), int(x))] = sorted(np.flatnonzero(dilated[y, x, :]).tolist())
        return lookup

    # ---- UI construction ----

    def _build_ui(self):
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        outer = QtWidgets.QVBoxLayout(central)

        main_row = QtWidgets.QHBoxLayout()

        left_widget = QtWidgets.QWidget()
        left_widget.setLayout(self._build_left_panel())
        left_widget.setMaximumWidth(330)

        # the left column (counts, profile, FOV picture, scatter plot,
        # sort/title/blur controls, save/load, SEUDO run panel) spans the
        # full window height and can be taller than the window -- scroll it,
        # rather than letting it get silently squashed
        left_scroll = QtWidgets.QScrollArea()
        left_scroll.setWidget(left_widget)
        left_scroll.setWidgetResizable(True)
        left_scroll.setMaximumWidth(350)
        left_scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        main_row.addWidget(left_scroll, 0)

        main_row.addLayout(self._build_right_panel(), 1)
        outer.addLayout(main_row, 1)

        outer.addLayout(self._build_cell_nav_bar())

        QtWidgets.QShortcut(QtCore.Qt.Key_PageDown, self, self._scroll_thumbnails_page_down)
        QtWidgets.QShortcut(QtCore.Qt.Key_PageUp, self, self._scroll_thumbnails_page_up)
        QtWidgets.QShortcut(QtCore.Qt.Key_A, self, self._toggle_artifact_shortcut)
        for digit in range(10):
            QtWidgets.QShortcut(
                QtCore.Qt.Key_0 + digit, self, lambda d=digit: self.blur_spin.setValue(d)
            )
        QtWidgets.QShortcut(QtCore.Qt.Key_QuoteLeft, self, lambda: self.blur_spin.setValue(0))

    def _build_counts_panel(self):
        group = QtWidgets.QGroupBox('Classification counts')
        grid = QtWidgets.QGridLayout(group)

        def make_label():
            lbl = QtWidgets.QLabel()
            font = lbl.font()
            font.setPointSize(11)
            font.setBold(True)
            lbl.setFont(font)
            return lbl

        self.count_labels = {}
        keys = ('true', 'false', 'mixed', 'unclassified', 'artifact', 'non_artifact')
        for i, key in enumerate(keys):
            lbl = make_label()
            self.count_labels[key] = lbl
            grid.addWidget(lbl, i // 2, i % 2)

        return group

    def _build_left_panel(self):
        col = QtWidgets.QVBoxLayout()

        col.addWidget(self._build_counts_panel())

        self.profile_fig = Figure(figsize=(3.0, 2.4))
        self.profile_canvas = FigureCanvas(self.profile_fig)
        self.profile_canvas.setMinimumSize(280, 220)
        col.addWidget(self.profile_canvas)

        self.fov_fig = Figure(figsize=(3.2, 3.2))
        self.fov_canvas = FigureCanvas(self.fov_fig)
        self.fov_canvas.setMinimumSize(300, 300)
        self.fov_canvas.mpl_connect('button_press_event', self._on_fov_click)
        col.addWidget(self.fov_canvas)

        self.scatter_fig = Figure(figsize=(2.8, 1.8))
        self.scatter_canvas = FigureCanvas(self.scatter_fig)
        self.scatter_canvas.setMinimumSize(280, 180)
        col.addWidget(self.scatter_canvas)

        self.artifact_checkbox = QtWidgets.QCheckBox('Is artifact')
        self.artifact_checkbox.stateChanged.connect(self._on_artifact_changed)
        col.addWidget(self.artifact_checkbox)

        sort_row = QtWidgets.QHBoxLayout()
        sort_row.addWidget(QtWidgets.QLabel('Sort by:'))
        self.sort_combo = QtWidgets.QComboBox()
        self.sort_combo.addItems(SORT_ORDERS)
        self.sort_combo.currentTextChanged.connect(self._on_sort_changed)
        sort_row.addWidget(self.sort_combo)
        col.addLayout(sort_row)

        title_row = QtWidgets.QHBoxLayout()
        title_row.addWidget(QtWidgets.QLabel('Title:'))
        self.title_combo = QtWidgets.QComboBox()
        self.title_combo.addItems(TITLE_MODES)
        self.title_combo.currentTextChanged.connect(self._on_title_changed)
        title_row.addWidget(self.title_combo)
        col.addLayout(title_row)

        blur_row = QtWidgets.QHBoxLayout()
        blur_row.addWidget(QtWidgets.QLabel('Blur (px std):'))
        self.blur_spin = QtWidgets.QDoubleSpinBox()
        self.blur_spin.setRange(0.0, 20.0)
        self.blur_spin.setSingleStep(0.5)
        self.blur_spin.setValue(0.0)
        self.blur_spin.valueChanged.connect(self._on_blur_changed)
        blur_row.addWidget(self.blur_spin)
        col.addLayout(blur_row)

        save_btn = QtWidgets.QPushButton('Save classification...')
        save_btn.clicked.connect(self._on_save)
        col.addWidget(save_btn)

        load_btn = QtWidgets.QPushButton('Load classification...')
        load_btn.clicked.connect(self._on_load)
        col.addWidget(load_btn)

        col.addWidget(self._build_run_seudo_panel())

        col.addStretch(1)
        return col

    def _build_run_seudo_panel(self):
        # port of the core idea in @seudo/parallelSEUDO.m: run SEUDO restricted
        # to each cell's own detected-transient frames (see run_seudo_on_transients.py)
        group = QtWidgets.QGroupBox('Run SEUDO')
        form = QtWidgets.QFormLayout(group)

        self.seudo_sigma2_spin = QtWidgets.QDoubleSpinBox()
        self.seudo_sigma2_spin.setDecimals(4)
        self.seudo_sigma2_spin.setRange(0.0001, 10.0)
        self.seudo_sigma2_spin.setSingleStep(0.001)
        self.seudo_sigma2_spin.setValue(0.0020)
        form.addRow('sigma2:', self.seudo_sigma2_spin)

        self.seudo_lambda_blob_spin = QtWidgets.QDoubleSpinBox()
        self.seudo_lambda_blob_spin.setRange(0.0, 1000.0)
        self.seudo_lambda_blob_spin.setSingleStep(1.0)
        self.seudo_lambda_blob_spin.setValue(10.0)
        form.addRow('lambdaBlob:', self.seudo_lambda_blob_spin)

        self.seudo_blob_radius_spin = QtWidgets.QDoubleSpinBox()
        self.seudo_blob_radius_spin.setRange(0.1, 20.0)
        self.seudo_blob_radius_spin.setSingleStep(0.5)
        self.seudo_blob_radius_spin.setValue(3.0)
        form.addRow('blobRadius:', self.seudo_blob_radius_spin)

        self.seudo_ds_time_spin = QtWidgets.QSpinBox()
        self.seudo_ds_time_spin.setRange(1, 50)
        self.seudo_ds_time_spin.setValue(3)
        form.addRow('dsTime:', self.seudo_ds_time_spin)

        self.seudo_pad_space_spin = QtWidgets.QSpinBox()
        self.seudo_pad_space_spin.setRange(0, 100)
        self.seudo_pad_space_spin.setValue(5)
        form.addRow('padSpace:', self.seudo_pad_space_spin)

        cpu_count = os.cpu_count() or 1
        self.seudo_n_jobs_spin = QtWidgets.QSpinBox()
        self.seudo_n_jobs_spin.setRange(1, max(1, cpu_count))
        self.seudo_n_jobs_spin.setValue(min(8, cpu_count))
        form.addRow('parallel jobs:', self.seudo_n_jobs_spin)

        native_note = ('native accelerator: available' if NATIVE_AVAILABLE
                        else 'native accelerator: not built (using pure Python -- '
                             '"parallel jobs" won\'t help much, see seudo/_native/build_native.sh)')
        form.addRow(QtWidgets.QLabel(native_note))

        self.run_seudo_one_btn = QtWidgets.QPushButton('Run SEUDO (this cell)')
        self.run_seudo_one_btn.clicked.connect(self._run_seudo_this_cell)
        form.addRow(self.run_seudo_one_btn)

        self.run_seudo_all_btn = QtWidgets.QPushButton('Run SEUDO (all cells)')
        self.run_seudo_all_btn.clicked.connect(self._run_seudo_all_cells)
        form.addRow(self.run_seudo_all_btn)

        self.seudo_status_label = QtWidgets.QLabel('')
        self.seudo_status_label.setWordWrap(True)
        form.addRow(self.seudo_status_label)

        return group

    def _build_right_panel(self):
        col = QtWidgets.QVBoxLayout()

        col.addLayout(self._build_time_course_section())

        self.thumb_fig = Figure()
        self.thumb_canvas = FigureCanvas(self.thumb_fig)
        self.thumb_canvas.mpl_connect('button_press_event', self._on_thumbnail_click)

        self.thumb_scroll = QtWidgets.QScrollArea()
        self.thumb_scroll.setWidget(self.thumb_canvas)
        self.thumb_scroll.setWidgetResizable(False)
        self.thumb_scroll.setMinimumSize(400, 300)
        col.addWidget(self.thumb_scroll, 1)

        return col

    def _build_time_course_section(self):
        col = QtWidgets.QVBoxLayout()

        self.tc_fig = Figure(figsize=(9, 3))
        self.tc_canvas = FigureCanvas(self.tc_fig)
        self.tc_canvas.setMinimumHeight(220)
        self.tc_canvas.mpl_connect('scroll_event', self._on_tc_scroll)
        self.tc_canvas.mpl_connect('button_press_event', self._on_tc_pan_press)
        self.tc_canvas.mpl_connect('motion_notify_event', self._on_tc_pan_motion)
        self.tc_canvas.mpl_connect('button_release_event', self._on_tc_pan_release)
        col.addWidget(self.tc_canvas)

        tc_controls_row = QtWidgets.QHBoxLayout()
        tc_controls_row.addWidget(QtWidgets.QLabel(
            'scroll to zoom, drag to pan (horizontal only)'
        ))
        tc_controls_row.addStretch(1)
        reset_zoom_btn = QtWidgets.QPushButton('Reset zoom')
        reset_zoom_btn.clicked.connect(self._reset_tc_zoom)
        tc_controls_row.addWidget(reset_zoom_btn)
        col.addLayout(tc_controls_row)

        return col

    def _build_cell_nav_bar(self):
        row = QtWidgets.QHBoxLayout()

        row.addWidget(QtWidgets.QLabel('Cell:'))

        self.cell_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.cell_slider.setMinimum(1)
        self.cell_slider.setMaximum(max(1, self.se.n_cells))
        self.cell_slider.valueChanged.connect(lambda v: self._go_to_cell(v - 1))
        row.addWidget(self.cell_slider, 1)

        self.cell_spin = QtWidgets.QSpinBox()
        self.cell_spin.setMinimum(1)
        self.cell_spin.setMaximum(max(1, self.se.n_cells))
        self.cell_spin.valueChanged.connect(lambda v: self._go_to_cell(v - 1))
        row.addWidget(self.cell_spin)

        prev_btn = QtWidgets.QPushButton('◀ prev unclassified')
        prev_btn.clicked.connect(self.go_to_prev_unclassified)
        row.addWidget(prev_btn)

        next_btn = QtWidgets.QPushButton('next unclassified ▶')
        next_btn.clicked.connect(self.go_to_next_unclassified)
        row.addWidget(next_btn)

        return row

    # ---- cell navigation ----

    def _go_to_cell(self, cell_idx):
        cell_idx = int(min(max(0, cell_idx), self.se.n_cells - 1))
        self.this_cell = cell_idx
        self.cell_slider.blockSignals(True)
        self.cell_spin.blockSignals(True)
        self.cell_slider.setValue(cell_idx + 1)
        self.cell_spin.setValue(cell_idx + 1)
        self.cell_slider.blockSignals(False)
        self.cell_spin.blockSignals(False)
        self._tc_xlim = None  # switching cells resets the time-course zoom/pan
        self.refresh()
        self.thumb_scroll.verticalScrollBar().setValue(0)

    def _cell_transient_info(self):
        return self.tc_struct['transient_info'][self.this_cell]

    def _scroll_thumbnails_page_down(self):
        bar = self.thumb_scroll.verticalScrollBar()
        bar.setValue(bar.value() + bar.pageStep())

    def _scroll_thumbnails_page_up(self):
        bar = self.thumb_scroll.verticalScrollBar()
        bar.setValue(bar.value() - bar.pageStep())

    def _unclassified_cells(self):
        info = self.tc_struct['transient_info']
        result = []
        for i, ti in enumerate(info):
            cls = ti['classification']
            if cls.size > 0 and np.all(np.isnan(cls)) and not ti['is_artifact']:
                result.append(i)
        return result

    def go_to_next_unclassified(self):
        candidates = [c for c in self._unclassified_cells() if c > self.this_cell]
        if candidates:
            self._go_to_cell(min(candidates))

    def go_to_prev_unclassified(self):
        candidates = [c for c in self._unclassified_cells() if c < self.this_cell]
        if candidates:
            self._go_to_cell(max(candidates))

    # ---- classification interaction ----

    def _cell_auto_class(self):
        return self.auto_class[self.this_cell]

    def _plot_order(self):
        ti = self._cell_transient_info()
        n_trans = ti['times'].shape[0]
        metric = self._cell_auto_class()
        if self.sort_order == 'correlation' and metric['corrs'].size == n_trans:
            return list(np.argsort(-metric['corrs']))
        if self.sort_order == 'residual ratio' and metric['resRatios'].size == n_trans:
            return list(np.argsort(metric['resRatios']))
        return list(range(n_trans))

    def _on_thumbnail_click(self, event):
        if event.inaxes not in self._ax_to_trans:
            return
        trans_idx = self._ax_to_trans[event.inaxes]
        ti = self._cell_transient_info()
        ti['classification'][trans_idx] = cycle_classification(ti['classification'][trans_idx])
        # only this one transient's color changed anywhere -- no need to
        # rebuild the FOV picture, the whole time course, or the whole grid
        self._update_transient_colors([trans_idx])

    def _on_fov_click(self, event):
        # port of chooseCellFromClick: repeated clicks at the same pixel
        # cycle through overlapping cells at that location
        if event.xdata is None or event.ydata is None:
            return
        x = int(min(max(0, round(event.xdata)), self.se.mov_x - 1))
        y = int(min(max(0, round(event.ydata)), self.se.mov_y - 1))

        cell_list = self.pixel_lookup.get((y, x), [])
        if not cell_list:
            return

        if self._last_click_yx == (y, x):
            self._last_click_count = (self._last_click_count % len(cell_list)) + 1
        else:
            self._last_click_yx = (y, x)
            self._last_click_count = 1

        new_cell = cell_list[self._last_click_count - 1]
        if new_cell != self.this_cell:
            self._go_to_cell(new_cell)

    def _on_artifact_changed(self, _state):
        ti = self._cell_transient_info()
        ti['is_artifact'] = self.artifact_checkbox.isChecked()
        # affects every transient's color for this cell, but nothing about
        # what's plotted (shapes, time course data, window) changes
        self._update_transient_colors(range(ti['times'].shape[0]))

    def _toggle_artifact_shortcut(self):
        self.artifact_checkbox.setChecked(not self.artifact_checkbox.isChecked())

    def _on_sort_changed(self, text):
        self.sort_order = text
        self._draw_thumbnails()

    def _on_title_changed(self, text):
        self.title_mode = text
        self._draw_thumbnails()

    def _on_blur_changed(self, value):
        self.blur_sigma = value
        self._draw_thumbnails()

    def _on_save(self):
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save classification')
        if path:
            save_classification(self.tc_struct, path)

    def _on_load(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Load classification')
        if path:
            payload = load_classification(path)
            self.tc_struct['transient_info'] = payload['transient_info']
            # transient boundaries/shapes may have changed -- auto_class (computed
            # against the previous transient_info) is now stale, so redo it
            self._recompute_auto_class()
            self.refresh()

    # ---- running SEUDO ----

    def _collect_seudo_params(self):
        return dict(
            sigma2=self.seudo_sigma2_spin.value(),
            lambda_blob=self.seudo_lambda_blob_spin.value(),
            blob_radius=self.seudo_blob_radius_spin.value(),
            ds_time=self.seudo_ds_time_spin.value(),
            pad_space=self.seudo_pad_space_spin.value(),
            n_jobs=self.seudo_n_jobs_spin.value(),
        )

    def _run_seudo(self, which_cells):
        if getattr(self, '_seudo_worker', None) is not None and self._seudo_worker.isRunning():
            return

        params = self._collect_seudo_params()
        self.run_seudo_one_btn.setEnabled(False)
        self.run_seudo_all_btn.setEnabled(False)
        self.seudo_status_label.setText(f'Running SEUDO on {len(which_cells)} cell(s)...')

        worker = _SeudoRunWorker(self.se, self.which_struct, which_cells, params, self)
        worker.progress.connect(self._on_seudo_progress)
        worker.finished_ok.connect(self._on_seudo_finished)
        worker.failed.connect(self._on_seudo_failed)
        self._seudo_worker = worker
        worker.start()

    def _run_seudo_this_cell(self):
        self._run_seudo([self.this_cell])

    def _run_seudo_all_cells(self):
        self._run_seudo(list(range(self.se.n_cells)))

    def _on_seudo_progress(self, done, total, cell_idx):
        self.seudo_status_label.setText(f'SEUDO: {done}/{total} cell(s) done (last: cell {cell_idx + 1})')

    def _on_seudo_finished(self, result):
        self.run_seudo_one_btn.setEnabled(True)
        self.run_seudo_all_btn.setEnabled(True)
        self.seudo_status_label.setText(f'SEUDO done ({result["tc"].shape[1]} cell(s)).')
        # only the time-course overlay depends on tc_seudo -- profiles, shapes,
        # and classification are all untouched by running SEUDO
        self._draw_time_course()

    def _on_seudo_failed(self, message):
        self.run_seudo_one_btn.setEnabled(True)
        self.run_seudo_all_btn.setEnabled(True)
        self.seudo_status_label.setText(f'SEUDO failed: {message}')

    def _latest_seudo_tc_for_cell(self, cell_idx):
        if not self.se.tc_seudo:
            return None
        result = self.se.tc_seudo[-1]
        which_cells = result.get('params', {}).get('which_cells')
        if which_cells is None or cell_idx not in which_cells:
            return None
        col = list(which_cells).index(cell_idx)
        tc = result['tc'][:, col]
        return None if np.all(np.isnan(tc)) else tc

    # ---- drawing ----

    def refresh(self):
        """Full redraw of everything -- only needed when the current cell
        changes or the underlying data changes (new transients loaded).
        For classification/artifact color changes on the current cell, use
        _update_transient_colors instead; it's dramatically cheaper."""
        ti = self._cell_transient_info()
        self.artifact_checkbox.blockSignals(True)
        self.artifact_checkbox.setChecked(bool(ti['is_artifact']))
        self.artifact_checkbox.blockSignals(False)

        self._draw_profile()
        self._draw_fov()
        self._draw_time_course()
        self._draw_scatter()
        self._draw_thumbnails()
        self._update_counts()

    def _update_transient_colors(self, trans_indices):
        """Recolor the given transients' thumbnail border / time-course
        segment / scatter point in place -- no plot is cleared or rebuilt,
        so this is fast even when a cell has hundreds of transients."""
        ti = self._cell_transient_info()
        touched_thumb = touched_tc = touched_scatter = False

        for trans_idx in trans_indices:
            color = (pick_color('artifact') if ti['is_artifact']
                     else classification_color(ti['classification'][trans_idx]))

            ax = self._trans_to_ax.get(trans_idx)
            if ax is not None:
                for spine in ax.spines.values():
                    spine.set_color(color)
                touched_thumb = True

            line = self._tc_transient_lines.get(trans_idx)
            if line is not None:
                line.set_color(color)
                touched_tc = True

            label = self._tc_transient_labels.get(trans_idx)
            if label is not None:
                label.set_color(color)
                touched_tc = True

            point = self._scatter_points.get(trans_idx)
            if point is not None:
                point.set_color(color)
                touched_scatter = True

        if touched_thumb:
            self.thumb_canvas.draw_idle()
        if touched_tc:
            self.tc_canvas.draw_idle()
        if touched_scatter:
            self.scatter_canvas.draw_idle()

        self._update_counts()

    def _draw_profile(self):
        ti = self._cell_transient_info()
        y0, y1, x0, x1 = ti['window']
        prof = self.se.profiles[y0:y1 + 1, x0:x1 + 1, self.this_cell]

        self.profile_fig.clf()
        ax = self.profile_fig.add_subplot(111)
        ax.imshow(prof, cmap='gray')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f'source {self.this_cell + 1} profile', fontsize=9)
        self.profile_canvas.draw_idle()

    def _fov_backdrop(self):
        # the "all profiles" composite backdrop only depends on se.profiles,
        # never on which cell is selected -- compute it once and reuse
        if self._fov_backdrop_rgb is None:
            profiles = self.se.profiles.astype(float)
            per_cell_max = profiles.max(axis=(0, 1), keepdims=True)
            per_cell_max = np.where(per_cell_max > 0, per_cell_max, 1.0)
            prof_max = (profiles / per_cell_max).max(axis=2)
            prof_max_peak = prof_max.max()
            if prof_max_peak > 0:
                prof_max = prof_max / prof_max_peak
            self._fov_backdrop_rgb = np.stack([prof_max] * 3, axis=-1)
        return self._fov_backdrop_rgb

    def _draw_fov(self):
        # port of the axesAllProfiles image blend in setNewCell: a dimmed
        # grayscale composite of every cell's (peak-normalized) profile,
        # with the current cell highlighted in teal and its analysis
        # window outlined. Click to jump to a cell (chooseCellFromClick).
        backdrop = self._fov_backdrop()

        this_profile = self.se.profiles[:, :, self.this_cell].astype(float)
        tp_peak = this_profile.max()
        if tp_peak > 0:
            this_profile = this_profile / tp_peak
        highlight = np.stack([this_profile * 0, this_profile * 1.0, this_profile * 0.5], axis=-1)

        im = 1 - (1 - backdrop * 0.7) * (1 - highlight)
        im = np.clip(im, 0, 1)

        self.fov_fig.clf()
        ax = self.fov_fig.add_subplot(111)
        ax.imshow(im)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f'all {self.se.n_cells} profiles', fontsize=10)

        ti = self._cell_transient_info()
        y0, y1, x0, x1 = ti['window']
        ax.plot(
            [x0 - 0.5, x1 + 0.5, x1 + 0.5, x0 - 0.5, x0 - 0.5],
            [y0 - 0.5, y0 - 0.5, y1 + 0.5, y1 + 0.5, y0 - 0.5],
            '-', linewidth=1.5, color=FOV_HIGHLIGHT_COLOR,
        )
        self.fov_canvas.draw_idle()

    def _draw_time_course(self):
        ti = self._cell_transient_info()
        tc = self.tc_struct['tc'][:, self.this_cell]
        self.tc_fig.clf()
        ax = self.tc_fig.add_subplot(111)
        ax.plot(tc, color=(0.9, 0.9, 0.9), linewidth=0.5)

        self._tc_transient_lines = {}
        self._tc_transient_labels = {}
        for tt in range(ti['times'].shape[0]):
            s, e = ti['times'][tt]
            color = classification_color(ti['classification'][tt])
            line, = ax.plot(range(s, e + 1), tc[s:e + 1], color=color)
            self._tc_transient_lines[tt] = line

            # label each transient with its number (matching the thumbnail
            # grid's default "transient id" title), placed at its peak
            peak_local = int(np.argmax(tc[s:e + 1]))
            peak_frame = s + peak_local
            peak_value = tc[peak_frame]
            label = ax.annotate(
                str(tt + 1), xy=(peak_frame, peak_value), xytext=(0, 3),
                textcoords='offset points', fontsize=7, ha='center', va='bottom', color=color,
            )
            self._tc_transient_labels[tt] = label

        seudo_tc = self._latest_seudo_tc_for_cell(self.this_cell)
        if seudo_tc is not None:
            valid = ~np.isnan(seudo_tc)
            scale = _scale_seudo_for_display(tc, seudo_tc, valid)
            seudo_tc_display = seudo_tc * scale
            label = 'SEUDO' if np.isclose(scale, 1.0) else f'SEUDO (×{scale:0.2f})'
            first = True
            for s, e in _contiguous_runs(valid):
                ax.plot(range(s, e + 1), seudo_tc_display[s:e + 1], '-', color=SEUDO_OVERLAY_COLOR,
                         linewidth=1.5, label=label if first else None)
                first = False
            ax.legend(fontsize=7, loc='upper right')

        ax.set_title(f'source {self.this_cell + 1}', fontsize=10)
        for spine in ('top', 'right'):
            ax.spines[spine].set_visible(False)

        if self._tc_xlim is not None:
            ax.set_xlim(*self._tc_xlim)
        else:
            ax.set_xlim(0, len(tc))
            self._tc_xlim = (0, len(tc))

        self.tc_fig.tight_layout()
        self.tc_canvas.draw_idle()

    # ---- time-course horizontal zoom/pan ----

    def _tc_axes(self):
        return self.tc_fig.axes[0] if self.tc_fig.axes else None

    def _on_tc_scroll(self, event):
        ax = self._tc_axes()
        if ax is None or event.inaxes != ax or event.xdata is None:
            return
        n_frames = self.se.mov_f
        cur_left, cur_right = ax.get_xlim()
        width = cur_right - cur_left
        scale = 0.8 if event.button == 'up' else 1.25
        new_width = max(5.0, min(width * scale, float(n_frames)))
        rel = (event.xdata - cur_left) / width if width else 0.5
        new_left = event.xdata - rel * new_width
        new_right = new_left + new_width
        if new_left < 0:
            new_left, new_right = 0.0, new_width
        elif new_right > n_frames:
            new_right = float(n_frames)
            new_left = new_right - new_width
        self._tc_xlim = (new_left, new_right)
        ax.set_xlim(new_left, new_right)
        self.tc_canvas.draw_idle()

    def _on_tc_pan_press(self, event):
        ax = self._tc_axes()
        if ax is None or event.inaxes != ax or event.button != 1:
            return
        self._tc_pan_start = (event.xdata, ax.get_xlim())

    def _on_tc_pan_motion(self, event):
        if self._tc_pan_start is None or event.xdata is None:
            return
        ax = self._tc_axes()
        if ax is None:
            return
        start_xdata, start_xlim = self._tc_pan_start
        n_frames = self.se.mov_f
        width = min(start_xlim[1] - start_xlim[0], float(n_frames))
        dx = start_xdata - event.xdata
        new_left = start_xlim[0] + dx
        new_right = new_left + width
        if new_left < 0:
            new_left, new_right = 0.0, width
        elif new_right > n_frames:
            new_right = float(n_frames)
            new_left = new_right - width
        self._tc_xlim = (new_left, new_right)
        ax.set_xlim(new_left, new_right)
        self.tc_canvas.draw_idle()

    def _on_tc_pan_release(self, _event):
        self._tc_pan_start = None

    def _reset_tc_zoom(self):
        self._tc_xlim = None
        self._draw_time_course()

    def _draw_scatter(self):
        # port of the axesCorrelation scatter plot in setNewCell (correlation
        # vs. autoClassifyTransients' resRatios, i.e. "contamination severity")
        ti = self._cell_transient_info()
        metric = self._cell_auto_class()
        self.scatter_fig.clf()
        ax = self.scatter_fig.add_subplot(111)

        self._scatter_points = {}
        n_trans = ti['times'].shape[0]
        if n_trans > 0 and metric['corrs'].size == n_trans and metric['resRatios'].size == n_trans:
            for tt in range(n_trans):
                point, = ax.plot(metric['corrs'][tt], metric['resRatios'][tt], 'o',
                                  color=classification_color(ti['classification'][tt]), markersize=4)
                self._scatter_points[tt] = point
        ax.set_xlim(-1, 1)
        ax.set_xlabel('correlation', fontsize=8)
        ax.set_ylabel('contam. severity', fontsize=8)
        ax.tick_params(labelsize=7)
        self.scatter_fig.tight_layout()
        self.scatter_canvas.draw_idle()

    def _thumbnail_title(self, trans_idx):
        metric = self._cell_auto_class()
        if self.title_mode == 'correlation' and metric['corrs'].size > trans_idx:
            return f"{metric['corrs'][trans_idx]:0.3f}"
        if self.title_mode == 'residual ratio' and metric['resRatios'].size > trans_idx:
            return f"{metric['resRatios'][trans_idx]:0.3f}"
        return str(trans_idx + 1)

    def _draw_thumbnails(self):
        # all of the current cell's transients are laid out in one tall grid
        # (n_trans_x columns) inside a QScrollArea -- scroll (wheel, drag, or
        # PageUp/PageDown) to move through them, rather than paging screen by
        # screen as in the MATLAB GUI.
        ti = self._cell_transient_info()
        self.thumb_fig.clf()
        self._ax_to_trans = {}
        self._trans_to_ax = {}

        n_trans = ti['times'].shape[0]
        if n_trans == 0:
            self.thumb_fig.set_size_inches(4, 3)
            self.thumb_canvas.setFixedSize(400, 300)
            ax = self.thumb_fig.add_subplot(111)
            ax.text(0.5, 0.5, 'This cell has no transients', ha='center', va='center')
            ax.axis('off')
            self.thumb_canvas.draw_idle()
            return

        plot_order = self._plot_order()
        n_cols = self.n_trans_x
        n_rows = -(-n_trans // n_cols)  # ceil div

        cell_px = 150
        dpi = 100
        self.thumb_fig.set_dpi(dpi)
        self.thumb_fig.set_size_inches(n_cols * cell_px / dpi, n_rows * cell_px / dpi)
        self.thumb_canvas.setFixedSize(n_cols * cell_px, n_rows * cell_px)

        axes = self.thumb_fig.subplots(n_rows, n_cols, squeeze=False)
        for i, ax in enumerate(axes.flat):
            if i >= n_trans:
                ax.axis('off')
                continue

            trans_idx = plot_order[i]
            shape = ti['shapes'][:, :, trans_idx]
            if self.blur_sigma > 0:
                shape = gaussian_filter(shape, sigma=self.blur_sigma)
            ax.imshow(shape, cmap='gray')
            ax.set_xticks([])
            ax.set_yticks([])

            color = pick_color('artifact') if ti['is_artifact'] else classification_color(ti['classification'][trans_idx])
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_color(color)
                spine.set_linewidth(3)

            ax.set_title(self._thumbnail_title(trans_idx), fontsize=8)
            self._ax_to_trans[ax] = trans_idx
            self._trans_to_ax[trans_idx] = ax

        # tight_layout()'s iterative constraint solving gets noticeably slow
        # with the large subplot counts a busy cell can have (100+); fixed
        # margins sized for the fixed 150px-per-cell grid look just as good.
        # top/bottom are fractions of the *total* figure height, which grows
        # with n_rows here (unlike left/right, since n_cols is fixed) -- so
        # they're computed from an absolute pixel margin, not a flat fraction,
        # or a tall grid ends up with a huge wasted gap above the first row.
        total_height_px = n_rows * cell_px
        top = 1 - 15 / total_height_px
        bottom = 5 / total_height_px
        self.thumb_fig.subplots_adjust(left=0.02, right=0.98, top=top, bottom=bottom, wspace=0.15, hspace=0.6)
        self.thumb_canvas.draw_idle()

    def _update_counts(self):
        info = self.tc_struct['transient_info']
        is_artifact = np.array([ti['is_artifact'] for ti in info])
        class_vals = np.concatenate([
            ti['classification'] for ti, art in zip(info, is_artifact) if not art
        ]) if np.any(~is_artifact) else np.array([])

        from .. import constants as C
        n_true = int(np.sum(class_vals == C.VAL_TRUE))
        n_false = int(np.sum(class_vals == C.VAL_FALSE))
        n_unc = int(np.sum(np.isnan(class_vals)))
        n_mix = int(class_vals.size - n_true - n_false - n_unc)
        n_art = int(np.sum(is_artifact))
        n_non_art = int(np.sum(~is_artifact))

        counts = {
            'true': (n_true, pick_color('true')),
            'false': (n_false, pick_color('false')),
            'mixed': (n_mix, pick_color('mixed')),
            'unclassified': (n_unc, pick_color('unclassified')),
            'artifact': (n_art, pick_color('artifact')),
            'non_artifact': (n_non_art, (0, 0, 0)),
        }
        labels = {
            'true': 'True', 'false': 'False', 'mixed': 'Mixed',
            'unclassified': 'Unclassified', 'artifact': 'Artifact', 'non_artifact': 'Non-artifact',
        }
        for key, (count, color) in counts.items():
            r, g, b = (int(255 * c) for c in color)
            self.count_labels[key].setText(f'{labels[key]}: {count}')
            self.count_labels[key].setStyleSheet(f'color: rgb({r},{g},{b})')
