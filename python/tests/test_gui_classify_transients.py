import types
from unittest import mock

import numpy as np
import pytest

from seudo.core import Seudo
from seudo.constants import VAL_FALSE, VAL_TRUE, VAL_UNC
from seudo.constants import pick_color

pytest.importorskip('PyQt5')

from seudo.gui.classify_transients import ClassifyTransientsWindow, SEUDO_OVERLAY_COLOR  # noqa: E402


def gaussian_patch(shape, center, sigma):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    return np.exp(-(((yy - center[0]) ** 2 + (xx - center[1]) ** 2) / (2 * sigma ** 2)))


def make_two_cell_seudo():
    y, x, f = 20, 40, 500
    prof0 = gaussian_patch((y, x), center=(10, 10), sigma=2.0)
    prof0[prof0 < 0.05 * prof0.max()] = 0.0
    prof1 = gaussian_patch((y, x), center=(10, 30), sigma=2.0)
    prof1[prof1 < 0.05 * prof1.max()] = 0.0

    act0 = np.zeros(f)
    act0[100:110] = 3.0
    act0[300:310] = 3.0
    act1 = np.zeros(f)
    act1[200:210] = 3.0

    rng = np.random.default_rng(3)
    movie = np.zeros((y, x, f))
    for ff in range(f):
        movie[:, :, ff] = (prof0 * act0[ff] + prof1 * act1[ff]
                            + 0.01 * rng.normal(size=(y, x)))

    profiles = np.stack([prof0, prof1], axis=2)
    tc = np.stack([act0, act1], axis=1)

    se = Seudo(movie, profiles, time_courses=tc)
    se.compute_transient_info('default', min_duration=3, blur_radius=1)
    return se


@pytest.fixture
def se():
    return make_two_cell_seudo()


def make_seudo_with_many_transients(n_bursts=30):
    """Cell 0 with many short transients -- enough rows in the thumbnail
    grid to require real scrolling."""
    y, x, f = 20, 20, n_bursts * 30
    prof0 = gaussian_patch((y, x), center=(10, 10), sigma=2.0)
    prof0[prof0 < 0.05 * prof0.max()] = 0.0

    act0 = np.zeros(f)
    for i in range(n_bursts):
        s = i * 30 + 10
        act0[s:s + 8] = 3.0

    rng = np.random.default_rng(6)
    movie = np.zeros((y, x, f))
    for ff in range(f):
        movie[:, :, ff] = prof0 * act0[ff] + 0.01 * rng.normal(size=(y, x))

    profiles = prof0[:, :, np.newaxis]
    tc = act0[:, np.newaxis]

    se = Seudo(movie, profiles, time_courses=tc)
    # bursts are too dense for the default percentile-based auto-threshold
    # (see test_transients.py) -- supply the transient frames explicitly
    transient_frames = (act0 > 0)[:, np.newaxis]
    se.compute_transient_info('default', transient_frames=transient_frames)
    return se


@pytest.fixture
def se_many_transients():
    return make_seudo_with_many_transients()


def test_window_builds_and_shows_first_cell(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        assert win.this_cell == 0
        assert win.se.n_cells == 2
        # cell 0 has 2 transients -> at least one thumbnail axes registered
        assert len(win._ax_to_trans) >= 2
    finally:
        win.close()


def test_click_cycles_classification(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        ti = win._cell_transient_info()
        assert np.all(np.isnan(ti['classification']))

        ax, trans_idx = next(iter(win._ax_to_trans.items()))
        fake_event = types.SimpleNamespace(inaxes=ax)

        win._on_thumbnail_click(fake_event)
        assert ti['classification'][trans_idx] == VAL_FALSE

        # after refresh(), axes objects are rebuilt -- re-fetch the mapping
        # for the same transient index before clicking again
        ax2 = next(a for a, i in win._ax_to_trans.items() if i == trans_idx)
        win._on_thumbnail_click(types.SimpleNamespace(inaxes=ax2))
        assert ti['classification'][trans_idx] == VAL_TRUE
    finally:
        win.close()


def test_artifact_toggle_updates_transient_info(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        assert win._cell_transient_info()['is_artifact'] is False
        win.artifact_checkbox.setChecked(True)
        assert win._cell_transient_info()['is_artifact'] is True
    finally:
        win.close()


def test_cell_navigation_via_spinbox(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        win.cell_spin.setValue(2)  # 1-indexed display -> cell index 1
        assert win.this_cell == 1
    finally:
        win.close()


def test_next_prev_unclassified_navigation(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        # both cells start fully unclassified
        assert win._unclassified_cells() == [0, 1]

        win.go_to_next_unclassified()
        assert win.this_cell == 1

        # classify cell 1's only transient -> no longer "unclassified"
        win._cell_transient_info()['classification'][:] = VAL_TRUE
        assert win._unclassified_cells() == [0]

        win.go_to_prev_unclassified()
        assert win.this_cell == 0
    finally:
        win.close()


def test_save_and_load_round_trip_through_buttons(qapp, se, tmp_path, monkeypatch):
    from PyQt5 import QtWidgets

    win = ClassifyTransientsWindow(se, 'default')
    try:
        ax, trans_idx = next(iter(win._ax_to_trans.items()))
        win._on_thumbnail_click(types.SimpleNamespace(inaxes=ax))
        expected = win._cell_transient_info()['classification'][trans_idx]

        path = str(tmp_path / 'saved.pkl')
        monkeypatch.setattr(
            QtWidgets.QFileDialog, 'getSaveFileName', staticmethod(lambda *a, **k: (path, ''))
        )
        win._on_save()

        # wipe classification, then reload from disk
        win._cell_transient_info()['classification'][:] = VAL_UNC
        monkeypatch.setattr(
            QtWidgets.QFileDialog, 'getOpenFileName', staticmethod(lambda *a, **k: (path, ''))
        )
        win._on_load()

        assert win._cell_transient_info()['classification'][trans_idx] == expected
    finally:
        win.close()


def test_auto_class_computed_on_open_without_touching_classification(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        assert len(win.auto_class) == se.n_cells
        metric0 = win.auto_class[0]
        assert metric0['corrs'].shape[0] == win._cell_transient_info()['times'].shape[0]
        assert metric0['resRatios'].shape[0] == metric0['corrs'].shape[0]
        # save_results=False -- must not have touched the real classification
        assert np.all(np.isnan(win._cell_transient_info()['classification']))
    finally:
        win.close()


def test_auto_classify_button_writes_classification_for_every_cell(qapp, se, monkeypatch):
    from PyQt5 import QtWidgets

    win = ClassifyTransientsWindow(se, 'default')
    try:
        monkeypatch.setattr(QtWidgets.QMessageBox, 'question',
                             staticmethod(lambda *a, **k: QtWidgets.QMessageBox.Yes))
        win.auto_classify_criterion_combo.setCurrentText('correlation')
        win.auto_classify_thresh_spin.setValue(-1.0)  # lenient -- everything should end up true

        win._on_auto_classify_clicked()

        for cell_idx in range(se.n_cells):
            ti = se.tc_default['transient_info'][cell_idx]
            if ti['times'].shape[0] > 0:
                assert np.all(ti['classification'] == VAL_TRUE)
    finally:
        win.close()


def test_auto_classify_res_ratio_criterion_uses_res_ratio_thresh(qapp, se, monkeypatch):
    from PyQt5 import QtWidgets

    win = ClassifyTransientsWindow(se, 'default')
    try:
        monkeypatch.setattr(QtWidgets.QMessageBox, 'question',
                             staticmethod(lambda *a, **k: QtWidgets.QMessageBox.Yes))
        win.auto_classify_criterion_combo.setCurrentText('SEUDO residual')
        win.auto_classify_thresh_spin.setValue(1e6)  # lenient in the res-ratio direction

        win._on_auto_classify_clicked()

        for cell_idx in range(se.n_cells):
            ti = se.tc_default['transient_info'][cell_idx]
            if ti['times'].shape[0] > 0:
                assert np.all(ti['classification'] == VAL_TRUE)
    finally:
        win.close()


def test_auto_classify_declined_confirmation_leaves_classification_untouched(qapp, se, monkeypatch):
    from PyQt5 import QtWidgets

    win = ClassifyTransientsWindow(se, 'default')
    try:
        monkeypatch.setattr(QtWidgets.QMessageBox, 'question',
                             staticmethod(lambda *a, **k: QtWidgets.QMessageBox.No))
        win._on_auto_classify_clicked()

        assert np.all(np.isnan(win._cell_transient_info()['classification']))
    finally:
        win.close()


def test_auto_classify_criterion_switch_resets_default_threshold(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        win.auto_classify_criterion_combo.setCurrentText('SEUDO residual')
        assert win.auto_classify_thresh_spin.value() == pytest.approx(0.5)
        win.auto_classify_criterion_combo.setCurrentText('correlation')
        assert win.auto_classify_thresh_spin.value() == pytest.approx(0.4)
    finally:
        win.close()


def test_sort_by_residual_ratio_orders_transients_ascending(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        win.sort_order = 'residual ratio'
        order = win._plot_order()
        res_ratios = win._cell_auto_class()['resRatios']
        assert list(order) == list(np.argsort(res_ratios))
    finally:
        win.close()


def test_sort_by_correlation_orders_transients_descending(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        win.sort_order = 'correlation'
        order = win._plot_order()
        corrs = win._cell_auto_class()['corrs']
        assert list(order) == list(np.argsort(-corrs))
    finally:
        win.close()


def test_title_mode_changes_thumbnail_labels(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        metric = win._cell_auto_class()

        win.title_mode = 'transient id'
        assert win._thumbnail_title(0) == '1'

        win.title_mode = 'correlation'
        assert win._thumbnail_title(0) == f"{metric['corrs'][0]:0.3f}"

        win.title_mode = 'residual ratio'
        assert win._thumbnail_title(0) == f"{metric['resRatios'][0]:0.3f}"
    finally:
        win.close()


def test_scatter_plot_draws_one_point_per_transient(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        n_trans = win._cell_transient_info()['times'].shape[0]
        lines = win.scatter_fig.axes[0].get_lines()
        assert len(lines) == n_trans
    finally:
        win.close()


def test_load_recomputes_stale_auto_class(qapp, se, tmp_path, monkeypatch):
    from PyQt5 import QtWidgets
    from seudo.classification_io import save_classification

    win = ClassifyTransientsWindow(se, 'default')
    try:
        # save cell 0's current (2-transient) classification state, then
        # smuggle in a transient_info with a different transient count for
        # cell 0 to prove auto_class gets recomputed (not left stale)
        path = str(tmp_path / 'saved.pkl')
        save_classification(win.tc_struct, path)

        payload = {
            'tc': win.tc_struct['tc'],
            'transient_info': [dict(win.tc_struct['transient_info'][0]), win.tc_struct['transient_info'][1]],
        }
        payload['transient_info'][0] = dict(payload['transient_info'][0])
        payload['transient_info'][0]['times'] = payload['transient_info'][0]['times'][:1]
        payload['transient_info'][0]['classification'] = payload['transient_info'][0]['classification'][:1]
        payload['transient_info'][0]['shapes'] = payload['transient_info'][0]['shapes'][:, :, :1]
        import pickle
        with open(path, 'wb') as fh:
            pickle.dump(payload, fh)

        monkeypatch.setattr(
            QtWidgets.QFileDialog, 'getOpenFileName', staticmethod(lambda *a, **k: (path, ''))
        )
        win._on_load()

        assert win.auto_class[0]['corrs'].shape[0] == 1
        assert win._cell_transient_info()['times'].shape[0] == 1
    finally:
        win.close()


def make_two_overlapping_cell_seudo():
    """Two cells whose profiles heavily overlap, so a single field-of-view
    pixel maps to both -- for testing the click-to-cycle behavior."""
    y, x, f = 20, 40, 500
    prof0 = gaussian_patch((y, x), center=(10, 15), sigma=2.0)
    prof0[prof0 < 0.05 * prof0.max()] = 0.0
    prof1 = gaussian_patch((y, x), center=(10, 17), sigma=2.0)
    prof1[prof1 < 0.05 * prof1.max()] = 0.0

    act0 = np.zeros(f)
    act0[100:110] = 3.0
    act1 = np.zeros(f)
    act1[200:210] = 3.0

    rng = np.random.default_rng(5)
    movie = np.zeros((y, x, f))
    for ff in range(f):
        movie[:, :, ff] = prof0 * act0[ff] + prof1 * act1[ff] + 0.01 * rng.normal(size=(y, x))

    profiles = np.stack([prof0, prof1], axis=2)
    tc = np.stack([act0, act1], axis=1)

    se = Seudo(movie, profiles, time_courses=tc)
    se.compute_transient_info('default', min_duration=3, blur_radius=1)
    return se


def test_fov_click_selects_cell_at_that_pixel(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        assert win.this_cell == 0
        win._on_fov_click(types.SimpleNamespace(xdata=30.0, ydata=10.0))  # cell 1's profile peak
        assert win.this_cell == 1
        win._on_fov_click(types.SimpleNamespace(xdata=10.0, ydata=10.0))  # cell 0's profile peak
        assert win.this_cell == 0
    finally:
        win.close()


def test_fov_click_on_background_does_nothing(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        assert win.this_cell == 0
        win._on_fov_click(types.SimpleNamespace(xdata=0.0, ydata=0.0))  # far from either profile
        assert win.this_cell == 0
    finally:
        win.close()


def test_fov_click_cycles_through_overlapping_cells():
    import PyQt5  # noqa: F401 (already gated by module-level importorskip)
    se_overlap = make_two_overlapping_cell_seudo()
    win = ClassifyTransientsWindow(se_overlap, 'default')
    try:
        y, x = 10, 16  # squarely inside the overlap between the two profiles
        cell_list = win.pixel_lookup[(y, x)]
        assert len(cell_list) == 2

        event = types.SimpleNamespace(xdata=float(x), ydata=float(y))
        win._on_fov_click(event)
        first = win.this_cell
        win._on_fov_click(event)
        second = win.this_cell

        assert {first, second} == set(cell_list)
        assert first != second
    finally:
        win.close()


def test_counts_bar_reflects_classification(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        assert win.count_labels['true'].text() == 'True: 0'
        assert win.count_labels['unclassified'].text() == 'Unclassified: 3'

        win._cell_transient_info()['classification'][0] = VAL_TRUE
        win.refresh()

        assert win.count_labels['true'].text() == 'True: 1'
        assert win.count_labels['unclassified'].text() == 'Unclassified: 2'
    finally:
        win.close()


def test_cell_nav_bar_widgets_present(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        assert win.cell_slider.maximum() == se.n_cells
        assert win.cell_spin.maximum() == se.n_cells
    finally:
        win.close()


def test_blur_control_changes_displayed_thumbnail_pixels(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        ax = next(iter(win._ax_to_trans))
        unblurred = ax.images[0].get_array().copy()

        win.blur_spin.setValue(3.0)
        assert win.blur_sigma == 3.0

        ax2 = next(iter(win._ax_to_trans))
        blurred = ax2.images[0].get_array()

        assert not np.array_equal(unblurred, blurred)
        assert blurred.shape == unblurred.shape
        # blurring should reduce the peak (energy spreads out)
        assert blurred.max() <= unblurred.max()
    finally:
        win.close()


def test_blur_zero_leaves_thumbnail_unchanged(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        ax = next(iter(win._ax_to_trans))
        original = ax.images[0].get_array().copy()

        win.blur_spin.setValue(0.0)  # already 0 -- shouldn't error or change anything
        ax2 = next(iter(win._ax_to_trans))
        still = ax2.images[0].get_array()

        assert np.array_equal(original, still)
    finally:
        win.close()


def test_no_pagination_widgets_remain(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        assert not hasattr(win, 'page_slider')
        assert not hasattr(win, 'page_label')
        assert not hasattr(win, 'page')
    finally:
        win.close()


def test_thumbnail_canvas_sized_to_show_all_transients_in_scroll_area(qapp, se_many_transients):
    win = ClassifyTransientsWindow(se_many_transients, 'default')
    try:
        n_trans = win._cell_transient_info()['times'].shape[0]
        assert n_trans == 30
        n_cols = win.n_trans_x
        n_rows = -(-n_trans // n_cols)

        assert win.thumb_canvas.width() == n_cols * 150
        assert win.thumb_canvas.height() == n_rows * 150
        assert win.thumb_scroll.widget() is win.thumb_canvas
        # canvas taller than the scroll area's viewport -> a real vertical scrollbar is needed
        assert win.thumb_canvas.height() > win.thumb_scroll.viewport().height()
    finally:
        win.close()


def test_page_up_down_scroll_thumbnails(qapp, se_many_transients):
    win = ClassifyTransientsWindow(se_many_transients, 'default')
    try:
        bar = win.thumb_scroll.verticalScrollBar()
        bar.setValue(0)
        win._scroll_thumbnails_page_down()
        assert bar.value() > 0
        moved = bar.value()
        win._scroll_thumbnails_page_up()
        assert bar.value() < moved
    finally:
        win.close()


def test_switching_cell_resets_scroll_position(qapp, se_many_transients):
    win = ClassifyTransientsWindow(se_many_transients, 'default')
    try:
        win.thumb_scroll.verticalScrollBar().setValue(200)
        assert win.thumb_scroll.verticalScrollBar().value() > 0

        # se_many_transients only has 1 cell, but switching to "another" cell
        # (itself) still goes through _go_to_cell -> scroll should reset
        win._go_to_cell(0)
        assert win.thumb_scroll.verticalScrollBar().value() == 0
    finally:
        win.close()


def _wait_for_seudo_worker(qapp, win, timeout_ms=30000):
    assert win._seudo_worker.wait(timeout_ms), 'SEUDO worker did not finish in time'
    for _ in range(20):
        qapp.processEvents()


def test_collect_seudo_params_reads_spinboxes(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        win.seudo_sigma2_spin.setValue(0.5)
        win.seudo_lambda_blob_spin.setValue(7.0)
        win.seudo_blob_radius_spin.setValue(2.5)
        win.seudo_ds_time_spin.setValue(2)
        win.seudo_pad_space_spin.setValue(4)
        win.seudo_n_jobs_spin.setValue(3)

        params = win._collect_seudo_params()
        assert params == dict(sigma2=0.5, lambda_blob=7.0, blob_radius=2.5, ds_time=2, pad_space=4, n_jobs=3)
    finally:
        win.close()


def test_run_seudo_this_cell_populates_tc_seudo_and_overlay(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        win._collect_seudo_params = lambda: dict(
            sigma2=0.02, lambda_blob=1.0, blob_radius=1.2, ds_time=1, pad_space=6, solver_max_iter=100,
        )
        assert win.run_seudo_one_btn.isEnabled()

        win._run_seudo_this_cell()
        assert not win.run_seudo_one_btn.isEnabled()
        assert not win.run_seudo_all_btn.isEnabled()

        _wait_for_seudo_worker(qapp, win)

        assert win.run_seudo_one_btn.isEnabled()
        assert win.run_seudo_all_btn.isEnabled()
        assert 'done' in win.seudo_status_label.text().lower()

        assert len(se.tc_seudo) == 1
        assert se.tc_seudo[0]['params']['which_cells'] == [0]
        assert win._latest_seudo_tc_for_cell(0) is not None
    finally:
        win.close()


def test_run_seudo_all_cells_covers_every_cell(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        win._collect_seudo_params = lambda: dict(
            sigma2=0.02, lambda_blob=1.0, blob_radius=1.2, ds_time=1, pad_space=6, solver_max_iter=100,
        )
        win._run_seudo_all_cells()
        _wait_for_seudo_worker(qapp, win)

        result = se.tc_seudo[-1]
        assert result['tc'].shape == (se.mov_f, se.n_cells)
        assert result['params']['which_cells'] == list(range(se.n_cells))
    finally:
        win.close()


def test_run_seudo_failure_reported_in_status_label(qapp, se, monkeypatch):
    import seudo.gui.classify_transients as mod

    def _boom(*args, **kwargs):
        raise RuntimeError('synthetic failure')

    monkeypatch.setattr(mod, 'run_seudo_restricted_to_transients', _boom)

    win = ClassifyTransientsWindow(se, 'default')
    try:
        win._collect_seudo_params = lambda: dict(
            sigma2=0.02, lambda_blob=1.0, blob_radius=1.2, ds_time=1, pad_space=6, solver_max_iter=100,
        )
        win._run_seudo_this_cell()
        _wait_for_seudo_worker(qapp, win)

        assert win.run_seudo_one_btn.isEnabled()
        assert 'failed' in win.seudo_status_label.text().lower()
        assert 'synthetic failure' in win.seudo_status_label.text()
    finally:
        win.close()


def test_contiguous_runs():
    from seudo.gui.classify_transients import _contiguous_runs

    mask = np.array([False, True, True, False, False, True, False, True, True, True])
    assert _contiguous_runs(mask) == [(1, 2), (5, 5), (7, 9)]
    assert _contiguous_runs(np.array([False, False])) == []
    assert _contiguous_runs(np.array([True, True])) == [(0, 1)]


def test_scroll_zooms_in_around_cursor(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        ax = win._tc_axes()
        assert ax.get_xlim() == (0, se.mov_f)

        win._on_tc_scroll(types.SimpleNamespace(inaxes=ax, xdata=100.0, button='up'))

        left, right = ax.get_xlim()
        assert right - left < se.mov_f  # narrower than before
        assert left < 100.0 < right     # still centered around the cursor
        assert win._tc_xlim == (left, right)
    finally:
        win.close()


def test_scroll_zoom_out_clamps_to_full_movie_bounds(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        ax = win._tc_axes()
        # zoom in first so there's room to zoom back out
        win._on_tc_scroll(types.SimpleNamespace(inaxes=ax, xdata=100.0, button='up'))
        for _ in range(20):
            win._on_tc_scroll(types.SimpleNamespace(inaxes=ax, xdata=100.0, button='down'))

        left, right = ax.get_xlim()
        assert left >= 0
        assert right <= se.mov_f
    finally:
        win.close()


def test_scroll_ignored_outside_time_course_axes(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        original_xlim = win._tc_axes().get_xlim()
        win._on_tc_scroll(types.SimpleNamespace(inaxes=None, xdata=100.0, button='up'))
        assert win._tc_axes().get_xlim() == original_xlim
    finally:
        win.close()


def test_pan_drag_shifts_view_without_changing_width(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        ax = win._tc_axes()
        win._on_tc_scroll(types.SimpleNamespace(inaxes=ax, xdata=100.0, button='up'))
        left0, right0 = ax.get_xlim()
        width0 = right0 - left0

        win._on_tc_pan_press(types.SimpleNamespace(inaxes=ax, xdata=100.0, button=1))
        win._on_tc_pan_motion(types.SimpleNamespace(xdata=80.0))  # drag left -> view shifts right...

        left1, right1 = ax.get_xlim()
        assert np.isclose(right1 - left1, width0)  # width preserved (no vertical/scale change)
        assert left1 > left0  # dragging the data left moves the visible window right

        win._on_tc_pan_release(None)
        assert win._tc_pan_start is None
    finally:
        win.close()


def test_pan_clamps_to_full_movie_bounds(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        ax = win._tc_axes()
        win._on_tc_scroll(types.SimpleNamespace(inaxes=ax, xdata=100.0, button='up'))

        win._on_tc_pan_press(types.SimpleNamespace(inaxes=ax, xdata=100.0, button=1))
        win._on_tc_pan_motion(types.SimpleNamespace(xdata=100000.0))  # drag way past the end

        left, right = ax.get_xlim()
        assert left >= 0
        assert right <= se.mov_f
    finally:
        win.close()


def test_reset_zoom_restores_full_range(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        ax = win._tc_axes()
        win._on_tc_scroll(types.SimpleNamespace(inaxes=ax, xdata=100.0, button='up'))
        assert win._tc_axes().get_xlim() != (0, se.mov_f)

        win._reset_tc_zoom()
        assert win._tc_axes().get_xlim() == (0, se.mov_f)
        assert win._tc_xlim == (0, se.mov_f)
    finally:
        win.close()


def test_switching_cell_resets_zoom(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        ax = win._tc_axes()
        win._on_tc_scroll(types.SimpleNamespace(inaxes=ax, xdata=100.0, button='up'))
        assert win._tc_axes().get_xlim() != (0, se.mov_f)

        win._go_to_cell(1)
        assert win._tc_axes().get_xlim() == (0, se.mov_f)
    finally:
        win.close()


def test_classifying_transient_preserves_zoom(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        ax = win._tc_axes()
        win._on_tc_scroll(types.SimpleNamespace(inaxes=ax, xdata=100.0, button='up'))
        zoomed_xlim = win._tc_axes().get_xlim()

        thumb_ax, _ = next(iter(win._ax_to_trans.items()))
        win._on_thumbnail_click(types.SimpleNamespace(inaxes=thumb_ax))

        assert win._tc_axes().get_xlim() == zoomed_xlim
    finally:
        win.close()


def test_seudo_overlay_draws_separate_segments_per_transient_block(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        win._collect_seudo_params = lambda: dict(
            sigma2=0.02, lambda_blob=1.0, blob_radius=1.2, ds_time=1, pad_space=6, solver_max_iter=100,
        )
        win._run_seudo_this_cell()
        _wait_for_seudo_worker(qapp, win)

        ax = win._tc_axes()
        seudo_lines = [ln for ln in ax.get_lines() if tuple(ln.get_color()) == SEUDO_OVERLAY_COLOR
                       or list(ln.get_color()) == list(SEUDO_OVERLAY_COLOR)]
        n_transients = win._cell_transient_info()['times'].shape[0]
        assert len(seudo_lines) == n_transients  # one contiguous segment per transient block

        # exactly one carries the legend label (may include a display-scale factor, e.g. "SEUDO (×0.85)")
        labeled = [ln for ln in seudo_lines if ln.get_label().startswith('SEUDO')]
        assert len(labeled) == 1
    finally:
        win.close()


def test_scale_seudo_for_display_matches_peak_amplitude():
    from seudo.gui.classify_transients import _scale_seudo_for_display

    tc = np.array([0.0, 0.0, 10.0, 0.0, 0.0])
    seudo_tc = np.array([np.nan, np.nan, 2.0, np.nan, np.nan])
    valid = ~np.isnan(seudo_tc)
    assert _scale_seudo_for_display(tc, seudo_tc, valid) == pytest.approx(5.0)


def test_scale_seudo_for_display_no_div_by_zero_when_seudo_is_flat_zero():
    from seudo.gui.classify_transients import _scale_seudo_for_display

    tc = np.array([0.0, 0.0, 10.0, 0.0, 0.0])
    seudo_tc = np.array([np.nan, np.nan, 0.0, np.nan, np.nan])
    valid = ~np.isnan(seudo_tc)
    assert _scale_seudo_for_display(tc, seudo_tc, valid) == 1.0


def test_scale_seudo_for_display_empty_mask_returns_one():
    from seudo.gui.classify_transients import _scale_seudo_for_display

    tc = np.array([1.0, 2.0])
    seudo_tc = np.array([np.nan, np.nan])
    valid = np.zeros(2, dtype=bool)
    assert _scale_seudo_for_display(tc, seudo_tc, valid) == 1.0


def test_seudo_overlay_rescaled_to_match_default_trace_height(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        ti = win._cell_transient_info()
        s0, e0 = ti['times'][0]
        tc = win.tc_struct['tc'][:, win.this_cell]
        true_peak = np.max(np.abs(tc[s0:e0 + 1]))

        # inject a fake SEUDO result on a deliberately different scale (1/10th height)
        fake_tc_col = np.full(se.mov_f, np.nan)
        fake_tc_col[s0:e0 + 1] = tc[s0:e0 + 1] / 10.0
        se.tc_seudo.append({
            'tc': fake_tc_col[:, None],
            'params': {'which_cells': [win.this_cell]},
        })

        win._draw_time_course()

        ax = win._tc_axes()
        seudo_line = next(ln for ln in ax.get_lines() if ln.get_label().startswith('SEUDO'))
        plotted_peak = np.max(np.abs(seudo_line.get_ydata()))

        # displayed height now matches the default trace's height, not the raw 1/10th-scale data
        assert plotted_peak == pytest.approx(true_peak, rel=0.05)
        assert seudo_line.get_label() == 'SEUDO (×10.00)'
    finally:
        win.close()


def test_profile_canvas_shows_current_cells_windowed_profile(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        ax = win.profile_fig.axes[0]
        shown = ax.images[0].get_array()
        ti = win._cell_transient_info()
        y0, y1, x0, x1 = ti['window']
        expected = se.profiles[y0:y1 + 1, x0:x1 + 1, win.this_cell]
        assert np.array_equal(shown, expected)

        win._go_to_cell(1)
        ax = win.profile_fig.axes[0]
        shown1 = ax.images[0].get_array()
        ti1 = win._cell_transient_info()
        y0, y1, x0, x1 = ti1['window']
        expected1 = se.profiles[y0:y1 + 1, x0:x1 + 1, 1]
        assert np.array_equal(shown1, expected1)
    finally:
        win.close()


def test_fov_backdrop_computed_once_and_reused_across_cells(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        backdrop_before = win._fov_backdrop()
        win._go_to_cell(1)
        backdrop_after = win._fov_backdrop()
        assert backdrop_before is backdrop_after  # same cached array, not recomputed
    finally:
        win.close()


def test_thumbnail_click_does_not_rebuild_fov_or_thumbnail_grid(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        with mock.patch.object(win, '_draw_fov', wraps=win._draw_fov) as fov_spy, \
                mock.patch.object(win, '_draw_thumbnails', wraps=win._draw_thumbnails) as thumb_spy, \
                mock.patch.object(win, '_draw_time_course', wraps=win._draw_time_course) as tc_spy, \
                mock.patch.object(win, '_draw_profile', wraps=win._draw_profile) as profile_spy:
            ax, _ = next(iter(win._ax_to_trans.items()))
            win._on_thumbnail_click(types.SimpleNamespace(inaxes=ax))

            fov_spy.assert_not_called()
            thumb_spy.assert_not_called()
            tc_spy.assert_not_called()
            profile_spy.assert_not_called()
    finally:
        win.close()


def test_thumbnail_click_only_recolors_the_clicked_transient(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        items = list(win._ax_to_trans.items())
        assert len(items) >= 2
        (ax0, idx0), (ax1, idx1) = items[0], items[1]

        win._on_thumbnail_click(types.SimpleNamespace(inaxes=ax0))

        color0 = tuple(ax0.spines['top'].get_edgecolor())[:3]
        color1 = tuple(ax1.spines['top'].get_edgecolor())[:3]
        assert np.allclose(color0, pick_color('false'))  # unclassified -> false is the first cycle step
        assert np.allclose(color1, pick_color('unclassified'))  # untouched
    finally:
        win.close()


def test_artifact_toggle_avoids_full_redraw_but_recolors_everything(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        with mock.patch.object(win, '_draw_fov', wraps=win._draw_fov) as fov_spy, \
                mock.patch.object(win, '_draw_thumbnails', wraps=win._draw_thumbnails) as thumb_spy, \
                mock.patch.object(win, '_draw_time_course', wraps=win._draw_time_course) as tc_spy:
            win.artifact_checkbox.setChecked(True)

            fov_spy.assert_not_called()
            thumb_spy.assert_not_called()
            tc_spy.assert_not_called()

        for ax in win._trans_to_ax.values():
            color = tuple(ax.spines['top'].get_edgecolor())[:3]
            assert np.allclose(color, pick_color('artifact'))
    finally:
        win.close()


def test_seudo_finished_only_redraws_time_course(qapp, se):
    win = ClassifyTransientsWindow(se, 'default')
    try:
        win._collect_seudo_params = lambda: dict(
            sigma2=0.02, lambda_blob=1.0, blob_radius=1.2, ds_time=1, pad_space=6, solver_max_iter=100,
        )
        with mock.patch.object(win, '_draw_fov', wraps=win._draw_fov) as fov_spy, \
                mock.patch.object(win, '_draw_thumbnails', wraps=win._draw_thumbnails) as thumb_spy, \
                mock.patch.object(win, '_draw_scatter', wraps=win._draw_scatter) as scatter_spy, \
                mock.patch.object(win, '_draw_profile', wraps=win._draw_profile) as profile_spy:
            win._run_seudo_this_cell()
            _wait_for_seudo_worker(qapp, win)

            fov_spy.assert_not_called()
            thumb_spy.assert_not_called()
            scatter_spy.assert_not_called()
            profile_spy.assert_not_called()

        assert win._latest_seudo_tc_for_cell(0) is not None
    finally:
        win.close()


def test_thumbnail_grid_top_margin_stays_a_constant_pixel_gap(qapp, se, se_many_transients):
    # regression test: subplots_adjust's top/bottom are fractions of the
    # *total* figure height, which grows with row count here -- a flat
    # fraction (e.g. top=0.90) would waste a huge, growing gap above the
    # first row on a busy cell. The margin should stay ~constant in pixels.
    win_few = ClassifyTransientsWindow(se, 'default')  # cell 0 has 2 transients -> 1 row
    win_many = ClassifyTransientsWindow(se_many_transients, 'default')  # 30 transients -> 6 rows
    try:
        def top_margin_px(win):
            fig_height_px = win.thumb_canvas.height()
            top_frac = win.thumb_fig.subplotpars.top
            return (1 - top_frac) * fig_height_px

        assert top_margin_px(win_few) == pytest.approx(15, abs=0.5)
        assert top_margin_px(win_many) == pytest.approx(15, abs=0.5)
    finally:
        win_few.close()
        win_many.close()
