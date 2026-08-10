import h5py
import numpy as np
import pytest

from seudo.gui import DebugStreamingWindow


@pytest.fixture
def capture_path(tmp_path):
    y, x = 20, 20
    n_frames, n_xtemp, n_xstab = 4, 3, 2

    path = tmp_path / 'debug_capture.h5'
    with h5py.File(path, 'w') as f:
        f.create_dataset('avg_frames', data=np.random.default_rng(0).normal(size=(n_frames, y, x)))
        f.create_dataset('avg_frame_indices', data=np.arange(10, 10 + n_frames, dtype=np.int64))

        xtemp_profiles = np.zeros((n_xtemp, y, x))
        for i in range(n_xtemp):
            xtemp_profiles[i, 5 + i, 5 + i] = 1.0
        f.create_dataset('xtemp_profiles', data=xtemp_profiles)
        f.create_dataset('xtemp_track_ids', data=np.arange(n_xtemp, dtype=np.int64))
        f.create_dataset('xtemp_creation_frames', data=np.array([1, 4, 9], dtype=np.int64))
        f.create_dataset('xtemp_last_seen_frames', data=np.array([3, 8, 12], dtype=np.int64))
        f.create_dataset('xtemp_match_counts', data=np.array([2, 3, 1], dtype=np.int64))
        # track 0: dropped; track 1: promoted (became cell 0); track 2: merged into cell 1
        f.create_dataset('xtemp_outcome_codes', data=np.array([0, 1, 2], dtype=np.int64))
        f.create_dataset('xtemp_outcome_cell_ids', data=np.array([-1, 0, 1], dtype=np.int64))
        f.create_dataset('xtemp_outcome_frames', data=np.array([-1, 5, 13], dtype=np.int64))

        xstab_profiles = np.zeros((y, x, n_xstab))
        xstab_profiles[3, 3, 0] = 1.0
        xstab_profiles[10, 10, 1] = 1.0
        f.create_dataset('xstab_profiles', data=xstab_profiles)
        f.create_dataset('xstab_promotion_frames', data=np.array([2, 15], dtype=np.int64))

    return str(path)


@pytest.fixture
def empty_capture_path(tmp_path):
    y, x = 20, 20
    path = tmp_path / 'debug_capture_empty.h5'
    with h5py.File(path, 'w') as f:
        f.create_dataset('avg_frames', data=np.zeros((0, y, x)))
        f.create_dataset('avg_frame_indices', data=np.zeros((0,), dtype=np.int64))
        f.create_dataset('xtemp_profiles', data=np.zeros((0, y, x)))
        f.create_dataset('xtemp_track_ids', data=np.zeros((0,), dtype=np.int64))
        f.create_dataset('xtemp_creation_frames', data=np.zeros((0,), dtype=np.int64))
        f.create_dataset('xtemp_last_seen_frames', data=np.zeros((0,), dtype=np.int64))
        f.create_dataset('xstab_profiles', data=np.zeros((y, x, 0)))
        f.create_dataset('xstab_promotion_frames', data=np.zeros((0,), dtype=np.int64))
    return str(path)


def test_window_builds_and_shows_first_item_of_each_panel(qapp, capture_path):
    win = DebugStreamingWindow(capture_path)
    try:
        assert win._frame_slider.maximum() == 3
        assert win._xtemp_slider.maximum() == 2
        assert win._xstab_slider.maximum() == 1
        assert 'frame 1/4' in win._frame_label.text()
        assert 'reports frame_index 10' in win._frame_label.text()
        assert 'added at frame 1' in win._xtemp_label.text()
        assert 'promoted at frame 2' in win._xstab_label.text()
    finally:
        win.close()


def test_frame_slider_updates_display(qapp, capture_path):
    win = DebugStreamingWindow(capture_path)
    try:
        win._frame_slider.setValue(2)
        assert 'frame 3/4' in win._frame_label.text()
        assert 'reports frame_index 12' in win._frame_label.text()
    finally:
        win.close()


def test_xtemp_slider_updates_display(qapp, capture_path):
    win = DebugStreamingWindow(capture_path)
    try:
        win._xtemp_slider.setValue(2)
        assert 'track 3/3' in win._xtemp_label.text()
        assert 'added at frame 9' in win._xtemp_label.text()
        assert 'last matched at frame 12' in win._xtemp_label.text()
        assert '1 match(es)' in win._xtemp_label.text()
        assert 'merged -> cell 1 at frame 13' in win._xtemp_label.text()
    finally:
        win.close()


def test_xtemp_panel_shows_each_outcome_kind(qapp, capture_path):
    win = DebugStreamingWindow(capture_path)
    try:
        win._xtemp_slider.setValue(0)
        assert 'dropped' in win._xtemp_label.text()
        assert '->' not in win._xtemp_label.text()  # dropped tracks have no destination cell

        win._xtemp_slider.setValue(1)
        assert 'promoted -> cell 0 at frame 5' in win._xtemp_label.text()
    finally:
        win.close()


def test_xtemp_panel_falls_back_gracefully_without_match_counts(qapp, tmp_path):
    # older capture files (before xtemp_match_counts existed) should still
    # load, with an explicit "unavailable" label rather than a crash or a
    # silently wrong number
    y, x = 20, 20
    path = tmp_path / 'debug_capture_no_match_counts.h5'
    with h5py.File(path, 'w') as f:
        f.create_dataset('avg_frames', data=np.zeros((0, y, x)))
        f.create_dataset('avg_frame_indices', data=np.zeros((0,), dtype=np.int64))
        xtemp_profiles = np.zeros((1, y, x))
        xtemp_profiles[0, 5, 5] = 1.0
        f.create_dataset('xtemp_profiles', data=xtemp_profiles)
        f.create_dataset('xtemp_track_ids', data=np.array([0], dtype=np.int64))
        f.create_dataset('xtemp_creation_frames', data=np.array([1], dtype=np.int64))
        f.create_dataset('xtemp_last_seen_frames', data=np.array([3], dtype=np.int64))
        # xtemp_match_counts deliberately omitted
        f.create_dataset('xstab_profiles', data=np.zeros((y, x, 0)))
        f.create_dataset('xstab_promotion_frames', data=np.zeros((0,), dtype=np.int64))

    win = DebugStreamingWindow(str(path))
    try:
        assert 'match count unavailable' in win._xtemp_label.text()
    finally:
        win.close()


def test_xstab_slider_updates_display(qapp, capture_path):
    win = DebugStreamingWindow(capture_path)
    try:
        win._xstab_slider.setValue(1)
        assert 'cell 2/2' in win._xstab_label.text()
        assert 'promoted at frame 15' in win._xstab_label.text()
    finally:
        win.close()


def test_handles_empty_captures_gracefully(qapp, empty_capture_path):
    win = DebugStreamingWindow(empty_capture_path)
    try:
        assert win._frame_slider.maximum() == 0
        assert win._frame_slider.isEnabled() is False
        assert 'no frames captured' in win._frame_label.text()
        assert 'no Xtemp tracks captured' in win._xtemp_label.text()
        assert 'no Xstab cells captured' in win._xstab_label.text()
    finally:
        win.close()
