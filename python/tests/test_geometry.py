import numpy as np

from seudo.geometry import compute_roi_coms


def test_com_and_bounds_single_pixel():
    roi = np.zeros((10, 10))
    roi[3, 7] = 5.0
    coms, bounds = compute_roi_coms(roi)
    assert np.allclose(coms[0], [7, 3])
    assert np.array_equal(bounds[0], [3, 3, 7, 7])


def test_com_symmetric_square():
    roi = np.zeros((10, 10))
    roi[4:7, 4:7] = 1.0
    coms, bounds = compute_roi_coms(roi)
    assert np.allclose(coms[0], [5, 5])
    assert np.array_equal(bounds[0], [4, 6, 4, 6])


def test_all_zero_roi():
    roi = np.zeros((5, 5))
    coms, bounds = compute_roi_coms(roi)
    assert np.isnan(coms[0]).all()
    assert np.array_equal(bounds[0], [0, 0, 0, 0])


def test_stack_of_rois():
    stack = np.zeros((10, 10, 2))
    stack[2, 2, 0] = 1.0
    stack[8, 8, 1] = 1.0
    coms, bounds = compute_roi_coms(stack)
    assert coms.shape == (2, 2)
    assert np.allclose(coms[0], [2, 2])
    assert np.allclose(coms[1], [8, 8])
