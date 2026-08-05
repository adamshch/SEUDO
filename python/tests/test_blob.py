import numpy as np

from seudo.blob import make_seudo_blob, matlab_fspecial_gaussian


def test_gaussian_kernel_properties():
    h = matlab_fspecial_gaussian(9, 1.5)
    assert h.shape == (9, 9)
    assert np.isclose(h.sum(), 1.0)
    # peak at center
    assert np.argmax(h) == h.size // 2


def test_seudo_blob_normalized_and_peaked():
    blob = make_seudo_blob(1.2)
    assert blob.shape[0] == blob.shape[1]
    assert blob.shape[0] % 2 == 1
    assert np.isclose(np.sum(blob ** 2), 1.0)
    center = blob.shape[0] // 2
    assert blob[center, center] == blob.max()
