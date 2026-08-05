import numpy as np

from seudo.classification_io import load_classification, save_classification


def test_save_load_round_trip(tmp_path):
    tc = np.random.default_rng(0).normal(size=(50, 2))
    transient_info = [
        {'times': np.array([[5, 10]]), 'classification': np.array([1.0]), 'is_artifact': False},
        {'times': np.zeros((0, 2), dtype=int), 'classification': np.array([]), 'is_artifact': True},
    ]
    tc_struct = {'tc': tc, 'transient_info': transient_info}

    path = tmp_path / 'classification.pkl'
    save_classification(tc_struct, path)
    loaded = load_classification(path)

    assert np.array_equal(loaded['tc'], tc)
    assert np.array_equal(loaded['transient_info'][0]['times'], [[5, 10]])
    assert loaded['transient_info'][0]['classification'][0] == 1.0
    assert loaded['transient_info'][1]['is_artifact'] is True
