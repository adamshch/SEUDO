import numpy as np

from seudo.core import Seudo
from seudo.transients import identify_transients


def test_identify_transients_finds_known_spike():
    rng = np.random.default_rng(0)
    f = 200
    tc = 0.05 * rng.normal(size=(f, 1))
    tc[80:100, 0] += 5.0  # clear transient, well above baseline noise

    detected = identify_transients(tc, blur_radius=1, min_duration=3)

    assert detected.shape == (f, 1)
    detected_frames = np.flatnonzero(detected[:, 0])
    assert detected_frames.size > 0
    # detected span should overlap the true spike
    assert detected_frames.min() <= 100
    assert detected_frames.max() >= 80


def test_identify_transients_drops_brief_blips():
    rng = np.random.default_rng(1)
    f = 200
    tc = 0.05 * rng.normal(size=(f, 1))
    tc[50, 0] += 5.0  # single-frame blip, shorter than min_duration

    detected = identify_transients(tc, blur_radius=0, min_duration=6)
    assert not detected.any()


def gaussian_patch(shape, center, sigma, amplitude=1.0):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    return amplitude * np.exp(-(((yy - center[0]) ** 2 + (xx - center[1]) ** 2) / (2 * sigma ** 2)))


def make_synthetic_seudo():
    # frame count / transient duration chosen so elevated frames are a small
    # minority (as in real data) -- the default detection threshold is based
    # on the 95th percentile, which only works as a discriminator when most
    # frames are baseline.
    y, x, f = 20, 20, 600
    profile = gaussian_patch((y, x), center=(10, 10), sigma=2.0)
    profile[profile < 0.05 * profile.max()] = 0.0

    activity = np.zeros(f)
    activity[100:110] = 3.0
    activity[300:310] = 3.0

    rng = np.random.default_rng(2)
    movie = np.zeros((y, x, f))
    for ff in range(f):
        movie[:, :, ff] = profile * activity[ff] + 0.01 * rng.normal(size=(y, x))

    profiles = profile[:, :, np.newaxis]
    tc = activity[:, np.newaxis]

    se = Seudo(movie, profiles, time_courses=tc)
    return se, profile, activity


def test_compute_transient_info_end_to_end():
    se, profile, activity = make_synthetic_seudo()

    info = se.compute_transient_info('default', min_duration=3, blur_radius=1)

    assert len(info) == 1
    cell_info = info[0]

    # two transients should be found, roughly at [100,110) and [300,310)
    assert cell_info['times'].shape[0] == 2
    assert cell_info['classification'].shape == (2,)
    assert np.all(np.isnan(cell_info['classification']))

    # shapes should closely resemble the true profile (high correlation)
    assert cell_info['shapes'] is not None
    assert cell_info['shapes'].shape[2] == 2
    assert np.all(cell_info['corr_with_profile'] > 0.9)


def test_compute_transient_info_preserves_classification_on_rerun():
    se, profile, activity = make_synthetic_seudo()

    info = se.compute_transient_info('default', min_duration=3, blur_radius=1)
    info[0]['classification'][:] = [1.0, -1.0]

    # rerun with identical settings -> times unchanged -> classification preserved
    info2 = se.compute_transient_info('default', min_duration=3, blur_radius=1)
    assert np.array_equal(info2[0]['classification'], [1.0, -1.0])


def test_compute_transient_info_with_explicit_transient_frames():
    se, profile, activity = make_synthetic_seudo()

    transient_frames = np.zeros((se.mov_f, se.n_cells), dtype=bool)
    transient_frames[100:110, 0] = True
    transient_frames[300:310, 0] = True

    info = se.compute_transient_info('default', transient_frames=transient_frames, t_pre=0, t_post=0)
    times = info[0]['times']
    assert np.array_equal(times, [[100, 109], [300, 309]])
