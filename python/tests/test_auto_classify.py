import numpy as np
import pytest

from seudo.auto_classify import auto_classify_transients
from seudo.constants import VAL_FALSE, VAL_TRUE
from seudo.core import Seudo


def gaussian_patch(shape, center, sigma):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    g = np.exp(-(((yy - center[0]) ** 2 + (xx - center[1]) ** 2) / (2 * sigma ** 2)))
    g[g < 0.05 * g.max()] = 0.0
    return g


def make_seudo_with_real_and_contaminated_transient():
    """Cell 0 has two detected transients: a real one (its own profile is
    truly active) and a contaminated one (cell 1's profile is what's
    actually active in the movie, but cell 0's own -- externally supplied,
    as if from a naive extraction -- time course shows a matching bump
    anyway, so compute_transient_info still detects and windows it as one
    of cell 0's transients). auto_classify_transients should tell them
    apart: high correlation / low residual ratio / True for the real one,
    the opposite for the contaminated one.
    """
    y, x, f = 20, 40, 600
    prof0 = gaussian_patch((y, x), center=(10, 15), sigma=2.0)
    prof1 = gaussian_patch((y, x), center=(10, 25), sigma=2.0)
    prof2 = np.zeros((y, x))  # silent cell, no transients at all

    true_act0 = np.zeros(f)
    true_act0[100:110] = 3.0
    act1 = np.zeros(f)
    act1[300:310] = 3.0

    rng = np.random.default_rng(4)
    movie = np.zeros((y, x, f))
    for ff in range(f):
        movie[:, :, ff] = prof0 * true_act0[ff] + prof1 * act1[ff] + 0.01 * rng.normal(size=(y, x))

    # externally-supplied tc for cell 0: correct at the real transient, but
    # also shows a matching bump during cell 1's activity (contamination)
    tc0 = true_act0.copy()
    tc0[300:310] = 3.0
    tc1 = act1.copy()
    tc2 = np.zeros(f)

    profiles = np.stack([prof0, prof1, prof2], axis=2)
    tc = np.stack([tc0, tc1, tc2], axis=1)

    se = Seudo(movie, profiles, time_courses=tc)
    se.compute_transient_info('default', min_duration=3, blur_radius=1)
    return se


@pytest.fixture
def se():
    return make_seudo_with_real_and_contaminated_transient()


def test_distinguishes_real_from_contaminated_transient(se):
    results = auto_classify_transients(se, 'default')

    ti0 = se.tc_default['transient_info'][0]
    assert ti0['times'].shape[0] == 2

    real_idx, contam_idx = (0, 1) if ti0['times'][0, 0] < ti0['times'][1, 0] else (1, 0)

    corrs = results[0]['corrs']
    res_ratios = results[0]['resRatios']
    classification = results[0]['classification']

    assert corrs[real_idx] > 0.9
    assert corrs[contam_idx] < corrs[real_idx]
    assert res_ratios[real_idx] < res_ratios[contam_idx]
    assert classification[real_idx] == VAL_TRUE
    assert classification[contam_idx] == VAL_FALSE

    # results were written back into the seudo object
    assert np.array_equal(ti0['classification'], classification)
    assert ti0['auto_class'] is not None


def test_cell_with_no_transients_returns_empty_metric(se):
    results = auto_classify_transients(se, 'default', which_cells=[2])
    assert np.isnan(results[0]['cfrac'])
    assert results[0]['corrs'].size == 0
    assert results[0]['classification'].size == 0
    assert se.tc_default['transient_info'][2]['classification'].size == 0


def test_progress_callback_invoked_per_cell(se):
    calls = []
    auto_classify_transients(
        se, 'default', which_cells=[0, 1, 2],
        progress_callback=lambda done, total, cell_id: calls.append((done, total, cell_id)))
    assert calls == [(1, 3, 0), (2, 3, 1), (3, 3, 2)]


def test_which_cells_restricts_processing(se):
    results = auto_classify_transients(se, 'default', which_cells=[1])
    assert len(results) == 1
    # cell 0 untouched -- still all-nan
    assert np.all(np.isnan(se.tc_default['transient_info'][0]['classification']))


def test_overwrite_protection(se):
    auto_classify_transients(se, 'default')  # first run, all cells unclassified -> fine

    with pytest.raises(ValueError, match='overwrite'):
        auto_classify_transients(se, 'default')

    # doesn't raise with overwrite=True
    auto_classify_transients(se, 'default', overwrite=True)


def test_corr_thresh_changes_classification(se):
    results_lenient = auto_classify_transients(se, 'default', which_cells=[0], save_results=False, corr_thresh=-1.0)
    assert np.all(results_lenient[0]['classification'] == VAL_TRUE)

    results_strict = auto_classify_transients(se, 'default', which_cells=[0], save_results=False, corr_thresh=1.1)
    assert np.all(results_strict[0]['classification'] == VAL_FALSE)


def test_res_ratio_criterion_changes_classification(se):
    results_lenient = auto_classify_transients(
        se, 'default', which_cells=[0], save_results=False, criterion='res_ratio', res_ratio_thresh=1e6)
    assert np.all(results_lenient[0]['classification'] == VAL_TRUE)

    results_strict = auto_classify_transients(
        se, 'default', which_cells=[0], save_results=False, criterion='res_ratio', res_ratio_thresh=-1.0)
    assert np.all(results_strict[0]['classification'] == VAL_FALSE)


def test_invalid_criterion_raises(se):
    with pytest.raises(ValueError, match='criterion'):
        auto_classify_transients(se, 'default', which_cells=[0], save_results=False, criterion='bogus')


def test_seudo_residual_distinguishes_real_from_contaminated_transient(se):
    # lambda_blob is scaled to the real demo dataset's pixel intensities by
    # default; this synthetic fixture's shapes peak around 1.0, so it needs
    # a much smaller value for any signal to survive the sparsity penalty.
    results = auto_classify_transients(
        se, 'default', which_cells=[0], save_results=False, criterion='seudo_residual',
        seudo_kwargs=dict(lambda_blob=0.1))

    ti0 = se.tc_default['transient_info'][0]
    real_idx, contam_idx = (0, 1) if ti0['times'][0, 0] < ti0['times'][1, 0] else (1, 0)

    fractions = results[0]['seudoResidual']
    assert fractions.shape[0] == 2
    # real transient: SEUDO's own-cell coefficient should closely match the
    # plain least-squares one (small fraction). Contaminated one: SEUDO
    # diverts most of the amplitude into the blob basis instead (fraction
    # much closer to 1).
    assert fractions[real_idx] < fractions[contam_idx]
    assert fractions[real_idx] < 0.5
    assert fractions[contam_idx] > 0.5


def test_seudo_residual_thresh_changes_classification(se):
    # high fraction is the "true" direction for this criterion (per
    # explicit user direction, opposite of res_ratio) -- so a very low
    # threshold is lenient (catches everything) and a very high one is strict.
    results_lenient = auto_classify_transients(
        se, 'default', which_cells=[0], save_results=False,
        criterion='seudo_residual', seudo_residual_thresh=-1e6)
    assert np.all(results_lenient[0]['classification'] == VAL_TRUE)

    results_strict = auto_classify_transients(
        se, 'default', which_cells=[0], save_results=False,
        criterion='seudo_residual', seudo_residual_thresh=1e6)
    assert np.all(results_strict[0]['classification'] == VAL_FALSE)
