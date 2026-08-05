import numpy as np

from seudo.constants import (
    VAL_FALSE, VAL_MIX, VAL_TRUE, VAL_UNC,
    classification_color, cycle_classification, get_class_name, pick_color,
)


def test_cycle_classification_order():
    v = VAL_UNC
    v = cycle_classification(v)
    assert v == VAL_FALSE
    v = cycle_classification(v)
    assert v == VAL_TRUE
    v = cycle_classification(v)
    assert v == VAL_MIX
    v = cycle_classification(v)
    assert np.isnan(v)


def test_get_class_name():
    assert get_class_name(VAL_TRUE) == 'true'
    assert get_class_name(VAL_FALSE) == 'false'
    assert get_class_name(VAL_MIX) == 'mixed'
    assert get_class_name(VAL_UNC) == 'unclassified'
    assert get_class_name(float('nan')) == 'unclassified'


def test_colors_are_normalized_rgb():
    for name in ('true', 'false', 'mixed', 'artifact', 'unclassified'):
        r, g, b = pick_color(name)
        assert 0 <= r <= 1 and 0 <= g <= 1 and 0 <= b <= 1

    assert classification_color(VAL_TRUE) == pick_color('true')
    assert classification_color(float('nan')) == pick_color('unclassified')
