"""Classification values and display colors, port of the small static
helpers scattered through @seudo/seudo.m: getClassValue/getClassName,
valTrue/valFalse/valMix/valUnc, pickColor -- plus the classification-cycle
rule from seudoClassifyTransients.m's updateClassification (lines ~990-1000).
"""

import numpy as np

VAL_TRUE = 1.0
VAL_FALSE = -1.0
VAL_MIX = 0.0
VAL_UNC = float('nan')

_COLORS_255 = {
    'true': (0, 8, 228),
    'false': (255, 4, 20),
    'mixed': (80, 200, 80),
    'lsq': (248, 156, 0),
    'artifact': (129, 81, 0),
    'unclassified': (128, 128, 128),
}


def pick_color(which):
    """Returns an (r, g, b) tuple in [0, 1], matching seudo.pickColor."""
    r, g, b = _COLORS_255[which]
    return (r / 255, g / 255, b / 255)


def get_class_name(value):
    if value is None or np.isnan(value):
        return 'unclassified'
    if value == VAL_FALSE:
        return 'false'
    if value == VAL_TRUE:
        return 'true'
    return 'mixed'


def classification_color(value):
    return pick_color(get_class_name(value))


def cycle_classification(value):
    """unclassified -> false -> true -> mixed -> unclassified."""
    if value is None or np.isnan(value):
        return VAL_FALSE
    if value == VAL_FALSE:
        return VAL_TRUE
    if value == VAL_TRUE:
        return VAL_MIX
    return VAL_UNC
