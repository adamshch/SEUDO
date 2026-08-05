"""Save/load a time-course struct's classification state to disk.

Port of the relevant part of seudoClassifyTransients.m's
pushbuttonSaveClassification_Callback -- MATLAB saves the whole tcStruct
(tc + transientInfo) as a .mat file; here we pickle the equivalent Python
dict, since this is a Python-to-Python round trip (no MATLAB interop is
attempted).
"""

import pickle


def save_classification(tc_struct, path):
    payload = {
        'tc': tc_struct['tc'],
        'transient_info': tc_struct.get('transient_info'),
    }
    with open(path, 'wb') as fh:
        pickle.dump(payload, fh)


def load_classification(path):
    with open(path, 'rb') as fh:
        return pickle.load(fh)


def load_classification_into(se, path, which_struct='default'):
    """Load a saved classification into se's chosen tc struct, replacing
    its 'transient_info' (and 'tc', if that struct didn't already have one).
    """
    payload = load_classification(path)
    tc_struct = se._resolve_tc_struct(which_struct)
    if tc_struct.get('tc') is None:
        tc_struct['tc'] = payload['tc']
    tc_struct['transient_info'] = payload['transient_info']
    return tc_struct
