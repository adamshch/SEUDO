"""Optional native (C++) accelerator for the SEUDO per-frame FISTA solve.

Not required for seudo to work: estimate_time_courses_with_seudo falls back
to the pure-Python solver (solver.fista_nonneg_weighted_l1) transparently
if this extension hasn't been built. To build it:

    cd seudo/_native && ./build_native.sh

(requires a C++14 compiler and `pip install pybind11`).

Wraps realSEUDO's Fista::Seudo class (vendor/, MIT licensed, see
vendor/LICENSE) -- the same nonneg-constrained weighted-L1 FISTA solve as
solver.fista_nonneg_weighted_l1, compiled and multithreaded.
"""

try:
    from ._seudo_native import seudo_solve
    NATIVE_AVAILABLE = True
except ImportError:
    seudo_solve = None
    NATIVE_AVAILABLE = False

__all__ = ['seudo_solve', 'NATIVE_AVAILABLE']
