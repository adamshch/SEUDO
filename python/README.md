# seudo

Python port (in progress) of [SEUDO](https://github.com/adamshch/SEUDO), a
toolbox for removing false transients from calcium imaging source
extraction. Original MATLAB implementation and paper: Gauthier & Charles,
*eLife* 2021.

This package is being ported incrementally from the MATLAB codebase
(`../code` in the source repository). Currently implemented:

- The core SEUDO solver (`estimate_time_courses_with_seudo`) — sparse
  regression separating each source's true activity from unmodeled
  ("blob") contamination, for a single frame at a time.
- Transient detection and per-transient spatial profile computation
  (`identify_transients`, `compute_transient_info`).
- Automatic true/false transient classification by profile correlation
  (`auto_classify_transients`).
- A PyQt5 GUI (`seudo.gui.ClassifyTransientsWindow`) for manually
  reviewing and classifying transients.

Not yet ported: the `rwl1df`/`bpdndf` dynamic SEUDO modes, the `'lbq'`
auto-classification method, parameter search / Bayesian optimization, and
several of the original MATLAB GUIs (parameter review, transient detail
viewer).

## Install

```bash
pip install seudo
```

This includes the PyQt5 classification GUI and h5py (for loading real
MATLAB v7.3 `.mat` files) as standard dependencies, since the GUI is a
core part of the workflow, not an add-on.

For local development, from this directory:

```bash
pip install -e .[dev]
pytest
```

## Quick start

```python
from seudo import Seudo

# movie: (Y, X, F) array; profiles: (Y, X, nCells) array
se = Seudo(movie, profiles, time_courses=time_courses)

se.compute_transient_info('default')
result = se.estimate_time_courses_with_seudo(sigma2=0.002, lambda_blob=10, blob_radius=3)
```

To open the manual classification GUI:

```python
from PyQt5 import QtWidgets
from seudo.gui import ClassifyTransientsWindow

app = QtWidgets.QApplication([])
win = ClassifyTransientsWindow(se, 'default')
win.show()
app.exec_()
```

## License

MIT, see `LICENSE`.
