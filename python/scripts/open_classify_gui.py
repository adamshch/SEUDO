"""Launch the transient classification GUI on the real demo dataset
(code/demoData1.mat), mirroring demo.m's:

    se = seudo(M,P,'timeCourses',T,'name','demo dataset')
    se.computeTransientInfo('default','transientFrames',TF)
    se.classifyTransients

Movie access is lazy (Hdf5Movie) so the ~1.3GB file is never fully loaded
into memory -- only the frames touched by an actual transient are read.
"""

import sys
import time
import traceback

import h5py
from PyQt5 import QtCore, QtWidgets

_orig_excepthook = sys.excepthook


def _excepthook(exc_type, exc_value, exc_tb):
    print('UNHANDLED EXCEPTION IN QT CALLBACK:', flush=True)
    traceback.print_exception(exc_type, exc_value, exc_tb)
    _orig_excepthook(exc_type, exc_value, exc_tb)


sys.excepthook = _excepthook

sys.path.insert(0, '.')
from seudo.core import Seudo  # noqa: E402
from seudo.gui.classify_transients import ClassifyTransientsWindow  # noqa: E402
from seudo.matlab_io import Hdf5Movie, load_profiles, load_time_courses, load_transient_frames  # noqa: E402

DEMO_PATH = '../code/demoData1.mat'


def main():
    print(f'opening {DEMO_PATH} ...')
    f = h5py.File(DEMO_PATH, 'r')

    profiles = load_profiles(f)
    T = load_time_courses(f, 'T')
    TF = load_transient_frames(f, 'TF')
    movie = Hdf5Movie(f['M'])

    se = Seudo(movie, profiles, time_courses=T, name='demo dataset')
    print(f'{se.n_cells} cells, {se.mov_f} frames, {se.mov_y}x{se.mov_x} pixels')

    print('computing transient info for all cells (using the provided TF)...')
    t0 = time.time()
    se.compute_transient_info('default', transient_frames=TF)
    n_trans = sum(ti['times'].shape[0] for ti in se.tc_default['transient_info'])
    print(f'done in {time.time() - t0:0.1f}s -- {n_trans} transients total')

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
    app.setQuitOnLastWindowClosed(True)

    win = ClassifyTransientsWindow(se, 'default')
    win.resize(1400, 900)
    win.show()
    print(f'GUI window opened -- isVisible={win.isVisible()}, entering event loop', flush=True)

    heartbeat = QtCore.QTimer()
    heartbeat.timeout.connect(lambda: print('heartbeat, still alive', flush=True))
    heartbeat.start(3000)

    exit_code = app.exec_()
    print(f'exec_() returned with code {exit_code}', flush=True)

    f.close()


if __name__ == '__main__':
    main()
