"""Launch the three-panel streaming-pipeline debug viewer on a captured
run. Run scripts/debug_capture_run.py first if scripts/output/
debug_capture.h5 doesn't exist yet -- see seudo/gui/debug_streaming.py's
module docstring for what's captured and why.
"""

import sys

from PyQt5 import QtWidgets

sys.path.insert(0, '.')
from seudo.gui import DebugStreamingWindow  # noqa: E402

H5_PATH = 'scripts/output/debug_capture.h5'


def main():
    app = QtWidgets.QApplication(sys.argv)
    win = DebugStreamingWindow(H5_PATH)
    win.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
