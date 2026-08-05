#!/usr/bin/env bash
# Builds the optional native (C++) SEUDO accelerator extension.
# Not required for the seudo package to work -- estimate_time_courses_with_seudo
# falls back to the pure-Python solver if this hasn't been built.
#
# Requires: a C++14 compiler (g++/clang++) and `pip install pybind11`.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

CXX="${CXX:-g++}"
if ! command -v "$CXX" >/dev/null 2>&1; then
    echo "error: no C++ compiler found ($CXX). Install one (e.g. 'sudo apt-get install build-essential') and retry." >&2
    exit 1
fi

if ! python3 -c "import pybind11" >/dev/null 2>&1; then
    echo "error: pybind11 not installed. Run: pip install pybind11" >&2
    exit 1
fi

PY_INCLUDES=$(python3 -m pybind11 --includes)
PY_EXT_SUFFIX=$(python3-config --extension-suffix)

SOURCES="vendor/fista.cpp vendor/fista_types.cpp vendor/fista_v1.cpp vendor/fista_v1_minimizer.cpp \
         vendor/fista_gradient.cpp vendor/par_fista.cpp vendor/draw.cpp vendor/seudo.cpp \
         vendor/strprintf.cpp vendor/ptwrap.cpp seudo_native_binding.cpp"

echo "building _seudo_native${PY_EXT_SUFFIX} ..."
"$CXX" -std=c++14 -O3 -fPIC -pthread -shared $PY_INCLUDES \
    -I. $SOURCES -o "_seudo_native${PY_EXT_SUFFIX}"
echo "done: $(pwd)/_seudo_native${PY_EXT_SUFFIX}"
