"""Builds the optional native (C++) SEUDO accelerator as part of the normal
`pip install` / `python -m build` process.

This is genuinely optional: if no C++ compiler is available (or the build
fails for any other reason), the extension is skipped with a warning and
`pip install` still succeeds -- seudo falls back to the pure-Python solver
at runtime (see seudo/_native/__init__.py). This keeps the package
installable everywhere while building the real accelerator wherever a
toolchain is present, which is what makes `pip install seudo` include
compiled acceleration for typical users without requiring it for everyone.

All the actual project metadata (name, version, dependencies, ...) lives in
pyproject.toml; this file only adds the compiled extension.
"""

from setuptools import setup

try:
    from pybind11.setup_helpers import Pybind11Extension, build_ext as _pybind11_build_ext
    _HAVE_PYBIND11 = True
except ImportError:
    from setuptools.command.build_ext import build_ext as _pybind11_build_ext
    _HAVE_PYBIND11 = False

_VENDOR_DIR = 'seudo/_native/vendor'
_VENDOR_SOURCES = [
    f'{_VENDOR_DIR}/fista.cpp',
    f'{_VENDOR_DIR}/fista_types.cpp',
    f'{_VENDOR_DIR}/fista_v1.cpp',
    f'{_VENDOR_DIR}/fista_v1_minimizer.cpp',
    f'{_VENDOR_DIR}/fista_gradient.cpp',
    f'{_VENDOR_DIR}/par_fista.cpp',
    f'{_VENDOR_DIR}/draw.cpp',
    f'{_VENDOR_DIR}/seudo.cpp',
    f'{_VENDOR_DIR}/strprintf.cpp',
    f'{_VENDOR_DIR}/ptwrap.cpp',
]

ext_modules = []
if _HAVE_PYBIND11:
    ext_modules = [
        Pybind11Extension(
            'seudo._native._seudo_native',
            sources=['seudo/_native/seudo_native_binding.cpp'] + _VENDOR_SOURCES,
            include_dirs=['seudo/_native'],
            cxx_std=14,
            extra_compile_args=['-O3', '-pthread'],
            extra_link_args=['-pthread'],
        ),
    ]


class optional_build_ext(_pybind11_build_ext):
    """Never let a failure to build the native accelerator fail the whole
    install -- catch it, warn, and continue (see module docstring)."""

    def build_extensions(self):
        try:
            super().build_extensions()
        except Exception as exc:  # noqa: BLE001 -- deliberately broad: any build failure should degrade, not abort
            print(f'warning: could not build the native SEUDO accelerator ({exc}).')
            print('seudo will still work, using the pure-Python solver.')

    def build_extension(self, ext):
        try:
            super().build_extension(ext)
        except Exception as exc:  # noqa: BLE001
            print(f'warning: could not build extension {ext.name} ({exc}).')
            print('seudo will still work, using the pure-Python solver.')


setup(
    ext_modules=ext_modules,
    cmdclass={'build_ext': optional_build_ext},
)
