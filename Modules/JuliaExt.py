"""
Lazy bridge to the Julia runtime.

This module is the single entry point for all the Julia calls of python-sscha.
The Julia runtime is NOT booted when this module (or sscha) is imported;
it is initialized on the first call to :func:`get_main`, so that users who do
not use ``fourier gradient`` never pay the startup cost.

Two backends are supported:

* ``juliacall`` (PythonCall.jl) — the default. It has no libpython coupling
  (works with any Python interpreter, no ``python-jl`` needed) and installs
  Julia automatically on a fresh machine through ``juliapkg``.
* ``pyjulia`` (PyCall.jl) — legacy. It is used only if juliacall is not
  installed, or if another package (e.g. an old python-sscha) has already
  booted the PyJulia runtime in this process: only one Julia runtime can
  exist per process, so in that case we must reuse it.

The backend can be forced with the environment variable
``SSCHA_JULIA_BACKEND`` set to ``juliacall``, ``pyjulia`` or ``none``.

The object returned by :func:`get_main` mimics ``julia.Main`` from PyJulia:
attribute access gives callables, ``eval`` evaluates a code string and
``include`` loads a file. Under juliacall, numpy array arguments are
converted to native Julia ``Array``s (PyJulia semantics), because the
sscha Julia kernels use strictly-typed signatures that do not dispatch
on the no-copy ``PyArray`` wrappers, and array results are converted
back to numpy.
"""

import importlib.util
import os
import sys
import threading

import numpy as np

# The Julia source files defining the sscha kernels, included at first use.
_JL_FILES = ["fourier_gradient.jl"]

_BACKEND_ENV = "SSCHA_JULIA_BACKEND"

_lock = threading.RLock()
_main = None
_init_error = None


class JuliaError(ImportError):
    pass


def _requested_backend():
    backend = os.environ.get(_BACKEND_ENV, "").strip().lower()
    if backend in ("juliacall", "pyjulia", "none"):
        return backend
    return ""


def available():
    """Check if a Julia backend is installed, WITHOUT booting the runtime.

    This is a cheap check (importlib.find_spec). The runtime may still fail
    to initialize at first use (e.g. broken installation); in that case
    :func:`get_main` raises a JuliaError with the details.
    """
    if _main is not None:
        return True
    if _init_error is not None:
        return False

    backend = _requested_backend()
    if backend == "none":
        return False
    if backend == "juliacall":
        return importlib.util.find_spec("juliacall") is not None
    if backend == "pyjulia":
        return importlib.util.find_spec("julia") is not None

    if "julia.Main" in sys.modules:
        return True
    return (importlib.util.find_spec("juliacall") is not None
            or importlib.util.find_spec("julia") is not None)


def get_main():
    """Return the (lazily initialized) Julia Main proxy.

    The first call boots the Julia runtime and includes the sscha Julia
    sources; subsequent calls return the cached proxy. Raises JuliaError if
    no working backend is available.
    """
    global _main, _init_error

    if _main is not None:
        return _main

    with _lock:
        if _main is not None:
            return _main
        if _init_error is not None:
            raise JuliaError(
                "The Julia extension failed to initialize earlier:\n{}".format(
                    _init_error))
        try:
            _main = _initialize()
        except Exception as e:
            _init_error = "{}: {}".format(type(e).__name__, e)
            raise JuliaError(
                "Could not initialize the Julia extension.\n"
                "Install it with: pip install juliacall\n"
                "(Julia itself is downloaded automatically at first use.)\n"
                "Original error: {}".format(_init_error))
        return _main


def _initialize():
    backend = _requested_backend()
    if backend == "none":
        raise JuliaError(
            "The Julia extension is disabled ({}=none).".format(_BACKEND_ENV))

    if not backend:
        if "julia.Main" in sys.modules:
            # PyJulia is already running in this process (e.g. booted by
            # python-sscha); a second runtime cannot be created, reuse it.
            backend = "pyjulia"
        elif importlib.util.find_spec("juliacall") is not None:
            backend = "juliacall"
        elif importlib.util.find_spec("julia") is not None:
            backend = "pyjulia"
        else:
            raise JuliaError("Neither juliacall nor pyjulia is installed.")

    if backend == "juliacall":
        main = _init_juliacall()
    else:
        main = _init_pyjulia()

    dirname = os.path.dirname(os.path.abspath(__file__))
    for fname in _JL_FILES:
        main.include(os.path.join(dirname, fname))
    return main


def _init_juliacall():
    # These must be set before the first "import juliacall".
    # PyJulia honored JULIA_NUM_THREADS and defaulted to a single thread:
    # keep exactly that behavior (some kernels use shared buffers inside
    # Threads.@threads loops and are NOT safe with more than one thread).
    n_threads = os.environ.get("JULIA_NUM_THREADS", "1")
    os.environ.setdefault("PYTHON_JULIACALL_THREADS", n_threads)
    if os.environ.get("PYTHON_JULIACALL_THREADS", "1") != "1":
        # Required for safe Julia multithreading from Python.
        os.environ.setdefault("PYTHON_JULIACALL_HANDLE_SIGNALS", "yes")

    from juliacall import Main, convert
    return _JuliaCallMain(Main, convert)


def _init_pyjulia():
    import julia

    if "julia.Main" not in sys.modules:
        try:
            import julia.Main
        except Exception:
            # Statically linked python or libpython mismatch: PyCall cannot
            # use its precompiled cache. This path is slow (it recompiles
            # PyCall at every launch): prefer installing juliacall.
            from julia.api import Julia
            Julia(compiled_modules=False)
            import julia.Main

    return _PyJuliaMain(sys.modules["julia.Main"])


class _PyJuliaMain(object):
    """PyJulia passthrough: julia.Main already speaks numpy."""

    def __init__(self, main):
        self._main = main

    def include(self, path):
        return self._main.include(path)

    def eval(self, code):
        return self._main.eval(code)

    def __getattr__(self, name):
        return getattr(self._main, name)


# numpy dtype -> Julia element type, for the nested Vector{Vector{T}} case
_JL_ELTYPE = {
    "float32": "Float32",
    "float64": "Float64",
    "int32": "Int32",
    "int64": "Int64",
    "complex64": "ComplexF32",
    "complex128": "ComplexF64",
    "bool": "Bool",
}


class _JuliaCallMain(object):
    """juliacall proxy restoring PyJulia argument/return conventions."""

    def __init__(self, main, convert):
        # Avoid __getattr__ recursion: set everything through __dict__.
        self.__dict__["_main"] = main
        self.__dict__["_convert"] = convert
        self.__dict__["_array_type"] = main.seval("Array")
        self.__dict__["_nested_types"] = {}

    def include(self, path):
        return self._main.include(path)

    def eval(self, code):
        return self._from_julia(self._main.seval(code))

    def __getattr__(self, name):
        func = getattr(self._main, name)

        def _call(*args, **kwargs):
            jl_args = [self._to_julia(a) for a in args]
            jl_kwargs = {k: self._to_julia(v) for k, v in kwargs.items()}
            return self._from_julia(func(*jl_args, **jl_kwargs))

        _call.__name__ = name
        return _call

    def _to_julia(self, x):
        if isinstance(x, np.ndarray):
            return self._convert(self._array_type, x)

        # Lists/tuples of 1d arrays with a common dtype are what the
        # sparse-symmetry initializers expect as Vector{Vector{T}}.
        if (isinstance(x, (list, tuple)) and len(x) > 0
                and all(isinstance(e, np.ndarray) and e.ndim == 1 for e in x)):
            eltype = _JL_ELTYPE.get(x[0].dtype.name)
            if eltype is not None and all(e.dtype == x[0].dtype for e in x):
                nested = self._nested_types.get(eltype)
                if nested is None:
                    nested = self._main.seval(
                        "Vector{{Vector{{{}}}}}".format(eltype))
                    self._nested_types[eltype] = nested
                return self._convert(nested, list(x))

        return x

    def _from_julia(self, x):
        if isinstance(x, tuple):
            return tuple(self._from_julia(e) for e in x)

        import juliacall
        if isinstance(x, juliacall.ArrayValue):
            # Buffer-protocol view on the Julia data; numpy keeps the
            # wrapper alive through ndarray.base.
            return np.asarray(x)
        return x
