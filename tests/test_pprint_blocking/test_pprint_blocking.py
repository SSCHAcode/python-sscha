"""
Regression tests for issue #196.

When a SSCHA minimization runs in parallel under an MPI launcher
(``mpirun``/``srun``), the standard output is frequently a pipe opened in
*non-blocking* mode. A large write (for instance the table of imaginary
frequencies printed by ``SchaMinimizer.check_imaginary_frequencies``) fills
the pipe buffer and a write on a non-blocking descriptor raises
``BlockingIOError`` ([Errno 11]) instead of waiting for the buffer to drain.
This used to abort the whole calculation.

The fix lives in ``sscha.Parallel.pprint`` (the function aliased as ``print``
across the package): it restores the blocking mode of stdout and, as a last
resort, never lets a log line crash the run.
"""
import os
import sys
import threading

import pytest

import sscha.Parallel

# Non-blocking pipe semantics and os.set_blocking are POSIX-only.
pytestmark = pytest.mark.skipif(
    sys.platform.startswith("win") or not hasattr(os, "set_blocking"),
    reason="requires POSIX non-blocking file descriptors (os.set_blocking)",
)

# Much larger than any pipe buffer (~64 KiB) or stdio buffer, so that the
# write cannot complete in a single non-blocking shot.
BIG_MESSAGE = "x" * (4 * 1024 * 1024)


def _drain(read_fd, sink=None):
    """Consume a pipe until EOF, optionally collecting the bytes."""
    while True:
        chunk = os.read(read_fd, 1 << 16)
        if not chunk:
            break
        if sink is not None:
            sink.extend(chunk)


def test_builtin_print_raises_on_nonblocking_stdout():
    """Reproduce the original failure: a plain ``print`` on a non-blocking
    stdout raises ``BlockingIOError`` once the buffer fills up. This is exactly
    what ``pprint`` used to do before the fix."""
    read_fd, write_fd = os.pipe()
    os.set_blocking(write_fd, False)
    stdout = os.fdopen(write_fd, "w")
    saved = sys.stdout
    sys.stdout = stdout
    try:
        with pytest.raises(BlockingIOError):
            print(BIG_MESSAGE)  # builtin print, no reader draining the pipe
            stdout.flush()
    finally:
        sys.stdout = saved
        # Drain the leftover buffered bytes so closing does not block/raise.
        os.set_blocking(write_fd, True)
        drainer = threading.Thread(target=_drain, args=(read_fd,))
        drainer.start()
        try:
            stdout.close()
        except OSError:
            pass
        drainer.join()
        os.close(read_fd)


def test_pprint_survives_nonblocking_stdout():
    """With the fix, ``pprint`` restores blocking mode and the large write
    completes successfully instead of raising."""
    read_fd, write_fd = os.pipe()
    os.set_blocking(write_fd, False)

    received = bytearray()
    reader = threading.Thread(target=_drain, args=(read_fd, received))
    reader.start()

    stdout = os.fdopen(write_fd, "w")
    saved = sys.stdout
    sys.stdout = stdout
    try:
        # Must not raise BlockingIOError.
        sscha.Parallel.pprint(BIG_MESSAGE)
        stdout.flush()
    finally:
        sys.stdout = saved
        stdout.close()
        reader.join()
        os.close(read_fd)

    assert BIG_MESSAGE.encode() in bytes(received)


def test_pprint_never_raises_when_blocking_cannot_be_set(monkeypatch):
    """Safety net: even if blocking mode cannot be enforced, ``pprint`` must
    swallow the error rather than abort the calculation."""
    monkeypatch.setattr(sscha.Parallel, "_force_stdout_blocking", lambda: None)

    read_fd, write_fd = os.pipe()
    os.set_blocking(write_fd, False)
    stdout = os.fdopen(write_fd, "w")
    saved = sys.stdout
    sys.stdout = stdout
    try:
        # stdout stays non-blocking and nobody reads: the internal write
        # raises BlockingIOError, which pprint must catch.
        sscha.Parallel.pprint(BIG_MESSAGE)
    finally:
        sys.stdout = saved
        os.set_blocking(write_fd, True)
        drainer = threading.Thread(target=_drain, args=(read_fd,))
        drainer.start()
        try:
            stdout.close()
        except OSError:
            pass
        drainer.join()
        os.close(read_fd)
