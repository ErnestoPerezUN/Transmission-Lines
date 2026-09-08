# Academic code for the Transmission Line courses at Universidad Nacional de Colombia.
# No warranty of any kind; not for real-world design. See DISCLAIMER.md.
# License: to be defined (open source intended); until then all rights reserved.
"""Make ``src/fdtd`` importable without packaging the repository.

The rest of this repo is written as standalone scripts rather than as an
installed package, so the tests follow the same convention and simply put the
module directory on ``sys.path``.
"""
import pathlib
import sys

_SRC = pathlib.Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(_SRC / "fdtd"))
sys.path.insert(0, str(_SRC / "line_profile"))
