# EQlib

[![CI](https://github.com/oberbichler/EQlib/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/oberbichler/EQlib/actions) [![PyPI](https://img.shields.io/pypi/v/eqlib)](https://pypi.org/project/eqlib)

Tools for finding equilibrium.

## Installation

```
pip install eqlib
```

## Development

Building from source requires a C++17 compiler. All other dependencies
(Eigen, fmt, spdlog, robin-map, pybind11) are fetched automatically during
the build.

```
pip install -e .
pip install pytest hyperjet
pytest
```

## Reference

If you use EQlib, please refer to the official GitHub repository:

```
@misc{EQlib,
  author = "Thomas Oberbichler",
  title = "EQlib",
  howpublished = "\url{http://github.com/oberbichler/EQlib}",
}
```
