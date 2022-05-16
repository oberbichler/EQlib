try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:
    import importlib_metadata

try:
    __version__ = importlib_metadata.version(__name__)
except:
    __version__ = 'dev'

from eqlib.eqlib_ext import Constraint, Node, Objective, Problem
