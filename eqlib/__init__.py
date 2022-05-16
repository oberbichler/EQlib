import importlib.metadata as importlib_metadata

__version__ = importlib_metadata.version(__name__)

from eqlib.eqlib_ext import Constraint, Node, Objective, Problem
