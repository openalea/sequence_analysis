"""Sequence Analysis init file"""
__revision__ = "$Id$"

from openalea.stat_tool import *
#import openalea.stat_tool._stat_tool

from importlib.resources import as_file, files
from os.path import join as pj
import os
from pathlib import Path
import sys

import openalea.stat_tool.interface as interface


from .correlation import *
from .simulate import *
from .compare import *

from .time_events import *
# from top_parameters import *
# from tops import *
from .sequences import *
from .hidden_semi_markov import *
from .hidden_variable_order_markov import *
from .semi_markov import *
from .data_transform  import *

from .estimate import *
from .nonhomogeneous_markov import *
from .renewal import *
from .variable_order_markov import *
from .distance_matrix import *
from .enums_seq import *

try:
    __version__ = version("openalea.stat_tool")
except PackageNotFoundError:
    # package is not installed
    pass

#if sys.platform.startswith("win"):
#    os.add_dll_directory(str(Path(__file__).parent.parent / "lib"))


def get_shared_data(file):
    import openalea.sequence_analysis

    datadir = files("openalea.sequence_analysis.data")
    with as_file(datadir / file) as f:
        return str(f)