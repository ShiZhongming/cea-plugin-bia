"""
This script calculates:
the Daily Light Integral (DLI) in [mol/m2/day] for each building envelope surface.
"""

from __future__ import division
from __future__ import print_function

import cea.config
import cea.inputlocator
import cea.plugin

import os
import time
from itertools import repeat
from math import *
from multiprocessing import Pool

import pandas as pd

import cea.utilities.parallel
from cea.constants import HOURS_IN_YEAR
from cea.resources.radiation_daysim import daysim_main, geometry_generator


__author__ = "Zhongming Shi"
__copyright__ = "Copyright 2022, Future Cities Laboratory, Singapore - ETH Zurich"
__credits__ = ["Zhongming Shi"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Zhongming Shi"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"


class BiaCropPlugin(cea.plugin.CeaPlugin):
    """
    Define the plugin class - unless you want to customize the behavior, you only really need to declare the class. The
    rest of the information will be picked up from ``default.config``, ``schemas.yml`` and ``scripts.yml`` by default.
    """
    pass