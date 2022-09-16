"""
This script creates:
.csv files as the inputs to be visualised in calendar heatmap graphs and dashboards.
.
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
from csv import writer

import pandas as pd

import cea.utilities.parallel
from cea.constants import HOURS_IN_YEAR
from cea.resources.radiation_daysim import daysim_main, geometry_generator
from bia.equation.bia_dli import calc_DLI
from bia.equation.bia_crop_cycle import calc_crop_cycle
from bia.equation.bia_metric import calc_bia_metric
from bia.equation.bia_select import calc_bia_crop_profile

__author__ = "Zhongming Shi"
__copyright__ = "Copyright 2022, Future Cities Laboratory, Singapore - ETH Zurich; " \
                "University of Calgary, Alberta, Canada"
__credits__ = ["Zhongming Shi"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Zhongming Shi"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"


