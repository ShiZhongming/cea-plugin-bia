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


# Generate .csv files to be used as the input for bia visualisation
# for each crop type
# and write to disk
def calc_bia_visual_each():

    """
     This function calculates the number of growth cycles,
     including both initial and subsequent ones, for the selected crop, based on the length of each season.

    :param cycl_i_day: the initial cycle of the selected crop in days
    :type cycl_i_day: integer
    :param cycl_s_day: the subsequent cycle of the selected crop in days
    :type cycl_s_day: integer
    :param n_cycl: the limit of number of growth cycles of the selected crop in days
    :type n_cycl: integer
    :param season_srf:  lengths of season(s) of the selected crop in days for each building surface of a whole year
    :type season_srf: list

    :return cycl_srf: number of cycles, including both initial and subsequent ones,
    for each building surface of a whole year
    :type cycl_srf: list
    :return cycl_i_srf: number of cycles, including initial ones only,
    for each building surface of a whole year
    :type cycl_i_srf: list
    :return cycl_s_srf: number of cycles, including subsequent ones only,
    for each building surface of a whole year
    :type cycl_s_srf: list
    """

# Generate .csv files to be used as the input for bia visualisation
# for all crop types combined
# and write to disk
def calc_bia_visual_all():

    """
     This function calculates the number of growth cycles,
     including both initial and subsequent ones, for the selected crop, based on the length of each season.

    :param cycl_i_day: the initial cycle of the selected crop in days
    :type cycl_i_day: integer
    :param cycl_s_day: the subsequent cycle of the selected crop in days
    :type cycl_s_day: integer
    :param n_cycl: the limit of number of growth cycles of the selected crop in days
    :type n_cycl: integer
    :param season_srf:  lengths of season(s) of the selected crop in days for each building surface of a whole year
    :type season_srf: list

    :return cycl_srf: number of cycles, including both initial and subsequent ones,
    for each building surface of a whole year
    :type cycl_srf: list
    :return cycl_i_srf: number of cycles, including initial ones only,
    for each building surface of a whole year
    :type cycl_i_srf: list
    :return cycl_s_srf: number of cycles, including subsequent ones only,
    for each building surface of a whole year
    :type cycl_s_srf: list
    """