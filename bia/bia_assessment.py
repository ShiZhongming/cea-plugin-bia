"""
This script calculates:
the Daily Light Integral (DLI) in [mol/m2/day], the crop yields (kg), carbon emission reduction (CO2-equivalent),
cost (CAPEX and OPEX in USD) for the selected crop type on each building envelope surface.
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
from bia.equation.bia_dli import calc_DLI
from bia.equation.bia_crop import calc_crop


__author__ = "Zhongming Shi"
__copyright__ = "Copyright 2022, Future Cities Laboratory, Singapore - ETH Zurich"
__credits__ = ["Zhongming Shi"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Zhongming Shi"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"


class BiaAssessmentPlugin(cea.plugin.CeaPlugin):
    """
    Define the plugin class - unless you want to customize the behavior, you only really need to declare the class. The
    rest of the information will be picked up from ``default.config``, ``schemas.yml`` and ``scripts.yml`` by default.
    """
    pass



def main(config):
    assert os.path.exists(config.scenario), 'Scenario not found: %s' % config.scenario
    locator = cea.inputlocator.InputLocator(config.scenario, config.plugins)

    # DLI calculations
    print('Running Day Light Integral calculations with scenario = %s' % config.scenario)
    print('Running Day Light Integral calculations with annual-radiation-threshold-kWh/m2 = %s'
          % config.agriculture.annual_radiation_threshold_BIA)
    print('Running Day Light Integral calculations with crop-on-roof = %s' % config.agriculture.crop_on_roof)
    print('Running Day Light Integral calculations with crop-on-wall = %s' % config.agriculture.crop_on_wall)

    building_names = locator.get_zone_building_names()
    num_process = config.get_number_of_processes()
    n = len(building_names)
    cea.utilities.parallel.vectorize(calc_DLI, num_process)(repeat(locator, n), repeat(config, n), building_names)

    # BIA assessment
    # print('Running BIA assessment for the selected crop of {crop_type}'.format(crop_type=config.agriculture.crop_type))
    #
    # building_names = locator.get_zone_building_names()
    # num_process = config.get_number_of_processes()
    # n = len(building_names)
    # cea.utilities.parallel.vectorize(calc_crop, num_process)(repeat(locator, n), repeat(config, n), building_names)


if __name__ == '__main__':
    main(cea.config.Configuration())
