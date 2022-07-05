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
from bia.equation.bia_crop_cycle import calc_crop_cycle
from bia.equation.bia_metric import calc_bia_metric, bia_result_aggregate_write

__author__ = "Zhongming Shi"
__copyright__ = "Copyright 2022, Future Cities Laboratory, Singapore - ETH Zurich; " \
                "University of Calgary, Alberta, Canada"
__credits__ = ["Zhongming Shi"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Zhongming Shi"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"


class BiaAssessmentPlugin(cea.plugin.CeaPlugin):

    pass


def main(config):
    assert os.path.exists(config.scenario), 'Scenario not found: %s' % config.scenario
    locator = cea.inputlocator.InputLocator(config.scenario, config.plugins)

    # BIA assessment
    print('Executing CEA Building-integrated agriculture assessment (BIA) for {type_crop}'
          .format(type_crop=config.agriculture.type_crop))
    building_names = locator.get_zone_building_names()
    num_process = config.get_number_of_processes()
    n = len(building_names)
    cea.utilities.parallel.vectorize(calc_crop_cycle, num_process)(repeat(config, n), building_names)
    cea.utilities.parallel.vectorize(calc_bia_metric, num_process)(repeat(config, n), building_names)

    # aggregate the results of each building as a 'total' file and write to disk
    bia_result_aggregate_write(locator, config, building_name)

if __name__ == '__main__':
    main(cea.config.Configuration())
