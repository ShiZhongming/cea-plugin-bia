"""
This script calculates:
the main script to activate all BIA-related assessments.
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

    type_crop = config.agriculture.type_crop
    print('Executing CEA Building-Integrated Agriculture (BIA) assessment for {type_crop}.'
          .format(type_crop=type_crop))
    building_names = locator.get_zone_building_names()
    num_process = config.get_number_of_processes()
    n = len(building_names)

    # DLI calculations for each surface of each building
    cea.utilities.parallel.vectorize(calc_DLI, num_process)\
        (repeat(locator, n), repeat(config, n), building_names)

    # BIA metrics for each surface of each building
    cea.utilities.parallel.vectorize(calc_bia_metric, num_process)\
        (repeat(locator, n), repeat(config, n), building_names, repeat(type_crop, n))

    # aggregate the results of each building as a 'total' file and write to disk
    # if the file exists, delete it
    bia_path = config.scenario + "/outputs/data/potentials/agriculture/BIA_assessment_total_{type_crop}.csv" \
        .format(type_crop=type_crop)
    if os.path.exists(bia_path):
        os.remove(bia_path)
    cea.utilities.parallel.vectorize(bia_result_aggregate_write, num_process)\
        (repeat(locator, n), repeat(config, n), building_names, repeat(type_crop, n))

if __name__ == '__main__':
    main(cea.config.Configuration())
