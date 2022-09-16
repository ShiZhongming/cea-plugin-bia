"""
This script calculates:
the main script to activate all equations for preparing the data for visualisation.
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

from bia.equation.bia_visual import calc_bia_visual_each
from bia.equation.bia_visual import calc_bia_visual_all



__author__ = "Zhongming Shi"
__copyright__ = "Copyright 2022, Future Cities Laboratory, Singapore - ETH Zurich; " \
                "University of Calgary, Alberta, Canada"
__credits__ = ["Zhongming Shi"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Zhongming Shi"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"


class BiaPlotterPlugin(cea.plugin.CeaPlugin):

    pass



def main(config):
    assert os.path.exists(config.scenario), 'Scenario not found: %s' % config.scenario
    locator = cea.inputlocator.InputLocator(config.scenario, config.plugins)

    # List of crop types included for the building-integrate agriculture (BIA) visualisation
    # At least two types
    types_crop = config.crop_plot.types_crop
    building_names = locator.get_zone_building_names()
    num_process = config.get_number_of_processes()
    n = len(building_names)

    # activate the function that generates
    # .csv files to be used as the input for bia visualisation
    # all the results are stored in the folder "agriculture\plots\"

    for type_crop in range(len(types_crop)):
        # activate the equations for generating .csv files to be used as the input for bia visualisation
        # for each building
        # for each crop type
        # and write to disk
        cea.utilities.parallel.vectorize(calc_bia_visual_each, num_process)\
            (repeat(locator, n), repeat(config, n), building_names, repeat(types_crop[type_crop], n))

    # Generate .csv files to be used as the input for bia visualisation
    # for each building
    # for all crop types combined
    # and write to disk
    cea.utilities.parallel.vectorize(calc_bia_visual_all, num_process)\
        (repeat(locator, n), repeat(config, n), building_names)

if __name__ == '__main__':
    main(cea.config.Configuration())
