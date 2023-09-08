"""
This script calculates:
the main script to activate all crop profiling equations.
"""

from __future__ import division
from __future__ import print_function
import cea.config
import cea.inputlocator
import cea.plugin
import os
from itertools import repeat
import cea.utilities.parallel
from bia.equation.bia_metric import calc_bia_metric
from bia.equation.bia_select import calc_bia_crop_profile

__author__ = "Zhongming Shi"
__copyright__ = "Copyright 2022, Future Cities Laboratory, Singapore - ETH Zurich; " \
                "University of Calgary, Alberta, Canada"
__credits__ = ["Zhongming Shi"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Zhongming Shi"
__email__ = "shi@arch.ethz.ch"
__status__ = "Production"


class BiaProfilerPlugin(cea.plugin.CeaPlugin):

    pass



def main(config):
    assert os.path.exists(config.scenario), 'Scenario not found: %s' % config.scenario
    locator = cea.inputlocator.InputLocator(config.scenario, config.plugins)

    # List of crop types considered for the building-integrate agriculture (BIA) crop profiling
    # At least two types
    types_crop = config.crop_profile.types_crop
    building_names = locator.get_zone_building_names()
    num_process = config.get_number_of_processes()
    n = len(building_names)

    # activate the function that calculates
    # the BIA metrics for each building surface for each candidate crop type
    # all the results are stored in the folder "agriculture\surface\"

    for type_crop in range(len(types_crop)):
        # activate the bia metric equations for every surface
        cea.utilities.parallel.vectorize(calc_bia_metric, num_process)\
            (repeat(locator, n), repeat(config, n), building_names, repeat(types_crop[type_crop], n))

    # Create the crop profiles for each building's surface and write to disk
    cea.utilities.parallel.vectorize(calc_bia_crop_profile, num_process)\
        (repeat(locator, n), repeat(config, n), building_names)

if __name__ == '__main__':
    main(cea.config.Configuration())
