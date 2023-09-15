"""
This script calculates:
the main script to activate DLI calculation.
"""

from __future__ import division
from __future__ import print_function
import cea.config
import cea.inputlocator
import cea.plugin
import os
from itertools import repeat
import cea.utilities.parallel
from bia.equation.bia_dli import calc_DLI
from bia.equation.bia_metric import calc_bia_metric, bia_result_aggregate_write

__author__ = "Zhongming Shi"
__copyright__ = "Copyright 2023, A/S Group, ITA, ETH Zurich"
__credits__ = ["Zhongming Shi"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Zhongming Shi"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"


class BiaDliPlugin(cea.plugin.CeaPlugin):

    pass


def main(config):
    assert os.path.exists(config.scenario), 'Scenario not found: %s' % config.scenario

    locator = cea.inputlocator.InputLocator(config.scenario, config.plugins)
    building_names = config.agriculture.buildings
    num_process = config.get_number_of_processes()
    n = len(building_names)

    # DLI calculations for each surface of each building
    cea.utilities.parallel.vectorize(calc_DLI, num_process)\
        (repeat(locator, n), repeat(config, n), building_names)

if __name__ == '__main__':
    main(cea.config.Configuration())
