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


class BiaAssessmentPlugin(cea.plugin.CeaPlugin):

    pass


def main(config):
    assert os.path.exists(config.scenario), 'Scenario not found: %s' % config.scenario
    locator = cea.inputlocator.InputLocator(config.scenario, config.plugins)

    # Check if the DLI calculations had been executed before proceeding to BIA Assessment
    building_names = config.agriculture.buildings
    num_process = config.get_number_of_processes()
    n = len(building_names)

    dir_dli = config.scenario + "/outputs/data/potentials/agriculture/dli"     # path of the directory
    dli_file = os.listdir(dir_dli)  # Getting the list of directories
    if len(dli_file) != n:      # if some or all the DLI results are not in-place
        print("(HINT: cancel, delete Folder outputs>data>potentials>agriculture, and start with DLI Calculation again)")
        exit()

    # BIA assessment
    type_crop = config.agriculture.type_crop
    print('Executing CEA Building-Integrated Agriculture (BIA) assessment for {type_crop}.'
          .format(type_crop=type_crop))

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
