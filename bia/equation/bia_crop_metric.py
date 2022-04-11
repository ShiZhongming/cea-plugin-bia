"""
This script calculates:
the crop yields (kg), carbon reduction (eCO2kg), CAPEX and OPEX (USD)
for the selected crop type on each building envelope surface.
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


# filter by crop on wall/roof/window defined in config.file

def calc_crop(locator, config, building_name):
    """
    This function calculates the series of metrics for the selected crop type.
    These metrics are crop yield in kg (existing),
    carbon emission reduction in CO2-equivalent (soon to come), CAPEX and OPEX in USD (soon to come).

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator

    :param building_name: list of building names in the case study
    :type building_name: Series

    :return: Building_crop.csv with BIA metrics for each building envelope surface.

    """

    t0 = time.perf_counter()

    # path to the results of DLI calculations using script: bia_dli.py
    dli_path = config.scenario + \
          "/outputs/data/potentials/agriculture/{building}_DLI.csv".format(building=building_name)

    # read the results of DLI calculations
    dli_results = pd.read_csv(dli_path)
    print('reading DLI simulation results done')

    # the selected crop type
    type_crop = config.agriculture.type_crop

    # read crop properties in the BIA database for the selected crop type
    crop_properties = calc_properties_crop_db(config)
    # print("gathering the properties of {type_crop}.".format(type_crop=type_crop))


    # spot the suitable surfaces for the selected crop type
    # this function calculates the surfaces and days for crop-growing 1 or 0, how many cycles per year on each surface

    cycle_total_annual = xxx['cycle'].sum()

    if not cycle_total_annual == 0:

        # Calculate the crop yields (kg) for the selected crop type on each building envelope surface.

        # Calculate the carbon emission (CO2-equivalent) for the selected crop type on each building envelope surface.

        # Calculate the cost (CAPEX and OPEX in USD) for the selected crop type on each building envelope surface.

        # merge the calculated BIA results
        sensors_metadata_clean_DLI_daily = pd.merge(sensors_metadata_clean, sensors_DLI_daily, left_index=True,
                                                    right_index=True, how="left")

        # write the BIA results
        output_path = config.scenario + \
          "/outputs/data/potentials/agriculture/{building}_BIA.csv".format(building=building_name)
        sensors_metadata_clean_DLI_daily.to_csv(output_path, index=True,
                                    float_format='%.2f',
                                    na_rep=0)

        print('BIA assessments for each surface on Building', building_name,
              'done - time elapsed: %.2f seconds' % (time.perf_counter() - t0))

    else:  # No surface meets the minimum DLI requirement of the selected crop type
        print("Unfortunately, {type_crop} is unlikely to grow on the BIA-permissible surfaces on the roof and facade "
              "of Building {building_name}".format(type_crop=type_crop, building_name=building_name))
        pass