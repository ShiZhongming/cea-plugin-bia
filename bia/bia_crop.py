"""
This script calculates:
the crop yields (kg), carbon emission reduction (CO2-equivalent), cost (CAPEX and OPEX in USD)
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


class BiaCropPlugin(cea.plugin.CeaPlugin):
    """
    Define the plugin class - unless you want to customize the behavior, you only really need to declare the class. The
    rest of the information will be picked up from ``default.config``, ``schemas.yml`` and ``scripts.yml`` by default.
    """
    pass


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
    if not os.path.exists(dli_path):
        raise ValueError("Consider execute Command bia-dli first!")

    # read the results of DLI calculations
    dli_results = pd.read_csv(dli_path)
    print('reading DLI simulation results done')

    # the selected crop type
    crop_type = config.agriculture.crop_type

    # read crop properties in the BIA database for the selected crop type
    bia_database_path = r"bia/bia_database.xlsx"
    crop_properties= calc_properties_crop_db(bia_database_path, config)
    print("gathering the properties of {crop_type}.".format(crop_type=crop_type))

    # spot the suitable surfaces for the selected crop type
    # this function calculates the surfaces and days for crop-growing 1 or 0, how many cycles per year on each surface
    xxx =
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
        print("Unfortunately, {crop_type} is unlikely to grow on the BIA-permissible surfaces on the roof and facade "
              "of Building {building_name}".format(crop_type=crop_type, building_name=building_name))
        pass


def calc_properties_crop_db(database_path, config):
    """
    To retrieve the crop properties stored in the BIA database for the selected crop type.

    :param crop: the selected crop type
    :type crop_type: string
    :return: dict with Properties of the crop type taken form the database
    """
    crop_type = config.agriculture.crop_type
    data = pd.read_excel(database_path, sheet_name="crop")
    crop_properties = data[data['crop_type'] == crop_type].reset_index().T.to_dict()[0]

    return crop_properties


def main(config):
    assert os.path.exists(config.scenario), 'Scenario not found: %s' % config.scenario
    locator = cea.inputlocator.InputLocator(config.scenario, config.plugins)

    print('Running BIA assessment with scenario = %s' % config.scenario)
    print('Running BIA assessment for the crop of {crop_type}'.format(crop_type=config.agriculture.crop_type))
    print('Running BIA assessment with crop-on-roof = %s' % config.agriculture.crop_on_roof)
    print('Running BIA assessment with crop-on-wall = %s' % config.agriculture.crop_on_wall)

    building_names = locator.get_zone_building_names()
    num_process = config.get_number_of_processes()
    n = len(building_names)
    cea.utilities.parallel.vectorize(calc_crop, num_process)(repeat(locator, n), repeat(config, n), building_names)


if __name__ == '__main__':
    main(cea.config.Configuration())

