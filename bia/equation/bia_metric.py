"""
This script calculates:
the crop yields (kg), GHG emissions (eCO2kg), CAPEX and OPEX (USD)
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

from bia.equation.bia_crop_cycle import calc_properties_crop_db
from bia.equation.bia_crop_cycle import calc_crop_cycle


__author__ = "Zhongming Shi"
__copyright__ = "Copyright 2022, Future Cities Laboratory, Singapore - ETH Zurich; " \
                "University of Calgary, Alberta, Canada"
__credits__ = ["Zhongming Shi"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Zhongming Shi"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"


# calculate the BIA metrics and write to disk
def calc_bia_metric(locator, config, building_name):

    """
    This function calculates the three categories of BIA metrics and write the results to disk.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator

    :param building_name: list of building names in the case study
    :type building_name: Series

    :return: no return

    """



# calculate crop yields for each surface
def calc_crop_yields(locator, config, building_name, pest_rate, ):

    """
    This function calculates crop yield in kg for each building envelope surface.
    At the moment, no pest harm is considered.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator

    :param building_name: list of building names in the case study
    :type building_name: Series

    :param pest_rate: ratio of yields devoured by pest, ranging from 0 to 1 with 1 being crops all eaten by pests.
    :type pest_rate: float

    :return: Dataframe with crop_yields in kg for each building envelope surface.

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

    return crop_yields


# calculate ghg emissions for each surface, crop as a produce
def calc_crop_ghg(locator, config, building_name):

    """
    This function calculates the GHG emissions (eCO2kg) for each surface.
    Two categories of GHG emissions are considered: transport and carbon fixation.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator

    :param building_name: list of building names in the case study
    :type building_name: Series

    :return: Dataframe with GHG emissions in eCO2kg for each building envelope surface.

    """

    return crop_ghg


# calculate hourly carbon fixation for each surface, crop as a plant
def calc_carbon_fixation(locator, config, building_name):
    """
    This function calculates the hourly carbon fixation (eCO2kg) for each surface.
    This fuction is dependent on Zhang Qianning and Huang Zhaolu's work. To be confirmed with NUS.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator

    :param building_name: list of building names in the case study
    :type building_name: Series

    :return: Dataframe with hourly carbon fixation in eCO2kg for each building envelope surface.

    """

    return carbon_fixation


# calculate CAPEX for each surface
def calc_crop_capex(locator, config, building_name):
    """
    This function calculates the CAPEX in USD (infrastructure) for each surface.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator

    :param building_name: list of building names in the case study
    :type building_name: Series

    :return: DataFrame with CAPEX (infrastructure) for each building envelope surface.

    """

    return crop_capex


# calculate OPEX for each surface
def calc_crop_cost(locator, config, building_name):
    """
    This function calculates the OPEX in USD (seed, pesticide, and maintenance) for each surface.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator

    :param building_name: list of building names in the case study
    :type building_name: Series

    :return: DataFrame with OPEXes (seed and pesticide) for each building envelope surface.

    """

    return crop_opex


# calculate carbon footprint (transport) for each kg of crop imported to Singapore and residences
def calc_crop_carbon_footprint(crop_yields, config):
    """
    This function calculates the GHG emissions (eCO2kg) for the crop yields if used in transport for each surface.
    Two categories of GHG emissions are considered: transport and carbon fixation.

    :param crop_yields: crop yields in kg for each building surface
    :type crop_yields: dataframe

    :return: Dataframe with GHG emissions in eCO2kg for the crop yields on each building envelope surface.

    """
    return crop_carbon_footprint



# filter by crop on wall/roof/window user-defined in config.file
def filter_crop_srf(locator, config, building_name):

    """
    This function filters out the calculate metrics of the surfaces that are not intended by the user to have BIA.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator

    :param building_name: list of building names in the case study
    :type building_name: Series

    :return: DataFrame with BIA metrics for each surfaces intended to have BIA.

    """

    # to retrieve the crop properties stored in the BIA database for the selected crop type
    crop_properties = calc_properties_crop_db(config)

    # get the user inputs
    type_crop = config.agriculture.type_crop
    bool_roof = config.agriculture.crop_on_roof
    bool_window = config.agriculture.crop_on_window
    bool_wall_u = config.agriculture.crop_on_wall_under_window
    bool_wall_b = config.agriculture.crop_on_wall_between_window

    # get the annual growth cycles for the selected crop type on each building envelope surface


    return crop_srf_filtered

