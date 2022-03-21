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
    crop_type = config.agriculture.crop_type

    # read crop properties in the BIA database for the selected crop type
    crop_properties = calc_properties_crop_db(config)
    # print("gathering the properties of {crop_type}.".format(crop_type=crop_type))


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
        print("Unfortunately, {crop_type} is unlikely to grow on the BIA-permissible surfaces on the roof and facade "
              "of Building {building_name}".format(crop_type=crop_type, building_name=building_name))
        pass


def calc_properties_crop_db(config):
    """
    To retrieve the crop properties stored in the BIA database for the selected crop type.

    :param database_path: the selected crop type
    :type database_path: string
    :return: dict with properties of the selected crop type retrieved form the database
    """

    # path to the bia database
    dir = os.path.dirname(__file__)
    database_path = os.path.join(dir, "bia_data.xlsx")

    crop_type = config.agriculture.crop_type
    data = pd.read_excel(database_path, sheet_name="crop")
    crop_properties = data[data['crop_type'] == crop_type].reset_index().T.to_dict()[0]

    return crop_properties

def calc_surface_crop_cycle(dli_results, crop_properties, config):

    """
     To spot the days that are suitable for the selected crop type based on its crop property for each building
     envelope surface.
     To calculate the number of growth cycles for the selected crop type for each building envelope surface.

    :param dli_results: the dli on each building envelope surface for 365 days
    :type dli_results: dataframe
    :param crop_properties: the crop property of the selected crop type
    :type crop_properties: dict
    :return: dataframe of each building envelope surface, the days to be utilised
    and the number of growth cycles per year
    """

    # unpack the properties of the selected crop type
    dscp = crop_properties.get('description')   # description of the selected crop type
    temp_opt_ger_l_c = crop_properties.get('temp_opt_ger_l_c')  # optimal germination Celsius temperature: lower bound
    temp_opt_ger_u_c = crop_properties.get('temp_opt_ger_u_c')  # optimal germination Celsius temperature: upper bound
    temp_opt_gro_l_c = crop_properties.get('temp_opt_gro_l_c')  # optimal growth Celsius temperature: lower bound
    temp_opt_gro_l_c = crop_properties.get('temp_opt_gro_l_c')  # optimal germination Celsius temperature: upper bound
    temp_por_l_c = crop_properties.get('temp_por_l_c')  # poor plant growth Celsius temperature: lower bound
    temp_por_u_c = crop_properties.get('temp_por_u_c')  # poor plant growth Celsius temperature: upper bound
    cycl_l_day = crop_properties.get('cycl_l_day')  # growth cycle in days: lower bound
    cycl_u_day = crop_properties.get('cycl_u_day')  # growth cycle in days: upper bound
    humd_l = crop_properties.get('humd_l')  # suitable humidity: lower bound
    humd_u = crop_properties.get('humd_u')  # suitable humidity: upper bound
    dli_l = crop_properties.get('dli_l')    # DLI requirement: lower bound
    dli_u = crop_properties.get('dli_u')    # DLI requirement: upper bound
    yld_grd_l_kg_sqm = crop_properties.get('yld_grd_l_kg_sqm')  # yield on ground in kilograms: lower bound
    yld_grd_u_kg_sqm = crop_properties.get('yld_grd_u_kg_sqm')  # yield on ground in kilograms: upper bound
    yld_bia_l_kg_sqm = crop_properties.get('yld_bia_l_kg_sqm')  # BIA yield in kilograms: lower bound
    yld_bia_u_kg_sqm = crop_properties.get('yld_bia_u_kg_sqm')  # BIA yield in kilograms: upper bound
    mkt_sg_sgd_kg = crop_properties.get('mkt_sg_sgd_kg')    # market price in Singapore: SGD 14.75/kg on 12 Dec 2021 at Lazada

    #

    return crop_surface