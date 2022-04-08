"""
This script calculates:
the annual growth cycles
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

    type_crop = config.agriculture.type_crop
    data = pd.read_excel(database_path, sheet_name="crop")
    crop_properties = data[data['type_crop'] == type_crop].reset_index().T.to_dict()[0]

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

    return crop_properties

def calc_chunk_dli_days(data):

    """
     To separate the chunks of dates into lists
     Then, the outputs can be used to calculate the growth cycles.

    :param data: the days' DLI meeting the requirement of the selected crop
    :type data: list
    :return: dataframe of each building envelope surface, the days to be utilised
    and the number of growth cycles per year
    """

    consecutive_list = []

    for chunks in range(len(data)):

        try:
            #check consecutiveness
            if data[chunks + 1] - data[chunks] == 1:

                #check if it's already in list
                if data[chunks] not in consecutive_list:
                    consecutive_list.append(data[chunks])

                #add last one too
                consecutive_list.append(data[chunks + 1])

            else:
                #yield here and empty list
                yield consecutive_list
                consecutive_list = []
        except Exception:
            pass
    yield consecutive_list

def calc_crop_cycle(config, building_name):

    """
     To spot the days that are suitable for the selected crop type based on its crop property for each building
     envelope surface.
     To calculate the number of growth cycles for the selected crop type for each building envelope surface.

    :param crop_properties: the crop property of the selected crop type
    :type crop_properties: dict
    :return: dataframe of each building envelope surface, the days to be utilised
    and the number of growth cycles per year
    """
    # get the properties of the selected crop
    crop_properties = calc_properties_crop_db(config)

    # unpack the properties of the selected crop type: cycle days
    cycl_l_day = int(crop_properties.get('cycl_l_day'))  # growth cycle in days: lower bound
    cycl_u_day = int(crop_properties.get('cycl_u_day'))  # growth cycle in days: upper bound

    # unpack the properties of the selected crop type: DLI
    dli_l = float(crop_properties.get('dli_l'))    # DLI requirement: lower bound
    dli_u = float(crop_properties.get('dli_u'))    # DLI requirement: upper bound

    # inputs to select the surface and calculate the number of growth cycles
    cycl_day = int((cycl_l_day + cycl_u_day) / 2)   # growth cycle for the selected crop in days
    dli_criteria = (dli_l + dli_u) / 2   # DLI requirement for the selected crop

    # read the daily DLI results
    dli_path = config.scenario + "/outputs/data/potentials/agriculture/{building}_DLI_daily.csv"\
        .format(building=building_name)
    cea_dli_results = pd.read_csv(dli_path)

    # slice the DLIs of the 365 days
    dli_365 = cea_dli_results.loc[:, "0": "364"]
    dli_365_one_cycl = pd.concat([dli_365, dli_365.iloc[:, 0:(cycl_day - 1)]], axis=1)

    # number of this building's surfaces
    n_surface = len(cea_dli_results.index)

    # spot the days that are eligible for growing the selected crop
    day = []    # to store the first days of potential cycles
    bool_df = pd.DataFrame()
    for i in range(365):
        # for a consecutive days equaling a growth cycle, true if the average DLI is below the required DLI
        bool = dli_365_one_cycl.iloc[:, [i, i + cycl_day - 1]].sum(axis=1) < dli_criteria * cycl_day
        # record the boolean values for each day of all surfaces
        bool_df = pd.concat([bool_df, bool], axis=1)

    # spot the first days of potential growth cycles
    day_365 = pd.DataFrame(columns=range(365))
    day_365.loc[0] = range(365)
    surface_day_365 = pd.concat([day_365] * n_surface, ignore_index=True)  # Ignores the index
    print("yyyyyyyy", surface_day_365)
    first_day = surface_day_365[bool_df].fillna(999)

    # spot the seasons (of days) suitable to grow the selected crop for each surface
    season = pd.DataFrame()
    for j in range(cycl_day):
        first_day += j
        season = pd.concat([season, first_day], axis=1)
    season_surface_with_duplicates = season.values.tolist()

    season_surface_without_duplicates = []
    for surface in len(n_surface):
        season_surface_without_duplicates.append(list(dict.fromkeys(season_surface_with_duplicates[surface])))

    # Calculate the number of cycles
    n_cycle = []
    tolerance_cycl = 0.05   # this tolerance allows some of the growth cycles a bit shorter than on the paper
    for surface in len(n_surface):
        season_surface_list = calc_chunk_dli_days(season_surface_without_duplicates[surface])
        # number of cycles (with tolerance) in each season
        n_cycle_season = math.floor([len(x) for x in season_surface_list] / cycl_day / (1 + tolerance_cycl))
        # total cycles throughout the year for each surface
        n_cycle.append(sum(n_cycle_season))

    return n_cycle


