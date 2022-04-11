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
import numpy as np

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

    # dscp = crop_properties.get('description')   # description of the selected crop type
    # temp_opt_ger_l_c = crop_properties.get('temp_opt_ger_l_c')  # optimal germination Celsius temperature: lower bound
    # temp_opt_ger_u_c = crop_properties.get('temp_opt_ger_u_c')  # optimal germination Celsius temperature: upper bound
    # temp_opt_gro_l_c = crop_properties.get('temp_opt_gro_l_c')  # optimal growth Celsius temperature: lower bound
    # temp_opt_gro_l_c = crop_properties.get('temp_opt_gro_l_c')  # optimal germination Celsius temperature: upper bound
    # temp_por_l_c = crop_properties.get('temp_por_l_c')  # poor plant growth Celsius temperature: lower bound
    # temp_por_u_c = crop_properties.get('temp_por_u_c')  # poor plant growth Celsius temperature: upper bound
    # cycl_l_day = crop_properties.get('cycl_l_day')  # growth cycle in days: lower bound
    # cycl_u_day = crop_properties.get('cycl_u_day')  # growth cycle in days: upper bound
    # humd_l = crop_properties.get('humd_l')  # suitable humidity: lower bound
    # humd_u = crop_properties.get('humd_u')  # suitable humidity: upper bound
    # dli_l = crop_properties.get('dli_l')    # DLI requirement: lower bound
    # dli_u = crop_properties.get('dli_u')    # DLI requirement: upper bound
    # yld_grd_l_kg_sqm = crop_properties.get('yld_grd_l_kg_sqm')  # yield on ground in kilograms: lower bound
    # yld_grd_u_kg_sqm = crop_properties.get('yld_grd_u_kg_sqm')  # yield on ground in kilograms: upper bound
    # yld_bia_l_kg_sqm = crop_properties.get('yld_bia_l_kg_sqm')  # BIA yield in kilograms: lower bound
    # yld_bia_u_kg_sqm = crop_properties.get('yld_bia_u_kg_sqm')  # BIA yield in kilograms: upper bound
    # mkt_sg_sgd_kg = crop_properties.get('mkt_sg_sgd_kg')    # market price in Singapore: SGD 14.75/kg on 12 Dec 2021 at Lazada

    return crop_properties

def calc_chunk_day_crop(data):

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
            if data[chunks + 1] - data[chunks] == 1:
                if data[chunks] not in consecutive_list:
                    consecutive_list.append(data[chunks])
                consecutive_list.append(data[chunks + 1])
            else:
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

    print("Running building", building_name)

    # get the properties of the selected crop
    crop_properties = calc_properties_crop_db(config)

    # unpack the properties of the selected crop type: cycle days
    cycl_i_day = int(crop_properties.get('cycl_i_day'))  # growth cycle in days: initial
    cycl_s_day = int(crop_properties.get('cycl_s_day'))  # growth cycle in days: subsequent
    n_cycl = int(crop_properties.get('n_cycl'))  # number of growth cycles: both initial and subsequent

    # unpack the properties of the selected crop type: DLI
    dli_l = float(crop_properties.get('dli_l'))    # DLI requirement: lower bound
    dli_u = float(crop_properties.get('dli_u'))    # DLI requirement: upper bound

    # inputs to select the surface and calculate the number of growth cycles
    dli_criteria = (dli_l + dli_u) / 2   # DLI requirement for the selected crop

    # read the daily DLI results
    dli_path = config.scenario + "/outputs/data/potentials/agriculture/{building}_DLI_daily.csv"\
        .format(building=building_name)
    cea_dli_results = pd.read_csv(dli_path)

    # slice the DLIs of the 365 days
    dli_365 = cea_dli_results.loc[:, "0": "364"]
    dli_365_add_one_cycl = pd.concat([dli_365, dli_365.iloc[:, 0:(cycl_i_day - 1)]], axis=1)

    # number of this building's surfaces
    n_surface = len(dli_365.index)

    # spot the first days of potential growth cycles, for the growth of at least one initial cycle only
    day = []    # to store the first days of potential cycles
    bool_df = pd.DataFrame()
    for i in range(365):
        # for a consecutive days equaling a growth cycle, true if the average DLI is below the required DLI
        bool = dli_365_add_one_cycl.iloc[:, i:(i + cycl_i_day)].sum(axis=1) >= dli_criteria * cycl_i_day
        bool = bool.to_frame(name=i)

        # record the boolean values for each day of all surfaces
        bool_df = pd.concat([bool_df, bool], axis=1)

    day_365 = pd.DataFrame(columns=range(365))
    day_365.loc[0] = range(365)
    surface_day_365 = pd.concat([day_365] * n_surface, ignore_index=True)  # Ignores the index
    # temp = surface_day_365[bool_df]
    # out = r'/Users/shi_zhongming/Dropbox/Calgary_U/Winter2022/Arch612/cea/Alberta_solar/BIA_test/outputs/data/
    # potentials/agriculture/{building}_first_day.csv'.format(building=building_name)
    # # temp.to_csv(out)
    # # print(building_name, temp)
    first_day = surface_day_365[bool_df].values.tolist()    # using bool_df as a mask; from dataframe to lists of a list
    first_day = [[x for x in y if not np.isnan(x)] for y in first_day]  # remove the nan in each list

    # spot the seasons (of days) suitable to grow the selected crop for each surface
    day_crop = first_day
    for j in range(cycl_i_day):
        other_day = [[x + j for x in y] for y in first_day]
        day_crop = [x + y for x, y in zip(day_crop, other_day)]     # with repeated days

    day_crop = [list(set(x)) for x in day_crop]     # remove the repeated days

    # Calculate the length of each season
    season_srf = []
    for surface in range(n_surface):
        # print("surface", surface)
        season_crop = calc_chunk_day_crop(day_crop[surface])
        len_season = [365 if len(x) > 365 else len(x) for x in season_crop]
        season_srf.append(len_season)

    # print('season_srf', season_srf)
    # print('len_season_srf', len(season_srf))
    # print('n_surface', n_surface)

    # Calculate the number of growth cycles for each building surface
    cycl_srf = calc_n_cycle_season(cycl_i_day, cycl_s_day, n_cycl, season_srf)

    return cycl_srf


def calc_n_cycle_season(cycl_i_day, cycl_s_day, n_cycl, season_srf):

    """
     This function calculates the number of growth cycles,
     including both initial and subsequent ones, for the selected crop.

    :param cycl_i_day: the initial cycle of the selected crop in days
    :type cycl_i_day: integer
    :param cycl_s_day: the subsequent cycle of the selected crop in days
    :type cycl_s_day: integer
    :param n_cycl: the limit of number of growth cycles of the selected crop in days
    :type n_cycl: integer
    :param season_srf:  lengths of season(s) of the selected crop in days for each building surface of a whole year
    :type season_srf: list

    :return season_srf: number of cycles, including both initial and subsequent ones,
    for each building surface of a whole year
    :type season_srf: list
    """

    tolerance_cycl = 0.05  # this tolerance allows some of the growth cycles a bit shorter than on the paper
    len_season_all = [item for sublist in season_srf for item in sublist]  # length of all seasons in a building

    # n_cycle_season = floor([len(x) for x in season_crop] / (cycl_day / (1 + tolerance_cycl)))

    return cycl_srf