"""
This script calculates:
the annual growth cycles
for the selected crop type on each building envelope surface.
"""

from __future__ import division
from __future__ import print_function
import os
import pandas as pd
import numpy as np
from cea.utilities import epwreader
from cea.inputlocator import InputLocator


__author__ = "Zhongming Shi"
__copyright__ = "Copyright 2023, A/S Group, ITA, ETH Zurich"
__credits__ = ["Zhongming Shi"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Zhongming Shi"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"


def calc_properties_crop_db(type_crop):
    """
    To retrieve the crop properties stored in the BIA database for the selected crop type.

    :param database_path: the selected crop type
    :type database_path: string
    :return: dict with properties of the selected crop type retrieved form the database
    :type dict
    """

    # path to the bia database
    dir = os.path.dirname(__file__)
    database_path = os.path.join(dir, "bia_data.xlsx")

    data = pd.read_excel(database_path, sheet_name="crop")
    crop_properties = data[data['type_crop'] == type_crop].reset_index().T.to_dict()[0]

    return crop_properties


def calc_properties_env_db(config):
    """
    To retrieve the environmental impacts of the selected crop type stored in the BIA database.

    :param database_path: the selected crop type
    :type database_path: string

    :return: DataFrame with properties of the environmental impact of the selected crop type retrieved form the database
    """

    # path to the bia database
    dir = os.path.dirname(__file__)
    database_path = os.path.join(dir, "bia_data.xlsx")

    type_crop = config.agriculture.type_crop
    data = pd.read_excel(database_path, sheet_name="env")
    # baseline scenario now only contains lettuce as an example; more to be added
    env_properties = data[data['type_crop'] == 'lettuce']

    return env_properties


def calc_properties_cost_db(config):
    """
    To retrieve the expenditures related the selected crop type stored in the BIA database.

    :param database_path: the selected crop type
    :type database_path: string

    :return: DataFrame with properties of expenditures related to the selected crop type retrieved form the database
    """

    # path to the bia database
    dir = os.path.dirname(__file__)
    database_path = os.path.join(dir, "bia_data.xlsx")

    type_crop = config.agriculture.type_crop
    data = pd.read_excel(database_path, sheet_name="cost")
    cost_properties = data[data['type_crop'] == type_crop]

    return cost_properties


def calc_chunk_day_crop(data):

    """
     To separate the chunks of dates into lists
     Then, the outputs can be used to calculate the growth cycles.

    :param data: the days' DLI meeting the requirement of the selected crop
    :type data: list
    :return: a list of chunks of dates (each chunk = one season) for each building surface
    :type consecutive_list: list
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


def calc_crop_cycle(config, building_name, type_crop):

    """
     To spot the days that are suitable for the selected crop type based on its crop property for each building
     envelope surface.
     To calculate the number of growth cycles for the selected crop type for each building envelope surface.

    :param crop_properties: the crop property of the selected crop type
    :type crop_properties: dict

    :return cycl_srf: number of cycles, including both initial and subsequent ones,
    for each building surface in each season of a whole year
    :type cycl_srf: list

    :return season_srf: number of days in each season (the beginning and the end of each year are looped as one single
    season and placed at the end of the list)
    :type season_srf: list

    :return date_srf: the days (0 to 364, in total 365 days in a non-leap year) that are eligible for growing the
    selected crop type
    :type date_srf: list

    :return cycl_i_srf: number of cycles, including initial ones only,
    for each building surface of a whole year
    :type cycl_i_srf: list

    :return cycl_s_srf: number of cycles, including subsequent ones only,
    for each building surface of a whole year
    :type cycl_s_srf: list

    """
    print("Reading and analysing the weather data.")
    locator = InputLocator(scenario=config.scenario)
    weather_data = epwreader.epw_reader(locator.get_weather_file())
    temp_C_h = weather_data['drybulb_C']

    # calculate from hourly temp_C_h to daily temp_C_d
    day = pd.Series(range(0, 365))
    hour_to_day = day.repeat(24).reset_index().pop('index')
    temp_C_h = pd.merge(temp_C_h, hour_to_day, left_index=True, right_index=True, how="left")
    temp_C_d = temp_C_h.groupby(['index']).sum().reset_index() / 24
    del temp_C_d['index']       # average Temperature_C per 365 days

    print("Calculating the planting information for {type_crop} on Building {building}."
          .format(type_crop=type_crop, building=building_name))

    # get the properties of the selected crop
    crop_properties = calc_properties_crop_db(type_crop)

    # unpack the properties of the selected crop type: cycle days
    cycl_i_day = int(crop_properties.get('cycl_i_day'))  # growth cycle in days: initial
    cycl_s_day = int(crop_properties.get('cycl_s_day'))  # growth cycle in days: subsequent
    n_cycl = int(crop_properties.get('n_cycl'))  # number of growth cycles: both initial and subsequent
    temp_opt_gro_l_c = int(crop_properties.get('temp_opt_gro_l_c'))     # preferred temperature for growth: lower bound
    temp_opt_gro_u_c = int(crop_properties.get('temp_opt_gro_u_c'))     # preferred temperature for growth: upper bound

    # unpack the properties of the selected crop type: DLI
    dli_l = float(crop_properties.get('dli_l'))    # DLI requirement: lower bound
    # dli_u = float(crop_properties.get('dli_u'))    # DLI requirement: upper bound

    # inputs to select the surface and calculate the number of growth cycles
    dli_criteria = dli_l   # DLI requirement for the selected crop
    temp_criteria_l = temp_opt_gro_l_c      # temperature requirement: lower bound
    temp_criteria_u = temp_opt_gro_u_c      # temperature requirement: upper bound

    # read the DLI results
    dli_path = config.scenario + "/outputs/data/potentials/agriculture/dli/{building}_DLI.csv"\
        .format(building=building_name)
    cea_dli_results = pd.read_csv(dli_path)

    # slice the DLIs/temps of the 365 days
    dli_365 = cea_dli_results.loc[:, "0": "364"]
    dli_365_add_one_cycl = pd.concat([dli_365, dli_365.iloc[:, 0:(cycl_i_day - 1)]], axis=1)
    temp_c_365_add_one_cycl = pd.concat([temp_C_d, temp_C_d.iloc[:, 0:(cycl_i_day - 1)]], axis=1)

    # number of this building's surfaces
    n_surface = len(dli_365.index)

    # spot the first days of potential growth cycles, for the growth of at least one initial cycle only
    day = []    # to store the first days of potential cycles
    bool_df = pd.DataFrame()
    for i in range(365):
        # for consecutive days equaling a growth cycle:
        # true if the average DLI is above the required DLI
        bool_dli = dli_365_add_one_cycl.iloc[:, i:(i + cycl_i_day)].sum(axis=1) >= dli_criteria * cycl_i_day
        # true if the average temp_C is above the lower bound of preferred temp
        bool_temp_l = temp_c_365_add_one_cycl.iloc[:, i:(i + cycl_i_day)].sum(axis=1) >= temp_criteria_l * cycl_i_day
        # true if the average temp_C is below the upper bound of preferred temp
        bool_temp_u = temp_c_365_add_one_cycl.iloc[:, i:(i + cycl_i_day)].sum(axis=1) <= temp_criteria_u * cycl_i_day
        # merge criteria
        mask = [all(tup) for tup in zip(bool_dli, bool_temp_l, bool_temp_u)]
        bool = pd.DataFrame(bool_dli, columns=['{i}'.format(i=i)])

        # record the boolean values for each day of all surfaces
        bool_df = pd.concat([bool_df, bool], axis=1)

    day_365 = pd.DataFrame(columns=range(365))
    day_365.loc[0] = range(365)
    surface_day_365 = pd.concat([day_365] * n_surface, ignore_index=True)  # Ignores the index

    first_day = surface_day_365[bool_df].values.tolist()    # using bool_df as a mask; from dataframe to lists of a list
    first_day = [[x for x in y if not np.isnan(x)] for y in first_day]  # remove the nan in each list

    # spot the seasons (of days) suitable to grow the selected crop for each surface
    day_crop = first_day
    for j in range(cycl_i_day):
        other_day = [[x + j for x in y] for y in first_day]
        day_crop = [x + y for x, y in zip(day_crop, other_day)]     # with repeated days

    day_crop = [list(set(x)) for x in day_crop]     # remove the repeated days

    # Calculate the length (in number of days) of each season
    season_srf = []
    date_srf = []

    for surface in range(n_surface):
        day_srf = day_crop[surface]
        day_srf.sort()  # sort the days if they are not in an ascending order

        if day_srf:
            excess_day = int(day_srf[-1]) - 364
            # this will be greater or equal to 0 if there are day numbers beyond 365 days of a year (non-leap year)

            # when the beginning and the end of a year need to be connected as one season
            if excess_day >= 0:
                # update the eligible days at the beginning of the year based on the excessive days
                day_srf = list(range(excess_day)) + day_srf

                # remove the dates that are equal or larger than 365
                day_srf = [x for x in day_srf if x < 365]

                day_srf = list(set(day_srf))    # remove the repeated days, if any
                season_crop = calc_chunk_day_crop(day_srf)
                len_season = [365 if len(x) > 365 else len(x) for x in season_crop]

                if len(len_season) >= 2:
                # remove the original first and last season: then add the merged season at the end of the list
                    merged_season = len_season[0] + len_season[-1] - excess_day
                    len_season.pop(-1)  # remove the original last season
                    len_season.pop(0)   # remove the original first season
                    len_season.append(merged_season)    # add the newly merged season at the end of the seasons

            else:
                season_crop = calc_chunk_day_crop(day_srf)
                len_season = [365 if len(x) > 365 else len(x) for x in season_crop]

        else:
            len_season = [0]

        season_srf.append(len_season)
        date_srf.append(day_srf)

    # Calculate the number of growth cycles for each building surface
    cycl_srf, cycl_i_srf, cycl_s_srf = calc_n_cycle_season(cycl_i_day, cycl_s_day, n_cycl, season_srf)

    return season_srf, cycl_srf, date_srf, cycl_i_srf, cycl_s_srf


def calc_n_cycle_season(cycl_i_day, cycl_s_day, n_cycl, season_srf):

    """
     This function calculates the number of growth cycles,
     including both initial and subsequent ones, for the selected crop, based on the length of each season.

    :param cycl_i_day: the initial cycle of the selected crop in days
    :type cycl_i_day: integer
    :param cycl_s_day: the subsequent cycle of the selected crop in days
    :type cycl_s_day: integer
    :param n_cycl: the limit of number of growth cycles of the selected crop in days
    :type n_cycl: integer
    :param season_srf:  lengths of season(s) of the selected crop in days for each building surface of a whole year
    :type season_srf: list

    :return cycl_srf: number of cycles, including both initial and subsequent ones,
    for each building surface of a whole year
    :type cycl_srf: list
    :return cycl_i_srf: number of cycles, including initial ones only,
    for each building surface of a whole year
    :type cycl_i_srf: list
    :return cycl_s_srf: number of cycles, including subsequent ones only,
    for each building surface of a whole year
    :type cycl_s_srf: list
    """

    tolerance_cycl = 0.05  # this tolerance allows some of the growth cycles a bit shorter than in the database
    crop_life = cycl_i_day + cycl_s_day * (n_cycl - 1)  # days of the full life of the selected crop
    n_season_srf = [len(x) for x in season_srf]     # number of seasons on each surface

    # calculate the cycles of the selected crop in the non-full-life times
    cycl_srf = []  # number of cycles, including both initial and subsequent ones
    cycl_i_srf = []  # number of cycles, including initial ones only
    cycl_s_srf = []  # number of cycles, including subsequent ones only

    for surface in range(len(season_srf)):
        day_season = season_srf[surface]
        n_season = n_season_srf[surface]

        # number of cycles when plant can grow full life
        n_full_life = [x // crop_life for x in day_season]
        n_cycl_season = [x * n_cycl for x in n_full_life]   # full life

        # remainder days for non-full life
        day_non_full_life = [x - crop_life * y for x, y in zip(day_season,
                                                               n_cycl_season)]

        for n in range(1, n_cycl):
            for season in range(n_season):
                if day_non_full_life[season] >= (cycl_i_day + (n - 1) * cycl_s_day)/crop_life/(1 + tolerance_cycl) \
                        and day_non_full_life[season] < (cycl_i_day + n * cycl_s_day)/crop_life/(1 + tolerance_cycl):
                    n_cycl_season[season] = n_cycl_season[season] + n

        cycl_srf.append(n_cycl_season)

        # number of initial cycles in each season
        n_i_cycl = [abs(-1 * (x // n_cycl)) for x in n_cycl_season]
        cycl_i_srf.append(n_i_cycl)

        # number of subsequent cycles in each season
        n_s_cycl = [x - y for x, y in zip(n_cycl_season, n_i_cycl)]
        cycl_s_srf.append(n_s_cycl)

    return cycl_srf, cycl_i_srf, cycl_s_srf
