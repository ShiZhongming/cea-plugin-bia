"""
This script creates:
the crop profile and planting calendar for each building surface, based on one of the following user-defined objectives:

environmental impacts including GHG Emissions (kg CO2-eq), energy (kWh) and water use (litre),
costs including capital and operational expenditures (USD)
for the selected crop type on each building envelope surface.
"""

from __future__ import division
from __future__ import print_function

import cea.config
import cea.inputlocator
import cea.plugin
import cea.utilities.parallel
from cea.constants import HOURS_IN_YEAR
from cea.resources.radiation_daysim import daysim_main, geometry_generator
from cea.utilities.standardize_coordinates import get_lat_lon_projected_shapefile
from cea.analysis.costs.equations import calc_capex_annualized, calc_opex_annualized

import os
import time
from itertools import repeat
from math import *
from multiprocessing import Pool
import pandas as pd
from geopandas import GeoDataFrame as gdf
import numpy as np

from bia.equation.bia_metric import calc_bia_metric
from bia.equation.bia_crop_cycle import calc_properties_crop_db, calc_chunk_day_crop, \
    calc_crop_cycle, calc_properties_env_db, calc_properties_cost_db, calc_n_cycle_season


__author__ = "Zhongming Shi"
__copyright__ = "Copyright 2022, Future Cities Laboratory, Singapore - ETH Zurich; " \
                "University of Calgary, Alberta, Canada"
__credits__ = ["Zhongming Shi"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Zhongming Shi"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"


# calculate the bia crop profiles for each surface and write to disk
def calc_bia_crop_profile(locator, config, building_name):

    """
    This function calculates the bia crop profiles for the one user-defined bia metric as an objective
     and write the results to disk.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator

    :param building_name: list of building names in the case study
    :type building_name: Series

    :return: no return; write the results to disc

    """

    t0 = time.perf_counter()

    # create the crop profile 365-calendar for each building surface
    bia_calendar_srf_df_to_filter = calc_crop_calendar(locator, config, building_name)

    # filter out the unwanted surfaces
    bia_calendar_srf_df_filtered = filter_crop_srf_profiler(locator, config, building_name,
                                                            bia_calendar_srf_df_to_filter)

    # write the results to disk
    # write to disk
    output_path = config.scenario + \
                  "/outputs/data/potentials/agriculture/{building}_BIA_crop_profile and planting calendar.csv" \
        .format(building=building_name)
    bia_calendar_srf_df_filtered.to_csv(output_path, index=False, na_rep=0)

    print('BIA crop profiles for each surface on Building', building_name,
          'have been created - time elapsed: %.2f seconds' % (time.perf_counter() - t0))


# ranks the crop types for each surface based on the user-defined BIA objective
def calc_crop_rank(locator, config, building_name):
    """
    This function ranks the crop types for each surface based on the user-defined BIA objective.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator

    :param building_name: list of building names in the case study
    :type building_name: Series

    :return: crop_rank_i_srf: the rank (index) of crop types for each surface based on the user defined BIA objective
    :type crop_rank_i_srf: list

    :return: crop_rank_srf: crop types (ranked) for each surface based on the user defined BIA objective
    :type crop_rank_srf: list

    """
    # get the user-defined list of crop types
    types_crop = config.crop_profile.types_crop
    # get the user-defined objective of bia-metric
    bia_metric_obj = config.crop_profile.bia_assessment_metric_objective

    # for each user-defined crop type
    bia_metric_obj_matrix = []
    for type_crop in range(len(types_crop)):

        # read the BIA metrics results for each building surface and each crop type
        bia_metrics_path = config.scenario + \
                           "/outputs/data/potentials/agriculture/surface/{building}_BIA_metrics_{type_crop}.csv"\
                               .format(building=building_name, type_crop=types_crop[type_crop])
        bia_metrics_srf = pd.read_csv(bia_metrics_path)

        # for each surface
        bia_metric_obj_srf = []
        for surface in range(len(bia_metrics_srf)):
            # calculate the user-fined BIA metric objective

            if bia_metric_obj == 'CropYield':
                the_bia_metric = bia_metrics_srf['yield_kg_per_year'].iloc[surface]
                bia_metric_obj_srf.append(the_bia_metric)

            elif bia_metric_obj == 'GHGEmission':
                the_bia_metric = bia_metrics_srf['ghg_kg_co2_bia'].iloc[surface]
                bia_metric_obj_srf.append(the_bia_metric)

            elif bia_metric_obj == 'EnergyUse':
                the_bia_metric = bia_metrics_srf['energy_kWh_bia'].iloc[surface]
                bia_metric_obj_srf.append(the_bia_metric)

            elif bia_metric_obj == 'WaterUse':
                the_bia_metric = bia_metrics_srf['water_l_bia'].iloc[surface]
                bia_metric_obj_srf.append(the_bia_metric)

            elif bia_metric_obj == 'AnnualisedCAPEX':
                the_bia_metric = bia_metrics_srf['capex_all_annualised_USD'].iloc[surface]
                bia_metric_obj_srf.append(the_bia_metric)

            elif bia_metric_obj == 'AnnualisedCAPEXPerKgYield':
                the_bia_metric = bia_metrics_srf['capex_all_annualised_USD'].iloc[surface] \
                                 / bia_metrics_srf['yield_kg_per_year'].iloc[surface]
                bia_metric_obj_srf.append(the_bia_metric)

            elif bia_metric_obj == 'AnnualOPEX':
                the_bia_metric = bia_metrics_srf['opex_all_USD_per_year'].iloc[surface]
                bia_metric_obj_srf.append(the_bia_metric)

            elif bia_metric_obj == 'AnnualOPEXPerKgYield':
                the_bia_metric = bia_metrics_srf['opex_all_USD_per_year'].iloc[surface] \
                                 / bia_metrics_srf['yield_kg_per_year'].iloc[surface]
                bia_metric_obj_srf.append(the_bia_metric)

            elif bia_metric_obj == 'AnnualCost':
                the_bia_metric = bia_metrics_srf['capex_all_annualised_USD'].iloc[surface] \
                                 + bia_metrics_srf['opex_all_USD_per_year'].iloc[surface]
                bia_metric_obj_srf.append(the_bia_metric)

            elif bia_metric_obj == 'AnnualCostPerKgYield':
                the_bia_metric = (bia_metrics_srf['capex_all_annualised_USD'].iloc[surface] \
                                 + bia_metrics_srf['opex_all_USD_per_year'].iloc[surface]) \
                                 / bia_metrics_srf['yield_kg_per_year'].iloc[surface]
                bia_metric_obj_srf.append(the_bia_metric)

            else:
                print('Error: user to define the objective BIA metric.')
                break

        bia_metric_obj_matrix.append(bia_metric_obj_srf)

    # store the results into DataFrame
    bia_metric_obj_matrix_df = pd.DataFrame(bia_metric_obj_matrix).T
    bia_metric_obj_matrix_df.columns = types_crop
    bia_metric_obj_matrix_df = bia_metric_obj_matrix_df.reset_index(drop=True)

    # create the rank of crop types for each surface
    # first, get the indices i that would sort the genres in descending (for BIA metrics that max is preferred)
    if bia_metric_obj == 'CropYield':
        i = np.argsort(bia_metric_obj_matrix_df.to_numpy() * -1, axis=1)
    else:   # ascending (for BIA metrics that min is preferred)
        i = np.argsort(bia_metric_obj_matrix_df.to_numpy() * 1, axis=1)
    # second, create a new DataFrame
    bia_crop_matrix_df = pd.DataFrame(bia_metric_obj_matrix_df.columns[i], columns=range(1, i.shape[1] + 1))
    # bia_crop_matrix_df.add_prefix('Rank')
    crop_rank_type_srf = bia_crop_matrix_df.values.tolist()
    crop_rank_i_srf = i.tolist()

    return crop_rank_i_srf, crop_rank_type_srf


# rank the crop types for each surface based on the user-defined BIA objective
def calc_crop_calendar(locator, config, building_name):
    """
    This function ranks the crop types for each surface based on the user-defined BIA objective.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator

    :param building_name: list of building names in the case study
    :type building_name: Series

    :return: bia_profile_srf_df: calendars of 365 days with each day specified the crop type to grow
    for (all) each building surface
    :type bia_profile_srf_df: DataFrame

    """

    print('Creating the crop profiles and planting calendars...')

    # get the user-defined list of crop types
    types_crop = config.crop_profile.types_crop
    # get the user-defined objective of bia-metric
    # bia_metric_obj = config.crop_profile.bia_assessment_metric_objective

    # to create the ranking crop types
    # for each surface based on the user-defined BIA objective
    crop_rank_i_srf, crop_rank_type_srf = calc_crop_rank(locator, config, building_name)  # a list and a list

    n_surface = len(crop_rank_i_srf)    # number of surface in this building
    n_crops_type = len(types_crop)      # number of candidate crop types

    # to create BIA calendars for each crop type on each building surface (without rank)
    date_srf_all = []
    for type_crop in range(n_crops_type):        # for each crop type

        # list: the days (0 to 364, in total 365 days in a non-leap year)
        # that are eligible for growing the crop type
        _, _, date_srf, _, _ = calc_crop_cycle(config, building_name, types_crop[type_crop])
        date_srf_all.append(date_srf)   # append all lists of dates into a single list

    # to create the lists of dates of each crop type for each building surface, not ranked
    crop_profile_srf_not_ranked = []
    for surface in range(n_surface):
        crop_profile_not_ranked = [x[surface] for x in date_srf_all]
        crop_profile_srf_not_ranked.append(crop_profile_not_ranked)

    # to re-order the lists of dates of each crop type for each building surface, ranked
    crop_profile_srf_ranked = []
    for surface in range(n_surface):
        crop_profile_not_ranked = crop_profile_srf_not_ranked[surface]
        rank_index = crop_rank_i_srf[surface]
        crop_profile_ranked = [x for _, x in sorted(zip(rank_index, crop_profile_not_ranked))]
        crop_profile_srf_ranked.append(crop_profile_ranked)

    # create a separate calendar for the Nth preferred crop type rank
    calendar_to_merge = []      # a list of DataFrame
    for n in range(n_crops_type):

        # create an empty DataFrame (calendar)
        calendar_to_fill = pd.DataFrame(columns=range(365))
        calendar_to_fill.insert(loc=0, column='srf', value=range(n_surface))

        for surface in range(n_surface):

            calendar_to_fill.at[surface, crop_profile_srf_ranked[surface][n]] = types_crop[n]
            calendar_to_fill.fillna('', inplace=True)


        calendar_to_merge.append(calendar_to_fill)  # a list

    # merge the calendars into a single DataFrame (calendar)
    # if a day is eligible for more than one crop type, the crop types are linked with a comma in the preferred rank
    bia_calendar_srf_df = pd.concat(calendar_to_merge, axis=0)
    bia_calendar_srf_df = bia_calendar_srf_df.groupby('srf')\
        .agg(lambda x: '/'.join(x.astype(str)))\
        .reset_index(drop=True)

    # indicate the no planting day/surface
    to_replace = '/' * (n_crops_type - 1)
    bia_calendar_srf_df = bia_calendar_srf_df.replace(to_replace, 'no planting')

    # get building surface info and join with the planting calendar
    dli_path = config.scenario + "/outputs/data/potentials/agriculture/{building}_DLI_daily.csv" \
        .format(building=building_name)
    cea_dli_results = pd.read_csv(dli_path)
    info_srf_df = cea_dli_results.loc[:, ['srf_index', 'AREA_m2', 'TYPE', 'wall_type']]
    bia_profile_srf_df = pd.concat([info_srf_df, bia_calendar_srf_df], axis=1)


    return bia_profile_srf_df


# filter by crop on wall/roof/window user-defined in config.file
def filter_crop_srf_profiler(locator, config, building_name, bia_metric_srf_df):
    """
    This function links the BIA metrics with surface info,
    then filters out the calculate metrics of the surfaces that are not intended by the user to have BIA.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator

    :param building_name: list of building names in the case study
    :type building_name: Series

    :return: DataFrame with BIA crop profiles and planting calendar for each surface intended to have BIA.

    """


    # get the user inputs
    bool_roof = config.crop_profile.crop_on_roof
    bool_window = config.crop_profile.crop_on_window
    bool_wall_u = config.crop_profile.crop_on_wall_under_window
    bool_wall_b = config.crop_profile.crop_on_wall_between_window

    # create the mask based the user input
    mask_df = pd.DataFrame(columns=['mask_roof', 'mask_window',
                                    'mask_wall_upper', 'mask_wall_lower',
                                    'mask_wall_b', 'mask'])
    mask_df['mask_roof'] = [bool_roof if x == 'roofs' else 0 for x in bia_metric_srf_df['TYPE'].tolist()]
    mask_df['mask_window'] = [bool_window if x == 'windows' else 0 for x in bia_metric_srf_df['TYPE'].tolist()]
    mask_df['mask_wall_upper'] = [bool_wall_u if x == 'upper' else 0 for x in bia_metric_srf_df['wall_type'].tolist()]
    mask_df['mask_wall_lower'] = [bool_wall_u if x == 'lower' else 0 for x in bia_metric_srf_df['wall_type'].tolist()]
    mask_df['mask_wall_b'] = [bool_wall_b if x == 'side' else 0 for x in bia_metric_srf_df['wall_type'].tolist()]
    mask_df['mask'] = [True if a == True or b == True or c == True or d == True or e == True else False
                       for a, b, c, d, e in zip(mask_df['mask_roof'],
                                                mask_df['mask_window'],
                                                mask_df['mask_wall_upper'],
                                                mask_df['mask_wall_lower'],
                                                mask_df['mask_wall_b']
                                                )]

    # remove the unwanted surfaces
    bia_metric_srf_df_filtered = bia_metric_srf_df[mask_df['mask']]

    return bia_metric_srf_df_filtered