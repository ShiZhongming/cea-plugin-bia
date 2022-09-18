"""
This script creates:
.csv files as the inputs to be visualised in calendar heatmap graphs and dashboards.
.
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
from csv import writer

import pandas as pd

import cea.utilities.parallel
from cea.constants import HOURS_IN_YEAR
from cea.resources.radiation_daysim import daysim_main, geometry_generator
from bia.equation.bia_dli import calc_DLI
from bia.equation.bia_crop_cycle import calc_crop_cycle
from bia.equation.bia_metric import calc_bia_metric
from bia.equation.bia_select import calc_bia_crop_profile

__author__ = "Zhongming Shi"
__copyright__ = "Copyright 2022, Future Cities Laboratory, Singapore - ETH Zurich; " \
                "University of Calgary, Alberta, Canada"
__credits__ = ["Zhongming Shi"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Zhongming Shi"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"


# generate .csv files to be used as the input for bia visualisation
def calc_bia_visual(locator, config, building_name):

    """
     This function generates .csv files to be used as the input for bia visualisation
     for each crop type's planting calendar
     for all crop types combined's planting calendar
     for all information to be included in the dashboard
     and write to disk

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator
    :param building_name: list of building names in the case study
    :type building_name: Series

    :return:n/a
    """

    t0 = time.perf_counter()

    # create the folder to store the BIA visual files
    dir = config.scenario + "/outputs/data/potentials/agriculture/plot"
    if not os.path.exists(dir):
        os.mkdir(dir)

    # generate each crop type's planting calendar, all crop types combined's planting calendar
    # then write to disk
    calc_crop_calendar_each_by_orie_floo(locator, config, building_name)

    # generate all information to be included in the dashboard and write to disk
    calc_crop_info_each_crop_by_orie_floo(locator, config, building_name)

    #print a message when all the .csv files have been generated
    print('All the .csv files for', building_name,
          'for BIA visualisation have been generated. - time elapsed: %.2f seconds' % (time.perf_counter() - t0))


# generate each crop type's planting calendar and write to disk
def calc_crop_calendar_each_by_orie_floo(locator, config, building_name):

    """
     This function generates .csv files to be used as the input for bia visualisation
     for each crop type's planting calendar
     and write to disk

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator
    :param building_name: list of building names in the case study
    :type building_name: Series

    :return:n/a
    """

    # read the daily DLI results
    dli_path = config.scenario + "/outputs/data/potentials/agriculture/{building}_DLI_daily.csv" \
        .format(building=building_name)
    cea_dli_results = pd.read_csv(dli_path)

    # gather the floor number and orientation information [e(ast), w(est), s(outh), n(orth)] for each building surface
    info_srf_df = cea_dli_results.loc[:, ['srf_index', 'orientation', 'n_floor']]

    # read the overall planting calendar
    crop_profile_path = config.scenario + \
                  "/outputs/data/potentials/agriculture/{building}_BIA_crop_profile and planting calendar.csv" \
        .format(building=building_name)
    crop_profile_df = pd.read_csv(crop_profile_path)

    # all the crop types to be plotted
    types_crop = config.crop_profile.types_crop

    # create an empty list to store the calendars for each crop type
    calendar_list = []

    # loop for each crop type
    for type_crop in range(len(types_crop)):

        # create a DataFrame to store the outcome
        # all surfaces are not grouped and kept as they are
        srf_all_df = info_srf_df

        # create a DataFrame to store the outcome
        # all surfaces are grouped by floor number and orientation
        srf_df = info_srf_df\
            .groupby(['orientation', 'n_floor'])\
            .agg({"srf_all": "sum"}).reset_index()

        for day in range(365):

            # label the days that are suitable to grow such crop type
            srf_all_df[day] = [1 if type_crop in x else 0 for x in crop_profile_df[day]]

            # create the DataFrame of eligible surfaces
            # grouped by floor number and orientation
            srf_eligible_df = srf_all_df[srf_all_df[day] == 1]\
                .groupby(['orientation', 'n_floor'])\
                .agg({"srf_eligible": "sum"}).reset_index()

            # add the column of eligible surfaces grouped by floor number and orientation
            srf_df = pd.merge(srf_df, srf_eligible_df, how='left',
                              left_on=['orientation', 'n_floor'],
                              right_on=['orientation', 'n_floor']
                              )

            # calculate the number of surfaces the surfaces that are suitable to grow such crop type
            # for the same floor number and orientation
            srf_df['suit'] = srf_all_df[day]\
                .groupby(['orientation', 'n_floor'])\
                .sum()

            # calculate the number of surfaces
            # for the same floor number and orientation
            srf_df['count'] = srf_all_df[day]\
                .groupby(['orientation', 'n_floor'])\
                .count()

            # calculate the percentage of surfaces suitable to grow such crop type
            # for the same floor number and orientation
            srf_df[day] = srf_df['suit'] / srf_df['count']

            # drop the 'suit' and 'count' columns
            srf_df = srf_df.drop(columns=['suit', 'count'])

        # append each crop type's calendar into a list
        calendar_list.append(srf_all_df)

        # process the DataFrame (each crop type) to the needed format
        date = pd.date_range('1/1/2022', periods=365, freq='D').strftime("%Y-%m-%d").tolist()
        formatted_srf_df = srf_df.T.reset_index(drop=True)
        formatted_srf_df.insert(0, 'date', date)

        # write the outcome (each crop type) to disk
        output_path = config.scenario + \
                           "/outputs/data/potentials/agriculture/plot/{building}_BIA_visual_{type_crop}.csv" \
                               .format(building=building_name, type_crop=types_crop[type_crop])
        formatted_srf_df.to_csv(output_path, index=False, na_rep=0)

    # create a combined calendar DataFrame
    all_srf_df = pd.concat(calendar_list).groupby('srf_index')[range(365)].sum().reset_index()

    # process the DataFrame (all crop types) to the needed format
    date = pd.date_range('1/1/2022', periods=365, freq='D').strftime("%Y-%m-%d").tolist()
    formatted_all_srf_df = all_srf_df.T.reset_index(drop=True)
    formatted_all_srf_df.insert(0, 'date', date)

    # write the outcome (all crop types) to disk
    output_path = config.scenario + \
                  "/outputs/data/potentials/agriculture/plot/{building}_BIA_visual_all_crop_types.csv" \
                      .format(building=building_name)
    formatted_all_srf_df.to_csv(output_path, index=False, na_rep=0)