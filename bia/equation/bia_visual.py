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

    # test if all the required inputs are in place
    # if not all in place, missing BIA profiler's only
    if check_bia_exist(locator, config, building_name) == 1:
        print("Some or all the input files for BIA plotter is/are missing. "
              "Consider to check if BIA Profiler has been executed successfully.")

    # if not all in place, missing BIA assessment's only
    elif check_bia_exist(locator, config, building_name) == 2:
        print("Some or all the input files for BIA plotter is/are missing. "
              "Consider to check if BIA Assessment has been executed successfully.")

    # if not all in place, missing both BIA profiler's and BIA assessment's
    elif check_bia_exist(locator, config, building_name) == 3:
        print("Some or all the input files for BIA plotter is/are missing. "
              "Consider to check if both BIA Profiler and BIA Assessment have been executed successfully.")

    # if all in place
    elif check_bia_exist(locator, config, building_name) == 0:
        # generate each crop type's planting calendar, all crop types combined's planting calendar
        # then write to disk
        visualise_crop_calendar_by_orie_floo(locator, config, building_name)

        # generate all information to be included in the dashboard and write to disk
        visualise_crop_assessment_by_orie_floo(locator, config, building_name)

        # print a message when all the .csv files have been generated
        print('All the .csv files for', building_name,
              'for BIA visualisation have been generated. - time elapsed: %.2f seconds' % (time.perf_counter() - t0))

    # Error, this scenario should not be existing.
    else:
        print('Error: please raise error at GitHub.')


# test if all the required inputs are in place
def check_bia_exist(locator, config, building_name):
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

    :return exist_test: Integer
    0 if all the required files are in place;
    1 if not all the required files are in place (missing BIA profiler's only);
    2 if not all the required files are in place (missing BIA assessment's only);
    3 if not all the required files are in place (missing both BIA profiler's and BIA assessment's)

    :type exist_test: Integer
    """
    # all the crop types to be plotted
    types_crop = config.crop_plot.types_crop

    # the path to the overall planting calendar
    crop_profile_path = config.scenario + \
                  "/outputs/data/potentials/agriculture/{building}_BIA_crop_profile and planting calendar.csv" \
        .format(building=building_name)

    # check through the files created by BIA Profiler/BIA Assessment
    test_list = []
    for type_crop in range(len(types_crop)):
        # the path to the BIA assessment results created by BIA profiler
        crop_assessment_path = config.scenario + \
                               "/outputs/data/potentials/agriculture/surface/{building}_BIA_metrics_{type_crop}.csv" \
                                   .format(building=building_name, type_crop=types_crop[type_crop])
        test_each = os.path.exists(crop_assessment_path)
        test_list.append(test_each)


    # if all the required input files are not in place (missing BIA profiler's only)
    if not os.path.exists(crop_profile_path) and all(test_list):
        exist_test = 1

    # if all the required input files are not in place (missing BIA assessment's only)
    elif os.path.exists(crop_profile_path) and not all(test_list):
        exist_test = 2

    # if all the required input files are not in place (missing both BIA profiler's and BIA assessment's)
    elif not os.path.exists(crop_profile_path) and not all(test_list):
        exist_test = 3

    # if all the required input files are in place
    else:
        exist_test = 0

    return exist_test


# generate each crop type's planting calendar and write to disk
def visualise_crop_calendar_by_orie_floo(locator, config, building_name):

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

    # the path to the overall planting calendar
    crop_profile_path = config.scenario + \
                  "/outputs/data/potentials/agriculture/{building}_BIA_crop_profile and planting calendar.csv" \
        .format(building=building_name)

    # all the crop types to be plotted
    types_crop = config.crop_plot.types_crop

    # read the overall planting calendar
    crop_profile_df = pd.read_csv(crop_profile_path)
    crop_calendar_df = pd.merge(info_srf_df, crop_profile_df, how='left',
                               left_on=['srf_index'],
                               right_on=['srf_index']
                               )

    # create an empty list to store the calendars for each crop type
    calendar_list = []

    # loop for each crop type
    for type_crop in range(len(types_crop)):

        # create a DataFrame to store the outcome
        # all surfaces are not grouped and kept as they are
        srf_all_df = info_srf_df

        # create a DataFrame to store the outcome
        # all surfaces are grouped by floor number and orientation
        flor_df = info_srf_df \
            .groupby(['orientation', 'n_floor']) \
            .agg({"srf_index": "sum"}).reset_index()

        for day in range(365):
            # list of eligible crops for each srf
            crops_list = crop_calendar_df[str(day)].tolist()

            # label the days that are suitable to grow such crop type
            srf_all_df[str(day)] = [1 if types_crop[type_crop] in str(x) else 0 for x in crops_list]

            # create the DataFrame of eligible surfaces
            # grouped by floor number and orientation
            # calculate the number of surfaces the surfaces that are suitable to grow such crop type
            # for the same floor number and orientation
            srf_eligible_df = srf_all_df[srf_all_df[str(day)] == 1] \
                .groupby(['orientation', 'n_floor']) \
                .agg({str(day): "sum"}).reset_index()

            # add the column of eligible surfaces grouped by floor number and orientation
            flor_df = pd.merge(flor_df, srf_eligible_df, how='left',
                               left_on=['orientation', 'n_floor'],
                               right_on=['orientation', 'n_floor']
                               )

            # calculate the number of surfaces
            # for the same floor number and orientation
            flor_df['count'] = srf_all_df \
                .groupby(['orientation', 'n_floor']) \
                .count().reset_index()[str(day)]

            # calculate the percentage of surfaces suitable to grow such crop type
            # for the same floor number and orientation
            flor_df[str(day)] = flor_df[str(day)] / flor_df['count']

            # drop the 'suit' and 'count' columns
            flor_df = flor_df.drop(columns=['count'])

        # append each crop type's calendar into a list
        calendar_list.append(srf_all_df.iloc[:, 3:])

        # process the DataFrame (each crop type) to the needed format
        date = pd.date_range('1/1/2022', periods=365, freq='D').strftime("%Y-%m-%d").tolist()
        handle = ['orientation', 'n_floor', 'srf_index'] + date
        formatted_flor_df = flor_df.T.reset_index(drop=True)
        formatted_flor_df.insert(0, 'date', handle)

        # write the outcome (each crop type) to disk
        output_path = config.scenario + \
                      "/outputs/data/potentials/agriculture/plot/{building}_BIA_planting_calendar_{type_crop}.csv" \
                          .format(building=building_name, type_crop=types_crop[type_crop])
        formatted_flor_df.to_csv(output_path, index=False, na_rep=0)

    # create a combined calendar DataFrame
    srf_all_df = sum(calendar_list)
    srf_all_df = pd.concat([info_srf_df.loc[:, ['srf_index', 'orientation', 'n_floor']], srf_all_df], axis=1)

    # create a DataFrame to store the outcome
    # all surfaces are grouped by floor number and orientation
    flor_all_df_info = info_srf_df \
        .groupby(['orientation', 'n_floor']) \
        .agg({"srf_index": "sum"}).reset_index()

    # calculate the number of surfaces the surfaces that are suitable to grow such crop type
    # for the same floor number and orientation
    flor_all_df_suit = srf_all_df \
        .groupby(['orientation', 'n_floor']) \
        .sum().reset_index()

    # calculate the number of surfaces
    # for the same floor number and orientation
    flor_all_df_count = srf_all_df \
        .groupby(['orientation', 'n_floor']) \
        .count().reset_index()

    # calculate the percentage of surfaces suitable to grow such crop type
    # for the same floor number and orientation
    flor_all_df_day = flor_all_df_suit.iloc[:, 2:].div(flor_all_df_count.iloc[:, 3:])
    flor_all_df = pd.concat([flor_all_df_info, flor_all_df_day], axis=1)

    # process the DataFrame (all crop types) to the needed format
    date = pd.date_range('1/1/2022', periods=365, freq='D').strftime("%Y-%m-%d").tolist()
    handle = ['orientation', 'n_floor', 'srf_index'] + date
    formatted_all_flor_df = flor_all_df.T.reset_index(drop=True)
    formatted_all_flor_df.insert(0, 'date', handle)

    # write the outcome (all crop types) to disk
    output_path = config.scenario + \
                  "/outputs/data/potentials/agriculture/plot/{building}_BIA_planting_calendar_all_crop_types.csv" \
                      .format(building=building_name)
    formatted_all_flor_df.to_csv(output_path, index=False, na_rep=0)


# generate each crop type's BIA assessment and write to disk
def visualise_crop_assessment_by_orie_floo(locator, config, building_name):

    """
     This function generates .csv files to be used as the input for bia visualisation
     for each crop type's BIA assessment
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
    info_srf_df = cea_dli_results.loc[:, ['srf_index', 'orientation', 'n_floor', 'wall_type', 'AREA_m2']]

    # all the crop types to be plotted
    types_crop = config.crop_plot.types_crop

    # loop for each crop type
    for type_crop in range(len(types_crop)):

        # read the BIA assessment outcome created by BIA profiler or BIA assessment
        crop_assessment_path = config.scenario + \
                  "/outputs/data/potentials/agriculture/surface/{building}_BIA_metrics_{type_crop}.csv"\
                      .format(building=building_name, type_crop=types_crop[type_crop])


        # read the BIA Assessment results for the crop type
        crop_assessment_df = pd.read_csv(crop_assessment_path)

        # merge the info DataFrame and BIA assessment DataFrame
        srf_df = pd.merge(info_srf_df, crop_assessment_df, how='left',
                          left_on=['srf_index'],
                          right_on=['srf_index']
                          )

        # create a DataFrame to store the outcome
        # all surfaces are grouped by floor number and orientation
        flor_df = info_srf_df \
            .groupby(['orientation', 'n_floor']) \
            .agg({"srf_index": "sum"}).reset_index()

        #### Calculations

        # surface area in square metre by floor number and orientation
        flor_area = srf_df \
            .groupby(['orientation', 'n_floor']) \
            .agg({"AREA_m2": "sum"}).reset_index()['AREA_m2'].tolist()

        # annual yield in kg by floor number and orientation
        flor_yield = srf_df \
            .groupby(['orientation', 'n_floor']) \
            .agg({"yield_kg_per_year": "sum"}).reset_index()['yield_kg_per_year'].tolist()
        flor_yield_sqm = [x / y for x, y in zip(flor_yield, flor_area)] # per sqm surface area

        # GHG in co2-eq kg by floor number and orientation
        flor_ghg = srf_df \
            .groupby(['orientation', 'n_floor']) \
            .agg({"ghg_kg_co2_bia": "sum"}).reset_index()['ghg_kg_co2_bia'].tolist()
        flor_ghg_sqm = [x / y for x, y in zip(flor_ghg, flor_area)] # per sqm surface area

        # GHG in co2-eq (mys) kg by floor number and orientation
        flor_ghg_mys = srf_df \
            .groupby(['orientation', 'n_floor']) \
            .agg({"ghg_kg_co2_mys": "sum"}).reset_index()['ghg_kg_co2_mys'].tolist()

        # GHG savings
        flor_ghg_saving = [x - y for x, y in zip(flor_ghg_mys, flor_ghg)]
        flor_ghg_saving_sqm = [x / y for x, y in zip(flor_ghg_saving, flor_area)] # per sqm surface area

        # energy use in kWh by floor number and orientation
        flor_energy = srf_df\
            .groupby(['orientation', 'n_floor']) \
            .agg({"energy_kWh_bia": "sum"}).reset_index()['energy_kWh_bia'].tolist()
        flor_energy_sqm = [x / y for x, y in zip(flor_energy, flor_area)] # per sqm surface area

        # energy use in kWh (mys) by floor number and orientation
        flor_energy_mys = srf_df\
            .groupby(['orientation', 'n_floor']) \
            .agg({"energy_kWh_mys": "sum"}).reset_index()['energy_kWh_mys'].tolist()

        # energy savings
        flor_energy_saving = [x - y for x, y in zip(flor_energy_mys, flor_energy)]
        flor_energy_saving_sqm = [x / y for x, y in zip(flor_energy_saving, flor_area)] # per sqm surface area

        # water use in litre by floor number and orientation
        flor_water = srf_df \
            .groupby(['orientation', 'n_floor']) \
            .agg({"water_l_bia": "sum"}).reset_index()['water_l_bia'].tolist()
        flor_water_sqm = [x / y for x, y in zip(flor_water, flor_area)] # per sqm surface area

        # water use in litre (mys) by floor number and orientation
        flor_water_mys = srf_df \
            .groupby(['orientation', 'n_floor']) \
            .agg({"water_l_mys": "sum"}).reset_index()['water_l_mys'].tolist()

        # water savings
        flor_water_saving = [x - y for x, y in zip(flor_water_mys, flor_water)]
        flor_water_saving_sqm = [x / y for x, y in zip(flor_water_saving, flor_area)] # per sqm surface area

        # CAPEX in USD by floor number and orientation
        flor_capex = srf_df \
            .groupby(['orientation', 'n_floor']) \
            .agg({"capex_all_USD": "sum"}).reset_index()['capex_all_USD'].tolist()
        flor_capex_sqm = [x / y for x, y in zip(flor_capex, flor_area)] # per sqm surface area

        # revenue in USD by floor number and orientation
        flor_revenue = srf_df \
            .groupby(['orientation', 'n_floor']) \
            .agg({"opex_sell_USD_per_year": "sum"}).reset_index()['opex_sell_USD_per_year'].tolist()
        flor_revenue_sqm = [x / y for x, y in zip(flor_revenue, flor_area)] # per sqm surface area
        flor_revenue_sqm_cr = [i * -1 for i in flor_revenue_sqm]

        # net annual revenue in USD by floor number and orientation
        flor_revenue_net = srf_df \
            .groupby(['orientation', 'n_floor']) \
            .agg({"opex_all_USD_per_year": "sum"}).reset_index()['opex_all_USD_per_year'].tolist()

        flor_revenue_net_sqm = [x / y for x, y in zip(flor_revenue_net, flor_area)] # per sqm surface area
        flor_revenue_net_sqm_cr = [i * -1 for i in flor_revenue_net_sqm]

        #### compile the results into one DataFrame
        flor_df['yield_kg_per_sqm_surface_area_per_year'] = flor_yield_sqm
        flor_df['ghg_kg_co2_per_sqm_bia'] = flor_ghg_sqm
        flor_df['ghg_kg_co2_per_sqm_saving'] = flor_ghg_saving_sqm
        flor_df['energy_kWh_per_sqm_bia'] = flor_energy_sqm
        flor_df['energy_kWh_per_sqm_saving'] = flor_energy_saving_sqm
        flor_df['water_l_per_sqm_bia'] = flor_water_sqm
        flor_df['water_l_per_sqm_saving'] = flor_water_saving_sqm
        flor_df['capex_per_sqm_USD'] = flor_capex_sqm
        flor_df['revenue_per_sqm_USD_per_year'] = flor_revenue_sqm_cr
        flor_df['net_revenue_per_sqm_USD_per_year'] = flor_revenue_net_sqm_cr

        # process the DataFrame (all crop types) to the needed format
        formatted_flor_df = flor_df

        # write the outcome (all crop types) to disk
        output_path = config.scenario + \
                      "/outputs/data/potentials/agriculture/plot/{building}_BIA_Assessment_{type_crop}.csv"\
                          .format(building=building_name, type_crop=types_crop[type_crop])
        formatted_flor_df.to_csv(output_path, index=False, na_rep=0)