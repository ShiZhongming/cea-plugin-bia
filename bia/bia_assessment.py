"""
This script calculates:
the main script to activate all BIA-related assessments.
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

__author__ = "Zhongming Shi"
__copyright__ = "Copyright 2022, Future Cities Laboratory, Singapore - ETH Zurich; " \
                "University of Calgary, Alberta, Canada"
__credits__ = ["Zhongming Shi"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Zhongming Shi"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"


class BiaAssessmentPlugin(cea.plugin.CeaPlugin):

    pass


# filter by crop on wall/roof/window user-defined in config.file
def filter_crop_srf(locator, config, building_name, bia_metric_srf_df):

    """
    This function filters out the calculate metrics of the surfaces that are not intended by the user to have BIA.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator

    :param building_name: list of building names in the case study
    :type building_name: Series

    :return: DataFrame with BIA metrics for each surfaces intended to have BIA.

    """

    # get the user inputs
    bool_roof = config.agriculture.crop_on_roof
    bool_window = config.agriculture.crop_on_window
    bool_wall_u = config.agriculture.crop_on_wall_under_window
    bool_wall_b = config.agriculture.crop_on_wall_between_window

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
                                             mask_df['mask_wall_b'])]

    # remove the unwanted surfaces
    bia_metric_srf_df_filterred = bia_metric_srf_df[mask_df['mask']]

    return bia_metric_srf_df_filterred

# create aggregated results for each building in one .csv file
def bia_result_aggregate_write(locator, config, building_name):
    """
    This function aggregates the results of each building as a 'total' file and write to disk.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator

    :param building_name: list of building names in the case study
    :type building_name: Series

    :return: write to disk: aggregated BIA metrics for each surfaces intended (defined by users) to have BIA.

    """

    # read the daily DLI results
    dli_path = config.scenario + "/outputs/data/potentials/agriculture/{building}_DLI_daily.csv"\
        .format(building=building_name)
    cea_dli_results = pd.read_csv(dli_path)

    # read the calculated BIA metrics results
    metric_path = config.scenario + "/outputs/data/potentials/agriculture/{building}_BIA_metrics.csv"\
        .format(building=building_name)
    cea_metric_results = pd.read_csv(metric_path)

    # merge the two results
    info_srf_df = cea_dli_results.loc[:, ['AREA_m2', 'total_rad_Whm2', 'TYPE', 'wall_type']]
    bia_metric_srf_df_all = pd.concat([info_srf_df, cea_metric_results], axis=1)

    # filter the ones not wanted by the user
    bia_to_write = filter_crop_srf(locator, config, building_name, bia_metric_srf_df_all)\
        .drop(['TYPE', 'wall_type', 'Unnamed: 0', 'yield_kg_per_sqm_per_year', 'total_rad_Whm2'], axis=1)\
        .sum(axis=0)\
        .to_frame()\
        .T

    # recalculate the yield per square metre
    yield_kg_per_sqm_per_year = bia_to_write['yield_kg_per_year'] / bia_to_write['AREA_m2']

    # aggregate each metric for the entire building
    bia_to_write.insert(loc=0, column='BUILDING', value=building_name)
    bia_to_write.insert(loc=3, column='yield_kg_per_sqm_per_year', value=yield_kg_per_sqm_per_year)

    # write to disk
    bia_path = config.scenario + "/outputs/data/potentials/agriculture/BIA_assessment_total.csv"
    # when the total file has not been created yet, create the file and export the DataFrame
    if not os.path.exists(bia_path):
        bia_to_write.to_csv(bia_path, index=False, float_format='%.2f', na_rep=0)

    # when the total file has already been created, open the file and add one row at the end
    else:
        row = bia_to_write.iloc[0].tolist()     # list of content to append
        row_hundredth = [round(float(num), 2) for num in row[1:]]
        row_hundredth.insert(0, row[0])

        with open(bia_path, 'a') as f_object:

            writer_object = writer(f_object)
            writer_object.writerow(row_hundredth)
            f_object.close()   # Close the file object


def main(config):
    assert os.path.exists(config.scenario), 'Scenario not found: %s' % config.scenario
    locator = cea.inputlocator.InputLocator(config.scenario, config.plugins)

    # BIA assessment
    print('Executing CEA Building-integrated agriculture assessment (BIA) for {type_crop}'
          .format(type_crop=config.agriculture.type_crop))
    building_names = locator.get_zone_building_names()
    num_process = config.get_number_of_processes()
    n = len(building_names)

    # # DLI calculations for each surface of each building
    # cea.utilities.parallel.vectorize(calc_DLI, num_process)\
    #     (repeat(locator, n), repeat(config, n), building_names)
    #
    # # BIA metrics for each surface of each building
    # cea.utilities.parallel.vectorize(calc_bia_metric, num_process)\
    #     (repeat(locator, n), repeat(config, n), building_names)

    # aggregate the results of each building as a 'total' file and write to disk
    # if the file exists, delete it
    bia_path = config.scenario + "/outputs/data/potentials/agriculture/BIA_assessment_total.csv"
    if os.path.exists(bia_path):
        os.remove(bia_path)
    cea.utilities.parallel.vectorize(bia_result_aggregate_write, num_process)\
        (repeat(locator, n), repeat(config, n), building_names)

if __name__ == '__main__':
    main(cea.config.Configuration())
