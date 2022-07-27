"""
This script creates:
the crop selection profile for each building surface, based on one of the following user-defined :objectives:

environmental impacts including GHG Emissions (kg CO2-eq), energy (kWh) and water use (litre),
costs including capital and oeprational expenditures (USD)
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

from bia.equation.bia_metric import calc_bia_metric, bia_result_aggregate_write, filter_crop_srf
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
def calc_bia_crop_select(locator, config, building_name):

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

    l_type_crop = config.crop_profile.types_crop
    building_names = locator.get_zone_building_names()
    n = len(building_names)

    # filter out the unwanted surfaces, aggregate the results (one .csv per scenario of crop type) and write to disk
    bia_result_aggregate_write(locator, config, building_name, type_crop)


print('BIA assessments for each surface on Building', building_name,
      'done - time elapsed: %.2f seconds' % (time.perf_counter() - t0))


