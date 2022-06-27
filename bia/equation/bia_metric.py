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
from geopandas import GeoDataFrame as gdf

import cea.utilities.parallel
from cea.constants import HOURS_IN_YEAR
from cea.resources.radiation_daysim import daysim_main, geometry_generator
from cea.utilities.standardize_coordinates import get_lat_lon_projected_shapefile

from bia.equation.bia_crop_cycle import calc_properties_crop_db, calc_chunk_day_crop, calc_crop_cycle


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
    t0 = time.perf_counter()

    # the selected crop type
    type_crop = config.agriculture.type_crop

    # read crop properties in the BIA database for the selected crop type
    crop_properties = calc_properties_crop_db(config)
    print("Gathering the properties of {type_crop}.".format(type_crop=type_crop))

    # read the daily DLI results
    dli_path = config.scenario + "/outputs/data/potentials/agriculture/{building}_DLI_daily.csv"\
        .format(building=building_name)
    cea_dli_results = pd.read_csv(dli_path)


    # activate the function that calculates
    # the yields [kg/year] for the selected crop type on each building envelope surface.
    yield_srf = calc_crop_yields(locator, config, building_name, cea_dli_results)


    # activate the function that calculates
    # the carbon emission reduction (CO2-equivalent) for the selected crop type on each building envelope surface.


    # activate the function that calculates
    # the CAPEX and OPEX for the selected crop type on each building envelope surface.


    # merge the results as a DataFrame


    # write the BIA results
    output_path = config.scenario + \
                  "/outputs/data/potentials/agriculture/{building}_BIA.csv".format(building=building_name)
    sensors_metadata_clean_DLI_daily.to_csv(output_path, index=True,
                                            float_format='%.2f',
                                            na_rep=0)


    print('BIA assessments for each surface on Building', building_name,
          'done - time elapsed: %.2f seconds' % (time.perf_counter() - t0))

# determine the location of the scenario (one of the four zones)
def locate_sce(latitude, longitude):

    """
    This function determines the location of the scenario in the four zones.
    The four zones: North, Tropic of Cancer, Tropic of Cancer, and South

    :param latitude: the latitude of the scenario
    :type latitude: float
    :param longitude: the longitude of the scenario
    :type longitude: float

    :return: zone_sce: one of the four zones on earth
    :type zone_sce: string

    """

    if latitude > 23.5:
        zone_sce = 'zone_north'
    if 0 <= latitude <= 23.5:
        zone_sce = 'zone_cancer'
    if -23.5 <= latitude < 0:
        zone_sce = 'zone_capricorn'
    if latitude < -23.5:
        zone_sce = 'zone_south'
    return zone_sce

# differentiate the cycles as facing and back from the sun for south/north building surfaces
# when they are located in the tropics of cancer and capricorn
def differentiate_cycl_srf_sun(latitude, longitude, date_srf, crop_properties):

    """
    This function determines the location of the scenario in the four zones.
    The four zones: North, Tropic of Cancer, Tropic of Cancer, and South

    :param latitude: the latitude of the scenario
    :type latitude: float
    :param longitude: the longitude of the scenario
    :type longitude: float

    :return: zone_sce: one of the four zones on earth
    :type zone_sce: string

    """

    cycl_i_srf = []
    cycl_o_srf = []

    # unpack the properties of the selected crop type: cycle days
    cycl_i_day = int(crop_properties.get('cycl_i_day'))  # growth cycle in days: initial
    cycl_s_day = int(crop_properties.get('cycl_s_day'))  # growth cycle in days: subsequent
    n_cycl = int(crop_properties.get('n_cycl'))  # number of growth cycles: both initial and subsequent

    # if the scenario is not located in the tropics
    if not -23.5 <= latitude <= 23.5:
        pass        # do nothing if the scenario is outside the tropics

    # if the scenario is located in the tropics
    # differentiate the cycles of each surface by facing and back from the sun
    else:
        # get the two tipping dates (non-leap year) based on the latitude
        summer_solstice = 171
        winter_solstice = 355
        degree_per_day = 0.25543478 # (23.5+23.5)/(355-171) = 0.25543478
        tipping_date_1 = int(winter_solstice - (latitude - (-23.5)) / degree_per_day)
        tipping_date_2 = int(summer_solstice - (23.5 - latitude) / degree_per_day)

        season_srf_1 = []
        season_srf_2 = []

        for surface in range(len(date_srf)):
            date_to_split = date_srf[surface]
            # between the two dates, north surfaces are facing the sun
            list_1 = [x for x in date_to_split if tipping_date_2 < x < tipping_date_1].sort()
            # not-between the two dates, south surface are facing the sun
            list_2 = [x for x in date_to_splitif if not tipping_date_2 < x < tipping_date_1].sort()

            # calculate the number of days in each season between the two tipping dates
            season_crop_1 = calc_chunk_day_crop(list_1)
            len_season_1 = [len(x) for x in season_crop_1]
            season_srf_1.append(len_season_1)

            # calculate the number of days in each season outside the two tipping dates
            season_crop_2 = calc_chunk_day_crop(list_2)
            len_season_2 = [len(x) for x in season_crop_2]

            # when the beginning and the end of each year are connected as a single season
            if list_2[0] == 0 and list_2[-1] == 364:
                # remove the original first and last season: then add the merged season at the end of the list
                merged_season_2 = len_season_2[0] + len_season_2[-1]
                len_season_2.pop(-1)  # remove the original last season
                len_season_2.pop(0)   # remove the original first season
                len_season_2.append(merged_season_2)    # add the newly merged season at the end of the seasons

            season_srf_2.append(len_season_2)

        # Calculate the number of growth cycles for each building surface
        cycl_i_srf = calc_n_cycle_season(cycl_i_day, cycl_s_day, n_cycl, season_srf_1)
        cycl_o_srf = calc_n_cycle_season(cycl_i_day, cycl_s_day, n_cycl, season_srf_2)

    return cycl_i_srf, cycl_o_srf


# Calculate the crop yields (kg) for the selected crop type on each building envelope surface
def calc_crop_yields(locator, config, building_name, crop_properties, cea_dli_results):

    """
    This function calculates crop yield in kg for each building envelope surface.
    At the moment, no pest harm is considered.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator

    :param building_name: list of building names in the case study
    :type building_name: Series

    :return: Dataframe with crop_yields in kg for each building envelope surface.
    """

    # activate the function that calculates
    # the eligible dates, seasons, and cycles for the selected crop type on each building surface
    season_srf, cycl_srf, date_srf = calc_crop_cycle(config, building_name)

    # the total number of cycles of an entire year
    cycl_bld_annual = sum([j for i in cycl_srf for j in i])

    # gather the orientation information [e(ast), w(est), s(outh), n(orth)] for each building surface
    orie_srf = cea_dli_results['orientation'].tolist()      # west, east, north, south
    type_srf = cea_dli_results['TYPE'].tolist()     #walls, windows, roofs
    wall_srf = cea_dli_results['wall_type'].tolist()        #upper, lower, side, non_wall
    area_srf = cea_dli_results['AREA_m2'].tolist()

    # read the T2 crop yield results from the database
    yld_bia_e_g_sqm_cycl = crop_properties.get('yld_bia_e_g_sqm_cycl')  # BIA yield in grams: east
    yld_bia_w_g_sqm_cycl = crop_properties.get('yld_bia_w_g_sqm_cycl')  # BIA yield in grams: west
    yld_bia_f_g_sqm_cycl = crop_properties.get('yld_bia_f_g_sqm_cycl')  # BIA yield in grams: facing the sun (north, south)
    yld_bia_b_g_sqm_cycl = crop_properties.get('yld_bia_b_g_sqm_cycl')  # BIA yield in grams: back from the sun (north, south)


    # number of this building's surfaces
    n_surface = len(orie_srf)

    # get the zone of the scenario
    # One of the four zones: North, Tropic of Cancer, Tropic of Cancer, and South
    zone_geometry_df = gdf.from_file(locator.get_zone_geometry())  # filepath to this scenario's zone.shp
    latitude, longitude = get_lat_lon_projected_shapefile(zone_geometry_df)     # get the lat and lon
    zone_sce = locate_sce(latitude, longitude)   # determine the location of the scenario (one of the four zones)

    # differentiate the cycles as facing and back from the sun for south/north building surfaces
    # when they are located in the tropics of cancer and capricorn
    cycl_i_srf, cycl_o_srf = differentiate_cycl_srf_sun(locator, config, building_name)

    yield_srf = []
    # at least 1 surface is eligible to grow the selected crop type
    if not cycl_bld_annual == 0:

        # annual yield of the selected crops in grams for each building surface
        for surface in range(n_surface):

            cycl = cycl_srf[surface]
            date = date_srf[surface]
            orie = orie_srf[surface]
            area = area_srf[surface]
            cycl_i = cycl_i_srf[surface]
            cycl_o = cycl_o_srf[surface]

            # when the surface is located north to the Tropic of Cancer
            if zone_sce == 'zone_north':
                if orie == 'east':
                    yield_g = yld_bia_e_g_sqm_cycl * area * sum(cycl)

                if orie == 'west':
                    yield_g = yld_bia_w_g_sqm_cycl * area * sum(cycl)

                if orie == 'south':
                    yield_g = yld_bia_f_g_sqm_cycl * area * sum(cycl)

                if orie == 'north':
                    yield_g = yld_bia_b_g_sqm_cycl * area * sum(cycl)

            # when the surface is located between the Tropic of Cancer and the Tropic of Capricorn
            if zone_sce == 'zone_cancer' or 'zone_capricorn':
                if orie == 'east':
                    yield_g = yld_bia_e_g_sqm_cycl * area * sum(cycl)

                if orie == 'west':
                    yield_g = yld_bia_w_g_sqm_cycl * area * sum(cycl)

                if orie == 'south':
                    yield_g_f = yld_bia_f_g_sqm_cycl * area * sum(cycl_o)
                    yield_g_b = yld_bia_b_g_sqm_cycl * area * sum(cycl_i)
                    yield_g = yield_g_f + yield_g_b

                if orie == 'north':
                    yield_g_f = yld_bia_f_g_sqm_cycl * area * sum(cycl_i)
                    yield_g_b = yld_bia_b_g_sqm_cycl * area * sum(cycl_o)
                    yield_g = yield_g_f + yield_g_b

            # when the surface is located south to the Tropic of Capricorn
            if zone_sce == 'zone_south':
                if orie == 'east':
                    yield_g = yld_bia_e_g_sqm_cycl * area * sum(cycl)

                if orie == 'west':
                    yield_g = yld_bia_w_g_sqm_cycl * area * sum(cycl)

                if orie == 'south':
                    yield_g = yld_bia_b_g_sqm_cycl * area * sum(cycl)

                if orie == 'north':
                    yield_g = yld_bia_f_g_sqm_cycl * area * sum(cycl)

            yield_srf.append(yield_g)

    # no surface is eligible to grow the selected crop type
    else:  # No surface meets the minimum DLI requirement of the selected crop type
        print("Unfortunately, {type_crop} is unlikely to grow on the BIA-permissible surfaces on the roof and facade "
              "of Building {building_name}".format(type_crop=type_crop, building_name=building_name))
        pass

    return yield_srf

# calculate ghg emissions for each surface, crop as a produce
def calc_crop_ghg(locator, config, yield_srf, building_name):

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
def calc_crop_opex(locator, config, building_name):
    """
    This function calculates the OPEX in USD (seed, pesticide, and maintenance) for each surface.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator

    :param building_name: list of building names in the case study
    :type building_name: Series

    :return: DataFrame with OPEXes (seed and pesticide) for each building envelope surface.

    """

    return crop_opex

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
