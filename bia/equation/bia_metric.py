"""
This script calculates:
crop yields (kg),
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


# calculate the BIA metrics and write to disk
def calc_bia_metric(locator, config, building_name):

    """
    This function calculates the three categories of BIA metrics and write the results to disk.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator

    :param building_name: list of building names in the case study
    :type building_name: Series

    :return: no return; write the results to disc

    """

    t0 = time.perf_counter()

    # read the daily DLI results
    dli_path = config.scenario + "/outputs/data/potentials/agriculture/{building}_DLI_daily.csv"\
        .format(building=building_name)
    cea_dli_results = pd.read_csv(dli_path)

    # activate the function that calculates
    # the eligible dates, seasons, and cycles for the selected crop type on each building surface
    season_srf, cycl_srf, date_srf, cycl_i_srf, cycl_s_srf = calc_crop_cycle(config, building_name)

    # activate the function that calculates
    # the yields [kg/year] for the selected crop type on each building envelope surface
    # plus yields per square metre [kg/sqm/year] building surface area
    print('Calculating yields (kg) for {type_crop}.'.format(type_crop=config.agriculture.type_crop))
    yield_srf_df = calc_crop_yields(locator, config, building_name, cea_dli_results, cycl_srf, date_srf)

    # activate the function that calculates
    # the environmental impacts (ghg emissions, energy consumption, water consumption)
    # for the selected crop type on each building envelope surface
    print('Calculating GHG emissions (kg CO2-equ), energy consumption (kWh), water consumption (L)) for {type_crop}'
          .format(type_crop=config.agriculture.type_crop))
    env_impacts_srf_df = calc_crop_environmental_impact(locator, config, building_name,
                                                        cea_dli_results, date_srf, yield_srf_df)

    # activate the function that calculates
    # the CAPEX (USD per kg vegetable), OPEX (USD per kg vegetable) in USD
    # for the selected crop type on each building envelope surface
    print('Calculating capital operational costs (USD) for {type_crop}'.format(type_crop=config.agriculture.type_crop))
    costs_srf_df = calc_crop_cost(locator, config, building_name,
                                  cea_dli_results, cycl_srf, cycl_i_srf, yield_srf_df, env_impacts_srf_df)

    # merge the results as a DataFrame
    bia_metric_srf_df = pd.concat([yield_srf_df, env_impacts_srf_df, costs_srf_df], axis=1)

    # write the BIA results (all non-filtered surface)
    output_path = config.scenario + \
                  "/outputs/data/potentials/agriculture/{building}_BIA_metrics.csv".format(building=building_name)
    bia_metric_srf_df.to_csv(output_path, index=True, float_format='%.2f', na_rep=0)

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
    elif 0 <= latitude <= 23.5:
        zone_sce = 'zone_cancer'
    elif -23.5 <= latitude < 0:
        zone_sce = 'zone_capricorn'
    else:
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
    :param date_srf: the days (0 to 364, in total 365 days in a non-leap year) that are eligible for growing the
    selected crop type
    :type date_srf: list
    :return: dict with properties of the selected crop type retrieved form the database
    :type dict

    :return: cycl_i_srf: when building located in the tropics, crop cycles between the two tipping dates
    :type cycl_i_srf: list
    :return: cycl_o_srf: when building located in the tropics, crop cycles before or after the two tipping dates
    :type cycl_o_srf: list

    """

    cycl_i_srf = []     # inside the two tipping days
    cycl_o_srf = []     # outside the two tipping days

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
        degree_per_day = 0.25543478     # (23.5+23.5)/(355-171) = 0.25543478
        tipping_date_1 = int(winter_solstice - (latitude - (-23.5)) / degree_per_day)
        tipping_date_2 = int(summer_solstice - (23.5 - latitude) / degree_per_day)

        season_srf_1 = []
        season_srf_2 = []

        for surface in range(len(date_srf)):
            date_to_split = date_srf[surface]

            # between the two dates, north surfaces are facing the sun
            list_1 = [x for x in date_to_split if tipping_date_2 < x < tipping_date_1]
            # not-between the two dates, south surface are facing the sun
            list_2 = [x for x in date_to_split if not tipping_date_2 < x < tipping_date_1]

            # true if between the two dates (north surfaces are facing the sun)
            # false if not-between the two dates (south surfaces are facing the sun)
            # mask = [True for x in date_to_split if tipping_date_2 < x < tipping_date_1].sort()

            # calculate the number of days in each season between the two tipping dates
            # if such number of days does not equal to zero
            if list_1:
                season_crop_1 = calc_chunk_day_crop(list_1)
                len_season_1 = [len(x) for x in season_crop_1]
                season_srf_1.append(len_season_1)

            else:
                season_srf_1.append([0])

            # calculate the number of days in each season outside the two tipping dates
            # if such number of days does not equal to zero
            if list_2:
                season_crop_2 = calc_chunk_day_crop(list_2)
                len_season_2 = [len(x) for x in season_crop_2]

                # when the beginning and the end of each year are connected as a single season
                if list_2[0] == 0 and list_2[-1] == 364:
                    # remove the original first and last season: then add the merged season at the end of the list
                    merged_season_2 = len_season_2[0] + len_season_2[-1]
                    len_season_2.pop(-1)  # remove the original last season
                    len_season_2.pop(0)  # remove the original first season
                    len_season_2.append(merged_season_2)  # add the newly merged season at the end of the seasons

                season_srf_2.append(len_season_2)

            else:
                season_srf_2.append([0])

        # Calculate the number of growth cycles for each building surface
        cycl_i_srf, cycl_i_i_srf, cycl_i_s_srf = calc_n_cycle_season(cycl_i_day, cycl_s_day, n_cycl, season_srf_1)
        cycl_o_srf, cycl_o_i_srf, cycl_o_s_srf = calc_n_cycle_season(cycl_i_day, cycl_s_day, n_cycl, season_srf_2)

    return cycl_i_srf, cycl_o_srf


# differentiate the eligible dates as in wet and dry seasons
def day_counts_srf_wet_dry_sgp(latitude, longitude, date_srf):

    """
    This function determines the location of the scenario in the four zones.
    The four zones: North, Tropic of Cancer, Tropic of Cancer, and South

    :param latitude: the latitude of the scenario
    :type latitude: float
    :param longitude: the longitude of the scenario
    :type longitude: float
    :param date_srf: the days (0 to 364, in total 365 days in a non-leap year) that are eligible for growing the
    selected crop type
    :type date_srf: list

    :return: n_day_wet_srf: number of days in the wet season for each building surface
    :type n_day_wet_srf: list
    :return: n_day_dry_srf: number of days in the dry season for each building surface
    :type n_day_dry_srf: list

    """
    n_day_wet_srf = []
    n_day_dry_srf = []

    # if the scenario is not located in the tropics
    if not ((-5 <= latitude <= 5) and (100 <= longitude <= 105)):
        print('The building is located far away from Singapore. '
              'No data available for the water usage in this location.')
    # do nothing if the scenario is not near Singapore

    # if the scenario is located near Singapore
    # count the days of wet/dry periods for each building surface
    else:
        # the four tipping dates (non-leap year) for changing the seasons (wet/dry)
        # reference: meteorological service singapore, SGP gov http://www.weather.gov.sg/climate-climate-of-singapore/
        print('The building is located near Singapore. Hence, it applies the same dry/wet seasons.')

        date_d_w_dec = 334  # december 1
        date_w_d_jan = 14   # January 15
        date_d_w_may = 150  # May 31
        date_w_d_sep = 272  # September 30

        n_day_wet_srf = []
        n_day_dry_srf = []

        for surface in range(len(date_srf)):
            date_to_split = date_srf[surface]
            # January 1 to January 15, wet
            list_1 = [x for x in date_to_split if 0 < x < date_w_d_jan]
            # January 16 to May 31, dry
            list_2 = [x for x in date_to_split if date_w_d_jan <= x < date_d_w_may]
            # June 1 to September 30, wet
            list_3 = [x for x in date_to_split if date_d_w_may <= x < date_w_d_sep]
            # October 1 to November 30, dry
            list_4 = [x for x in date_to_split if date_w_d_sep <= x < date_d_w_dec]
            # December 1 to December 31, wet
            list_5 = [x for x in date_to_split if date_d_w_dec <= x < 364]

            # calculate the number of days in wet and dry periods
            n_day_wet = len(list_1) + len(list_3) + len(list_5)
            n_day_dry = len(list_2) + len(list_4)

            # record the two numbers for each surface
            n_day_wet_srf.append(n_day_wet)
            n_day_dry_srf.append(n_day_dry)

    return n_day_wet_srf, n_day_dry_srf


# Calculate the crop yields (kg) for the selected crop type on each building envelope surface
def calc_crop_yields(locator, config, building_name, cea_dli_results, cycl_srf, date_srf):

    """
    This function calculates crop yield in kg for each building envelope surface.
    At the moment, no pest harm is considered.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator
    :param building_name: list of building names in the case study
    :type building_name: Series
    :param cea_dli_results: dli results stored in the csv file via
    "/outputs/data/potentials/agriculture/{building}_DLI_daily.csv"
    :type cea_dli_results: DataFrame
    :param cycl_srf: number of cycles, including both initial and subsequent ones,
    for each building surface of a whole year
    :type cycl_srf: list
    :param date_srf: the days (0 to 364, in total 365 days in a non-leap year) that are eligible for growing the
    selected crop type
    :type date_srf: list

    :return: yield_srf_df: each building surface's annual yield (kg/year)
    and yield per square metre surface area (kg/sqm/year)
    :type yield_srf_df: DataFrame
    """

    # the selected crop type
    type_crop = config.agriculture.type_crop

    # read crop properties in the BIA database for the selected crop type
    crop_properties = calc_properties_crop_db(config)
    print("Gathering the properties of {type_crop}.".format(type_crop=type_crop))

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

    # get the zone of the scenario
    # One of the four zones: North, Tropic of Cancer, Tropic of Cancer, and South
    zone_geometry_df = gdf.from_file(locator.get_zone_geometry())  # filepath to this scenario's zone.shp
    latitude, longitude = get_lat_lon_projected_shapefile(zone_geometry_df)     # get the lat and lon
    zone_sce = locate_sce(latitude, longitude)   # determine the location of the scenario (one of the four zones)
    print('The scenario is located in {zone_sce}.'.format(zone_sce=zone_sce))

    # differentiate the cycles as facing and back from the sun for south/north building surfaces
    # when they are located in the tropics of cancer and capricorn
    cycl_i_srf, cycl_o_srf = differentiate_cycl_srf_sun(latitude, longitude, date_srf, crop_properties)

    yield_srf = []
    yield_srf_per_sqm = []
    # at least 1 surface is eligible to grow the selected crop type

    if not cycl_bld_annual == 0:

        # annual yield of the selected crops in grams for each building surface
        for surface in range(len(orie_srf)):

            cycl = cycl_srf[surface]
            orie = orie_srf[surface]
            area = area_srf[surface]

            # when the surface is located north to the Tropic of Cancer
            if zone_sce == 'zone_north':
                if orie == 'east':
                    yield_g = yld_bia_e_g_sqm_cycl * area * sum(cycl)

                elif orie == 'west':
                    yield_g = yld_bia_w_g_sqm_cycl * area * sum(cycl)

                elif orie == 'south':
                    yield_g = yld_bia_f_g_sqm_cycl * area * sum(cycl)

                else:
                    yield_g = yld_bia_b_g_sqm_cycl * area * sum(cycl)

            # when the surface is located south to the Tropic of Capricorn
            elif zone_sce == 'zone_south':
                if orie == 'east':
                    yield_g = yld_bia_e_g_sqm_cycl * area * sum(cycl)

                elif orie == 'west':
                    yield_g = yld_bia_w_g_sqm_cycl * area * sum(cycl)

                elif orie == 'south':
                    yield_g = yld_bia_b_g_sqm_cycl * area * sum(cycl)

                else:
                    yield_g = yld_bia_f_g_sqm_cycl * area * sum(cycl)

            # when the surface is located between the Tropic of Cancer and the Tropic of Capricorn
            else:
                cycl_i = cycl_i_srf[surface]
                cycl_o = cycl_o_srf[surface]

                if orie == 'east':
                    yield_g = yld_bia_e_g_sqm_cycl * area * sum(cycl)

                elif orie == 'west':
                    yield_g = yld_bia_w_g_sqm_cycl * area * sum(cycl)

                elif orie == 'south':
                    yield_g_f = yld_bia_f_g_sqm_cycl * area * sum(cycl_o)
                    yield_g_b = yld_bia_b_g_sqm_cycl * area * sum(cycl_i)
                    yield_g = yield_g_f + yield_g_b

                else:
                    yield_g_f = yld_bia_f_g_sqm_cycl * area * sum(cycl_i)
                    yield_g_b = yld_bia_b_g_sqm_cycl * area * sum(cycl_o)
                    yield_g = yield_g_f + yield_g_b

            yield_srf.append(yield_g/1000)  # record the results and convert to kilograms
            yield_srf_per_sqm.append(yield_g/1000/area)

        # merged the results for each surface in DataFrame
        data = np.array_split(yield_srf + yield_srf_per_sqm, 2)

        yield_srf_df = pd.DataFrame(data).T.fillna(0)
        column_name = ['yield_kg_per_year', 'yield_kg_per_sqm_per_year']
        yield_srf_df.columns = column_name

    # no surface is eligible to grow the selected crop type
    else:  # No surface meets the minimum DLI requirement of the selected crop type
        print("Unfortunately, {type_crop} is unlikely to grow on any of the surfaces "
              "of Building {building_name}".format(type_crop=type_crop, building_name=building_name))
        pass

    return yield_srf_df

# calculate the three environmental impacts
def calc_crop_environmental_impact(locator, config, building_name, cea_dli_results, date_srf, yield_srf_df):

    """
    This function calculates the GHG Emissions (kg CO2-eq), energy (kWh) and water use (litre)
    for each building surface.
    Three categories stages related to 'farm-to-table' are considered: production, processing and transport.
    BIA saves this.

    See Page 11 of the Deloitte report.
    See Page 12 for details on GHG Emissions.

    This function also calculates the energy consumption (kWh per kg) in the three stages.
    BIA save this.

    This function also calculates the water consumption (litres per kg) in the three stages.
    BIA wastes this in Singapore.


    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator
    :param building_name: list of building names in the case study
    :type building_name: Series
    :param cea_dli_results: dli results stored in the csv file via
    "/outputs/data/potentials/agriculture/{building}_DLI_daily.csv"
    :type cea_dli_results: DataFrame
    :param date_srf: the days (0 to 364, in total 365 days in a non-leap year) that are eligible for growing the
    selected crop type
    :type date_srf: list
    :param: yield_srf_df: each building surface's annual yield (kg/year)
    and yield per square metre surface area (kg/sqm/year)
    :type yield_srf_df: DataFrame

    :return env_impacts_srf_df: DataFrame with the three environmental impacts for
    the 6 Deloitte scenarios + BIA scenario.
    :type env_impacts_srf_df: DataFrame

    """

    # the selected crop type
    type_crop = config.agriculture.type_crop

    # read crop properties in the BIA database for the selected crop type
    env_properties = calc_properties_env_db(config)

    # gather the orientation information [e(ast), w(est), s(outh), n(orth)] for each building surface
    orie_srf = cea_dli_results['orientation'].tolist()  # west, east, north, south
    area_srf = cea_dli_results['AREA_m2'].tolist()

    # get the list of scenario names (mys, idn, sgp X 4)
    crop_grow_sce = env_properties['scenario'].tolist()

    # the three metrics for Deloitte/TEMASEK/A-star report
    ghg_srf_sce_all = []
    energy_srf_sce_all = []
    water_srf_sce_all = []

    # get the yield results
    yield_srf = yield_srf_df['yield_kg_per_year'].tolist()

    # the six baseline scenarios
    for sce in crop_grow_sce:
        ghg_kg_co2_per_kg_sce = env_properties[env_properties['scenario'] == sce]['ghg_total'].values[0]
        energy_kWh_per_kg_sce = env_properties[env_properties['scenario'] == sce]['energy_total'].values[0]
        water_litres_per_kg_sce = env_properties[env_properties['scenario'] == sce]['water_total'].values[0]

        # GHG-related calculations
        ghg_srf_sce = [srf * ghg_kg_co2_per_kg_sce for srf in yield_srf]     #scenarios in the database

        # energy-related calculations
        energy_srf_sce = [srf * energy_kWh_per_kg_sce for srf in yield_srf]      #scenarios in the database

        # water-related calculations
        water_srf_sce = [srf * water_litres_per_kg_sce for srf in yield_srf]     #scenarios in the database

        # record the results
        ghg_srf_sce_all.append(ghg_srf_sce)
        energy_srf_sce_all.append(energy_srf_sce)
        water_srf_sce_all.append(water_srf_sce)

    # the three metrics for BIA scenarios
    # GHG-related calculations
    ghg_kg_co2_per_kg_bia = 0.271572    # assumption: we use the Singapore soil-cultivated non-greenhouse's
                                        # GHG emissions in the production stage to represent BIA's GHG emissions
                                        # Deloitte/TEMASEK/A-star report: source
    ghg_srf_bia = [srf * ghg_kg_co2_per_kg_bia for srf in yield_srf]

    # energy-related calculations
    energy_kWh_per_kg_bia = 0  # assumption: no energy (electricity) required in BIA
    energy_srf_bia = [srf * energy_kWh_per_kg_bia for srf in yield_srf]

    # water-related calculations
    water_litres_per_sqm_w = 0.1  # based on T2 lab experiments, wet season
    water_litres_per_sqm_d = 0.2  # based on T2 lab experiments, dry season
    water_srf_bia = []

    # differentiate the eligible dates as in wet and dry seasons, if in/near Singapore
    zone_geometry_df = gdf.from_file(locator.get_zone_geometry())  # filepath to this scenario's zone.shp
    latitude, longitude = get_lat_lon_projected_shapefile(zone_geometry_df)     # get the lat and lon
    n_day_wet_srf, n_day_dry_srf = day_counts_srf_wet_dry_sgp(latitude, longitude, date_srf)


    for surface in range(len(orie_srf)):
        # if near Singapore and we know the dry and wet seasons
        if n_day_dry_srf and n_day_wet_srf:
            orie = orie_srf[surface]
            area = area_srf[surface]
            n_day_wet = n_day_wet_srf[surface]
            n_day_dry = n_day_dry_srf[surface]

            if orie == 'east':
                water_litres = water_litres_per_sqm_w * (n_day_wet + n_day_dry) * area

            elif orie == 'west':
                water_litres = 0

            else:
                water_litres = water_litres_per_sqm_d * n_day_dry * area

            water_srf_bia.append(water_litres)  # record the results and convert to kilograms

        else:
            water_srf_bia.append(0)

    # merged the (6+1) scenarios for each surface in DataFrame
    data = ghg_srf_sce_all + [ghg_srf_bia] + energy_srf_sce_all + [energy_srf_bia] + water_srf_sce_all + [water_srf_bia]
    env_impacts_srf_df = pd.DataFrame(data).T.fillna(0)

    # create the column names for the DataFrame
    crop_grow_sce.append('bia')
    column_name = []
    for metrics in ['ghg_kg_co2', 'energy_kWh', 'water_l']:
        head = []
        for sce in crop_grow_sce:
            name = metrics + '_' + sce
            head.append(name)
        column_name.append(head)

    column_name = [x for xs in column_name for x in xs]
    env_impacts_srf_df.columns = column_name

    return env_impacts_srf_df


# calculate the costs (CAPEX, OPEX, market price of the yields) for each surface
def calc_crop_cost(locator, config, building_name,
                   cea_dli_results, cycl_srf, cycl_i_srf, yield_srf_df, env_impacts_srf_df):
    """
    This function calculates the CAPEX in USD (infrastructure) for each surface.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator
    :param building_name: list of building names in the case study
    :type building_name: Series
    :param cea_dli_results: dli results stored in the csv file via
    "/outputs/data/potentials/agriculture/{building}_DLI_daily.csv"
    :type cea_dli_results: DataFrame
    :param cycl_srf: number of cycles, including both initial and subsequent ones,
    for each building surface of a whole year
    :type cycl_srf: list
    :param cycl_i_srf:  number of cycles, including initial ones only for each building surface of a whole year
    :type cycl_i_srf: list
    :param: yield_srf_df: each building surface's annual yield (kg/year)
    and yield per square metre surface area (kg/sqm/year)
    :type yield_srf_df: DataFrame
    :param env_impacts_srf_df: DataFrame with the three environmental impacts for
    the 6 Deloitte scenarios + BIA scenario.
    :type env_impacts_srf_df: DataFrame

    :return: DataFrame with CAPEX (infrastructure + soil)
    and OPEX (seed, pesticide, fertilizer, water, -market price) for each building envelope surface.
    :type env_impacts_srf_df: DataFrame

    """

    # the selected crop type
    type_crop = config.agriculture.type_crop

    # gather the orientation information for each building surface
    area_srf = cea_dli_results['AREA_m2'].tolist()

    # sgd to usd
    ex_sgd_usd = 1.40       # July 9, 2022

    # water price in Singapore, assumption monthly usage above 40 cubic metres
    water_price_USD_per_l = (3.69 / 1000) / ex_sgd_usd  # https://www.pub.gov.sg/watersupply/waterprice

    # read cost properties in the BIA database for the selected crop type
    cost_properties = calc_properties_cost_db(config)
    Inv_IR_perc = cost_properties.get('IR_%').values[0]      # interest rate
    Inv_LT = cost_properties.get('LT_yr').values[0]      # lifetime in years
    shelf_USD_per_sqm = cost_properties.get(
        'shelf_sgd_sqm').values[0] / ex_sgd_usd     # infrastructure cost per square meter surface area in USD/sqm
    soil_USD_per_sqm = cost_properties.get(
        'soil_sgd_sqm').values[0] / ex_sgd_usd     # soil cost per square meter surface area in USD/sqm
    seed_USD_per_sqm = cost_properties.get(
        'seed_sgd_sqm').values[0] / ex_sgd_usd     # seed cost per square meter surface area in USD/sqm
    pesticide_USD_per_sqm = cost_properties.get(
        'pesticide_sgd_sqm').values[0] / ex_sgd_usd    # pesticide cost per square meter surface area in USD/sqm
    fertilizer_USD_per_sqm = cost_properties.get(
        'fertilizer_sgd_sqm').values[0] / ex_sgd_usd   # fertilizer cost per square meter surface area in USD/sqm
    market_price_USD_per_kg = cost_properties.get(
        'mkt_sg_sgd_kg').values[0] / ex_sgd_usd   # fertilizer cost per square meter surface area in USD/sqm

    # BIA's CAPEX and OPEX for each building surface
    capex_infrastructure_srf = []
    capex_soil_srf = []
    capex_all_srf = []
    capex_all_a_srf = []
    opex_seed_a_srf = []
    opex_pesticide_a_srf = []
    opex_fertilizer_a_srf = []
    opex_water_a_srf = []
    opex_sell_a_srf = []
    opex_all_a_srf = []

    for surface in range(len(area_srf)):
        cycl = sum(cycl_srf[surface])        # number of total cycles annually
        cycl_i = sum(cycl_i_srf[surface])     # number of initial cycles
        area = area_srf[surface]        # area of the surface in sqm
        yield_bia = yield_srf_df['yield_kg_per_year'].tolist()[surface]
        water_bia = env_impacts_srf_df['water_l_bia'].tolist()[surface]

        capex_infrastructure = shelf_USD_per_sqm * area
        capex_soil = soil_USD_per_sqm * area
        capex_all = capex_infrastructure + capex_soil   # initial capital cost in USD
        capex_all_a = calc_capex_annualized(capex_all, Inv_IR_perc, Inv_LT)     # annualised capital cost in USD/year

        opex_seed_a = seed_USD_per_sqm * area * cycl_i
        opex_pesticide_a = pesticide_USD_per_sqm * area * cycl
        opex_fertilizer_a = fertilizer_USD_per_sqm * area * cycl
        opex_water_a = water_price_USD_per_l * water_bia
        opex_sell_a = -1 * market_price_USD_per_kg * yield_bia
        opex_all_a = opex_seed_a + opex_pesticide_a + opex_fertilizer_a + \
                     opex_water_a + opex_sell_a  # annual operational cost in USD/year

        capex_infrastructure_srf.append(capex_infrastructure)
        capex_soil_srf.append(capex_soil)
        capex_all_srf.append(capex_all)
        capex_all_a_srf.append(capex_all_a)
        opex_seed_a_srf.append(opex_seed_a)
        opex_pesticide_a_srf.append(opex_pesticide_a)
        opex_fertilizer_a_srf.append(opex_fertilizer_a)
        opex_water_a_srf.append(opex_water_a)
        opex_sell_a_srf.append(opex_sell_a)
        opex_all_a_srf.append(opex_all_a)

    # merged the results for each surface in DataFrame
    data = np.array_split(capex_infrastructure_srf + capex_soil_srf +
                          capex_all_a_srf + capex_all_a_srf +
                          opex_seed_a_srf + opex_pesticide_a_srf + opex_fertilizer_a_srf + opex_water_a_srf +
                          opex_sell_a_srf +
                          opex_all_a_srf,
                          10)
    costs_srf_df = pd.DataFrame(data).T.fillna(0)
    column_name = ['capex_infrastructure_USD', 'capex_soil_USD',
                   'capex_all_USD', 'capex_all_annualised_USD',
                   'opex_seed_USD_per_year', 'opex_pesticide_USD_per_year',
                   'opex_fertilizer_USD_per_year', 'opex_water_USD_per_year',
                   'opex_sell_USD_per_year',
                   'opex_all_USD_per_year'
                   ]
    costs_srf_df.columns = column_name

    return costs_srf_df




