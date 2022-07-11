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
from cea.analysis.costs.equations import calc_capex_annualized, calc_opex_annualized

from bia.equation.bia_crop_cycle import calc_properties_crop_db, calc_chunk_day_crop, \
    calc_crop_cycle, calc_properties_env_db, calc_properties_cost_db


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
    print('Calculating yields (kg) for {type_crop}'.format(type_crop=config.agriculture.type_crop))
    yield_srf, yield_srf_per_sqm = calc_crop_yields(locator, config, building_name, cea_dli_results, cycl_srf)

    # activate the function that calculates
    # the environmental impacts (ghg emissions, energy consumption, water consumption)
    # for the selected crop type on each building envelope surface
    print('Calculating GHG emissions (kg CO2-equ), energy consumption (kWh), water consumption (L)) for {type_crop}'
          .format(type_crop=config.agriculture.type_crop))
    env_impacts_srf_df = calc_crop_environmental_impact(locator, config, building_name,
                                                        cea_dli_results, date_srf, yield_srf)

    # activate the function that calculates
    # the CAPEX (USD per kg vegetable), OPEX (USD per kg vegetable) in USD
    # for the selected crop type on each building envelope surface
    print('Calculating capital operational costs (USD) for {type_crop}'.format(type_crop=config.agriculture.type_crop))
    capex_srf, capex_a_srf, opex_srf = calc_crop_cost(locator, config, building_name,
                                                      cea_dli_results, cycl_srf, cycl_i_srf, yield_srf)

    # merge the results as a DataFrame
    bia_metric_srf_df = env_impacts_srf_df
    bia_metric_srf_df['yield_kg'] = yield_srf
    bia_metric_srf_df['yield_kg_per_sqm'] = yield_srf_per_sqm
    bia_metric_srf_df['capex_usd'] = capex_srf
    bia_metric_srf_df['capex_a_usd'] = capex_a_srf
    bia_metric_srf_df['opex_usd'] = opex_srf

    # write the BIA results (all non-filtered surface)
    output_path = config.scenario + \
                  "/outputs/data/potentials/agriculture/{building}_BIA_metrics.csv".format(building=building_name)
    bia_metric_srf_df_filtered.to_csv(output_path, index=True,
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


# differentiate the eligible dates as in wet and dry seasons
def day_counts_srf_wet_dry_sgp(latitude, longitude, date_srf):

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
    n_day_wet_srf = []
    n_day_dry_srf = []

    # if the scenario is not located in the tropics
    if not -5 <= latitude <= 5 and 100 <= longitude <= 105:
        pass        # do nothing if the scenario is not near Singapore

    # if the scenario is located near Singapore
    # count the days of wet/dry periods for each building surface
    else:
        # the four tipping dates (non-leap year) for changing the seasons (wet/dry)
        # reference: meteorological service singapore, SGP gov http://www.weather.gov.sg/climate-climate-of-singapore/

        date_d_w_dec = 334  # december 1
        date_w_d_jan = 14   # January 15
        date_d_w_may = 150  # May 31
        date_w_d_sep = 272  # September 30

        n_day_wet_srf = []
        n_day_dry_srf = []

        for surface in range(len(date_srf)):
            date_to_split = date_srf[surface]
            # January 1 to January 15, wet
            list_1 = [x for x in date_to_split if 0 < x < date_w_d_jan].sort()
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
def calc_crop_yields(locator, config, building_name, cea_dli_results, cycl_srf):

    """
    This function calculates crop yield in kg for each building envelope surface.
    At the moment, no pest harm is considered.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator

    :param building_name: list of building names in the case study
    :type building_name: Series

    :return: a list with crop_yields in kilograms for each building envelope surface.
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

    # differentiate the cycles as facing and back from the sun for south/north building surfaces
    # when they are located in the tropics of cancer and capricorn
    cycl_i_srf, cycl_o_srf = differentiate_cycl_srf_sun(locator, config, building_name)

    yield_srf = []
    yield_srf_per_sqm = []
    # at least 1 surface is eligible to grow the selected crop type
    if not cycl_bld_annual == 0:

        # annual yield of the selected crops in grams for each building surface
        for surface in range(len(orie_srf)):

            cycl = cycl_srf[surface]
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

            yield_srf.append(yield_g/1000)  # record the results and convert to kilograms

    # no surface is eligible to grow the selected crop type
    else:  # No surface meets the minimum DLI requirement of the selected crop type
        print("Unfortunately, {type_crop} is unlikely to grow on the BIA-permissible surfaces on the roof and facade "
              "of Building {building_name}".format(type_crop=type_crop, building_name=building_name))
        pass

    return yield_srf


# calculate ghg emissions for each surface, crop as a produce
def calc_crop_environmental_impact(locator, config, building_name, cea_dli_results, date_srf, yield_srf):

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

    :return: DataFrame with the three environmental impacts for the 6 Deloitte scenarios + BIA scenario.

    """

    # the selected crop type
    type_crop = config.agriculture.type_crop

    # read crop properties in the BIA database for the selected crop type
    env_properties = calc_properties_env_db(config)
    print("Gathering the environmental impacts data of {type_crop}.".format(type_crop=type_crop))

    # gather the orientation information [e(ast), w(est), s(outh), n(orth)] for each building surface
    orie_srf = cea_dli_results['orientation'].tolist()  # west, east, north, south
    area_srf = cea_dli_results['AREA_m2'].tolist()

    # get the list of scenario names (mys, idn, sgp X 4)
    crop_grow_sce = env_properties['scenario'].tolist()

    # the three metrics for Deloitte/TEMASEK/A-star report
    for sce in crop_grow_sce:
        ghg_kg_co2_per_kg_sce = env_properties[env_properties['scenario'] == sce]['ghg_total']
        energy_kWh_per_kg_sce = env_properties[env_properties['scenario'] == sce]['energy_total']
        water_litres_per_kg_sce = env_properties[env_properties['scenario'] == sce]['water_total']

        ghg_srf_sce_all = []
        energy_srf_sce_all = []
        water_srf_sce_all = []

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
    water_litres_per_sqm_w = 0.2  # based on T2 lab experiments, wet season
    water_litres_per_sqm_d = 0.1  # based on T2 lab experiments, dry season
    water_srf_bia = []

    # differentiate the eligible dates as in wet and dry seasons
    n_day_wet_srf, n_day_dry_srf = day_counts_srf_wet_dry_sgp(latitude, longitude, date_srf)

    for surface in range(len(orie_srf)):
        orie = orie_srf[surface]
        area = area_srf[surface]

        if orie == 'east':
            water_litres = water_litres_per_sqm_w * (n_day_wet_srf + n_day_dry_srf) * area

        if orie == 'west':
            water_litres = 0

        if orie == 'south' or 'north':
            water_litres = n_day_dry_srf * water_litres_per_sqm_d * area

        water_srf_bia.append(water_litres)  # record the results and convert to kilograms

    # merged the (6+1) scenarios for each surface in DataFrame
    data = list(chuncks(ghg_srf_sce_all + energy_srf_sce_all + water_srf_sce_all +
                        ghg_srf_bia + energy_srf_bia + water_srf_bia),
                len(orie_srf))

    env_impacts_srf_df = pd.DataFrame(data)

    # create the column names for the DataFrame
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
def calc_crop_cost(locator, config, building_name, cea_dli_results, cycl_srf, cycl_i_srf, yield_srf):
    """
    This function calculates the CAPEX in USD (infrastructure) for each surface.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator

    :param building_name: list of building names in the case study
    :type building_name: Series

    :return: DataFrame with CAPEX (infrastructure) for each building envelope surface.

    """

    # the selected crop type
    type_crop = config.agriculture.type_crop

    # gather the orientation information for each building surface
    area_srf = cea_dli_results['AREA_m2'].tolist()

    # sgd to usd
    ex_sgd_usd = 1.40       # July 9, 2022

    # read cost properties in the BIA database for the selected crop type
    cost_properties = calc_properties_cost_db(config)
    print("Gathering the expenditure information for {type_crop}.".format(type_crop=type_crop))
    Inv_IR_perc = cost_properties.get('IR_%')      # interest rate
    Inv_LT = cost_properties.get('LT_yr')       # lifetime in years
    shelf_USD_per_sqm = cost_properties.get(
        'shelf_sgd_sqm') / ex_sgd_usd     # initial investment (shelf + soil) per square meter surface area in USD/sqm
    soil_USD_per_sqm = cost_properties.get(
        'soil_sgd_sqm') / ex_sgd_usd     # initial investment (shelf + soil) per square meter surface area in USD/sqm
    seed_USD_per_sqm = cost_properties.get(
        'seed_sgd_sqm') / ex_sgd_usd     # seed cost per square meter surface area in USD/sqm
    pesticide_USD_per_sqm = cost_properties.get(
        'pesticide_sgd_sqm') / ex_sgd_usd    # pesticide cost per square meter surface area in USD/sqm
    fertilizer_USD_per_sqm = cost_properties.get(
        'fertilizer_sgd_sqm') / ex_sgd_usd   # fertilizer cost per square meter surface area in USD/sqm
    market_price_USD_per_kg = cost_properties.get(
        'mkt_sg_sgd_kg') / ex_sgd_usd   # fertilizer cost per square meter surface area in USD/sqm

    # CAPEX and OPEX of the selected crops in grams for each building surface
    for surface in range(len(area_srf)):
        cycl = cycl_srf[surface]        # number of total cycles annually
        cycl_i = cycl_i_srf[surface]      # number of initial cycles
        area = area_srf[surface]        # area of the surface in sqm

        capex = shelf_USD_per_sqm * area + soil_USD_per_sqm * area   # initial capital expenditure in USD
        capex_a = calc_capex_annualized(capex, Inv_IR_perc, Inv_LT)       # annualised capital expenditure in USD/year

        opex = seed_USD_per_sqm * area * cycl_i + \
               pesticide_USD_per_sqm * area * cycl + \
               fertilizer_USD_per_sqm * area * cycl - \
               market_price_USD_per_kg * yield_srf

        capex_srf.append(capex)
        capex_a_srf.append(capex_a)
        opex_srf.append(opex)

    return capex_srf, capex_a_srf, opex_srf




