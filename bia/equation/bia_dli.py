"""
This script calculates:
the Daily Light Integral (DLI) in [mol/m2/day] for each building envelope surface.
"""

from __future__ import division
from __future__ import print_function
import os
import time
import pyarrow.feather as feather
import pandas as pd
import geopandas as gpd


__author__ = "Zhongming Shi"
__copyright__ = "Copyright 2023, A/S Group, ITA, ETH Zurich"
__credits__ = ["Zhongming Shi"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Zhongming Shi"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"



def calc_DLI(locator, config, building_name):
    """
    This function first determines the surface area with sufficient solar radiation, and then calculates the optimal
    tilt angles of panels at each surface location. The panels are categorized into groups by their surface azimuths,
    tilt angles, and global irradiation. In the last, electricity generation from PV panels of each group is calculated.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator
    :param radiation_path: solar insulation data on all surfaces of each building (path
    :type radiation_path: String
    :param metadata_csv: data of sensor points measuring solar insulation of each building
    :type metadata_csv: .csv
    :param latitude: latitude of the case study location
    :type latitude: float
    :param longitude: longitude of the case study location
    :type longitude: float
    :param weather_path: path to the weather data file of the case study location
    :type weather_path: .epw
    :param building_name: list of building names in the case study
    :type building_name: Series
    :return: Building_PV.csv with PV generation potential of each building, Building_sensors.csv with sensor data of
        each PV panel.

    """

    t0 = time.perf_counter()
    radiation_path = locator.get_radiation_building_sensors(building_name)
    metadata_csv_path = locator.get_radiation_metadata(building_name)

    # select sensor point with sufficient solar radiation
    max_annual_radiation, annual_radiation_threshold, sensors_rad_clean, sensors_metadata_clean = \
        filter_low_potential(radiation_path, metadata_csv_path, config)


    if not sensors_metadata_clean.empty:
        # convert solar radiation to DLI
        sensors_DLI = calc_Whperm2_molperm2(sensors_rad_clean).T

        # label the sensors by their #floor and wall type (lower, upper, and sideX2)
        sensors_wall_type = calc_sensor_wall_type(locator, sensors_metadata_clean, building_name)
        sensors_metadata_clean['wall_type'] = sensors_wall_type

        # merge the calculated results
        sensors_metadata_clean_DLI = pd.merge(sensors_metadata_clean, sensors_DLI, left_index=True,
                                                    right_index=True, how="left")
        sensors_metadata_clean_DLI.reset_index(inplace=True)
        sensors_metadata_clean_DLI = sensors_metadata_clean_DLI.rename(columns={'index': 'srf_index'})

        # write the DLI results
        dir = config.scenario + "/outputs/data/potentials/agriculture/dli"
        if not os.path.exists(dir):
            os.makedirs(dir)
        output_path = dir + "/{building}_DLI.csv".format(building=building_name)
        sensors_metadata_clean_DLI.to_csv(output_path, index=False,
                                          float_format='%.2f',
                                          na_rep=0)  # write sensors metadata and DLI

        print('Calculations of DLI for each sensor on Building', building_name, 'done - time elapsed: %.2f seconds'
              % (time.perf_counter() - t0))

    else:  # This loop is activated when a building has not sufficient solar potential
        print("Unfortunately, Building", building_name, "has no BIA potential.")
        dir = config.scenario + "/outputs/data/potentials/agriculture/dli"
        if not os.path.exists(dir):
            os.mkdir(dir)

        pass


def calc_Whperm2_molperm2(radiation_Whperm2):
    """
    To calculate the total number of photons (in the 400-700 range) that reach the building surface

    :param absorbed_radiation_Wperm2: absorbed radiation [W/m2]
    :type absorbed_radiation_Wperm2: float

    :return dli_output_mol_h: dli per m2 per day [mol/m2/day]
    :rtype dli_output_mol_h: dataframe

    references:
    ..https://www.controlledenvironments.org/wp-content/uploads/sites/6/2017/06/Ch01.pdf.

    """
    wl_400_700nm = 0.45     # only about 45% of the energy of solar radiation is actually in the 400 - 700 nm range
    conversion = 4.57 * 3600   # from wh/m2 to μmol/m2/h
    dli_output_mol_h = radiation_Whperm2 * wl_400_700nm * conversion / 10**6

    # label the 8760 hours by 365 days
    day = pd.Series(range(0, 365))
    hour_to_day = day.repeat(24).reset_index().pop('index')
    dli_output_mol_h = pd.merge(dli_output_mol_h, hour_to_day, left_index=True, right_index=True, how="left")

    # calculate from hourly dli to dli
    dli_output_mol_day = dli_output_mol_h.groupby(['index']).sum().reset_index()
    del dli_output_mol_day['index']

    return dli_output_mol_day


def calc_building_height_info(locator):
    """
    To get the data of each building's height information, including building name, unit floor height,
    and total number of floors.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator

    :return: floor_number: data of each building's unit floor height (floor_to_floor_height)
            and total number of floors (n_floors).
    :rtype sensors_metadata_clean: dataframe

    """

    zone_path = locator.get_zone_geometry()
    zone_df = gpd.GeoDataFrame.from_file(zone_path)

    height = zone_df['height_ag'].astype(float)
    nfloors = zone_df['floors_ag'].astype(int)
    floor_to_floor_height = height / nfloors

    zone_df['floor_to_floor_height'] = floor_to_floor_height

    return zone_df


def calc_sensor_floor_number(locator, sensors_metadata_clean, building_name):

    """
    To get the floor number of each sensor.
    Attention!!! The floor numbers start from 1, not 0. The surfaces on the roof is labeled as total floor number + 1.

    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator
    :param sensors_metadata_clean: data of filtered sensor points measuring solar insulation of each building
    :type sensors_metadata_clean: dataframe
    :param building_name: list of building names in the case study
    :type building_name: Series

    :return: floor_number: list of numbers indicating the floor number of each sensor.
    :rtype sensors_metadata_clean: series

    """

    building_height_info = calc_building_height_info(locator)  # get the floor heights

    # filter the "top" surfaces and keep the sensors on facade only.
    # Facade sensors include both window and wall sensors.
    facades = sensors_metadata_clean[sensors_metadata_clean['orientation'] != 'top']


    # get the total floor numbers of the building being calculated
    n_floors = int(building_height_info['floors_ag'][building_height_info['Name'] == building_name])

    # calculate the number of facade sensors on each floor
    # n_sensors_each_floor = int(len(facades) // n_floors)

    # label the sensors with the floor number, which starts from 1
    facades_sorted = facades.sort_values(by=['Zcoor'])
    facades_sorted = pd.cut(facades_sorted['Zcoor'], bins=n_floors, labels=False)+1
    facades_sorted_df = facades_sorted.to_frame().rename(columns={'Zcoor': 'n_floor'})

    # reorder the '#floor' column to match the original order in sensors_metadata_clean by 'SURFACE'
    sensors_metadata_clean = pd.merge(sensors_metadata_clean, facades_sorted_df, left_index=True, right_index=True,
                                      how="left")
    sensors_metadata_clean['n_floor'].fillna(999, inplace=True)    # label the top sensors with 999

    floor_number = sensors_metadata_clean['n_floor']

    return floor_number


def calc_sensor_wall_type(locator, sensors_metadata_clean, building_name):
    """
    To get the floor number of each sensor. Attention! The floor numbers start from 1, not 0.
    :param locator: An InputLocator to locate input files
    :type locator: cea.inputlocator.InputLocator
    :param sensors_metadata_clean: data of filtered sensor points measuring solar insulation of each building
    :type sensors_metadata_clean: dataframe
    :param building_name: list of building names in the case study
    :type building_name: Series

    :return: sensors_wall_type: list of each sensor's wall type (lower, upper, and sideX2)
    :rtype sensors_metadata_clean: series

    """

    # label the sensors by their floor number
    sensors_floor_number = calc_sensor_floor_number(locator, sensors_metadata_clean, building_name)
    sensors_metadata_clean['n_floor'] = sensors_floor_number

    # label surface type by the four orientations
    orientation = ['north', 'east', 'south', 'west']
    results = []
    for i in orientation:
        # filter out surfaces that are not the selected orientation
        surfaces = sensors_metadata_clean[sensors_metadata_clean['orientation'] == i]

        # get the total number of floors and number of surfaces
        n_floors = max(sensors_metadata_clean['n_floor'])  # roof top is labeled as total floor number + 1

        # label surface type by the floor numbers
        results_n = []

        for j in range(1, int(n_floors)+1):
            # separate sensors on the walls and the windows of this floor
            walls = surfaces[(surfaces['n_floor'] == j) & (surfaces['TYPE'] == 'walls')].copy()
            windows = surfaces[(surfaces['n_floor'] == j) & (surfaces['TYPE'] == 'windows')].copy()

            if not windows.empty:
                # get the z-coordinates of window sensors
                sensors_windows_Zcoor = windows['Zcoor'].median()
                # label wall sensors by comparing their z-coordinates with that of the window sensors

                walls.loc[walls['Zcoor'] > sensors_windows_Zcoor, 'wall_type'] = 'upper'
                walls.loc[walls['Zcoor'] < sensors_windows_Zcoor, 'wall_type'] = 'lower'
                walls['wall_type'].fillna('side', inplace=True)

                results_n.append(walls)

        merged_results_n_df = pd.concat(results_n)
        results.append(merged_results_n_df)

    merged_results_df = pd.concat(results)

    # reorder the 'wall_type' column to match the original order in sensors_metadata_clean by 'SURFACE'
    sensors_wall_type = \
    pd.merge(sensors_metadata_clean, merged_results_df, left_index=True, right_index=True, how="left")['wall_type'] \

    # label the wall_type of windows as "non_wall"
    sensors_wall_type.fillna('non_wall', inplace=True)

    return sensors_wall_type


def filter_low_potential(radiation_path, metadata_csv_path, config):
    """
    To filter the sensor points/hours with low radiation potential.

    #. keep sensors above min radiation
    #. eliminate points when hourly production < 50 W/m2
    #. augment the solar radiation due to differences between panel reflectance and original reflectances used in daysim

    :param radiation_csv: solar insulation data on all surfaces of each building
    :type radiation_csv: .csv
    :param metadata_csv: solar insulation sensor data of each building
    :type metadata_csv: .csv

    :return max_annual_radiation: yearly horizontal radiation [Wh/m2/year]
    :rtype max_annual_radiation: float
    :return annual_radiation_threshold: minimum yearly radiation threshold for sensor selection [Wh/m2/year]
    :rtype annual_radiation_threshold: float
    :return sensors_rad_clean: radiation data of the filtered sensors [Wh/m2]
    :rtype sensors_rad_clean: dataframe
    :return sensors_metadata_clean: data of filtered sensor points measuring solar insulation of each building
    :rtype sensors_metadata_clean: dataframe

    Following assumptions are made:

    #. Sensor points with low yearly radiation are deleted. The threshold (minimum yearly radiation) is a percentage
       of global horizontal radiation. The percentage threshold (min_radiation) is a global variable defined by users.
    #. For each sensor point kept, the radiation value is set to zero when radiation value is below 50 W/m2.
    #. Unlike BIPV, window surfaces are not removed for BIA simulations
    """

    # read radiation file
    sensors_rad = feather.read_feather(radiation_path)
    sensors_metadata = pd.read_csv(metadata_csv_path)

    # join total radiation to sensor_metadata
    sensors_rad_sum = sensors_rad.sum(0).to_frame('total_rad_Whm2')  # add new row with yearly radiation
    sensors_metadata.set_index('SURFACE', inplace=True)
    sensors_metadata = sensors_metadata.merge(sensors_rad_sum, left_index=True, right_index=True)  # [Wh/m2]

    # keep sensors if allow bia installation on walls or on roofs
    if config.agriculture.crop_on_roof is False:
        sensors_metadata = sensors_metadata[sensors_metadata.TYPE != 'roofs']
    if config.agriculture.crop_on_wall_under_window is False:
        sensors_metadata = sensors_metadata[sensors_metadata.wall_type != 'upper' or 'lower']
    if config.agriculture.crop_on_wall_between_window is False:
        sensors_metadata = sensors_metadata[sensors_metadata.wall_type != 'side']

    # set min yearly radiation threshold for sensor selection
    # keep sensors above min production in sensors_rad
    max_annual_radiation = sensors_rad_sum.max().values[0]
    annual_radiation_threshold_Whperm2 = 0
    sensors_metadata_clean = sensors_metadata[sensors_metadata.total_rad_Whm2 >= annual_radiation_threshold_Whperm2]
    sensors_rad_clean = sensors_rad[sensors_metadata_clean.index.tolist()]  # keep sensors above min radiation

    return max_annual_radiation, annual_radiation_threshold_Whperm2, sensors_rad_clean, sensors_metadata_clean


