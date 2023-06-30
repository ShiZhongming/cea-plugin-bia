from setuptools import setup, find_packages

__author__ = "Zhongming Shi"
__copyright__ = "Copyright 2023, Future Cities Laboratory, Singapore - ETH Zurich; " \
                "University of Calgary, Alberta, Canada"
__credits__ = ["Zhongming Shi"]
__license__ = "MIT"
__version__ = "1.4"
__maintainer__ = "Zhongming Shi"
__email__ = "shi@arch.ethz.ch"
__status__ = "Production"


setup(name='cea_plugin_bia',
      version=__version__,
      description="A plugin for the City Energy Analyst: building-integrated agriculture extensions",
      license='MIT',
      author='Zhongming Shi',
      author_email='shi@arch.ethz.ch',
      url='https://github.com/shizhongming/cea-plugin-bia',
      long_description="A plugin for the City Energy Analyst: building-integrated agriculture (BIA) extensions."
                       "It has two main functions. One calculates crop yields (kg), environmental impacts "
                       "including GHG emissions (kg CO2-eq), energy (kWh) and water use (litre), costs including "
                       "capital and operational expenditures (USD) for the selected crop type on the selected "
                       "building envelope surface. The other one produces the crop profile and planting calendar"
                       " for each building surface, based on one of the user-defined objectives."
                       "As of August 1, 2022, this plugin works the best for Singapore or its adjacent regions "
                       "as the planting data have been acquired from the Tropical Technologies (T2) Laboratory "
                       "affiliated to the National University of Singapore in Clementi, Singapore.",
      py_modules=[''],
      packages=find_packages(),
      package_data={},
      include_package_data=True)
