from setuptools import setup, find_packages

__author__ = "Zhongming Shi"
__copyright__ = "Copyright 2022, Future Cities Laboratory, Singapore - ETH Zurich"
__credits__ = ["Zhongming Shi"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Zhongming Shi"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"


setup(name='cea_plugin_bia',
      version=__version__,
      description="A plugin for the City Energy Analyst: building-integrated agriculture extensions",
      license='MIT',
      author='Zhongming Shi',
      author_email='cea@arch.ethz.ch',
      url='https://github.com/shizhongming/cea-plugin-bia',
      long_description="A plugin for the City Energy Analyst: building-integrated agriculture (BIA) extensions."
                       "As of March 14, 2022, this plug-in calculates the DLI of each building envelope surface"
                       "Soon-to-come functionality includes BIA crop yield calculation",
      py_modules=[''],
      packages=find_packages(),
      package_data={},
      include_package_data=True)
