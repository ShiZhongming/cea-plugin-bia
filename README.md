# cea-plugin-bia
A repository for a CEA plugin that calculates crop yields (kg), environmental impacts including GHG Emissions (kg CO2-eq), energy (kWh) and water use (litre), costs including capital and oeprational expenditures (USD) for the selected crop type on the selected building envelope surface.

To install, clone this repo to a desired path (you would need to have `git` installed to run this command. Alternatively, you can also run this command in the CEA console, which comes with `git` pre-installed):

```git clone https://github.com/shizhongming/cea-plugin-bia.git DESIRED_PATH```


Open CEA console and enter the following command to install the plugin to CEA:

```pip install -e /Users/your_name/Documents/GitHub/cea-plugin-bia```


In the CEA console, enter the following command to enable the plugin in CEA:

```cea-config write --general:plugins bia.bia_assessment.BiaAssessmentPlugin```

Now you should be able to enter the following command to run the plugin:

```cea bia-assessment```
