# cea-plugin-bia
A repository for a CEA plugin that calculates metrics in building-integrated agriculture.

To install, clone this repo to a desired path (you would need to have `git` installed to run this command. Alternatively you can also run this command in the CEA console, which
comes with `git` pre-installed):

```git clone https://github.com/shizhongming/cea-plugin-bia.git DESIRED_PATH```


Open CEA console and enter the following command to install the plugin to CEA:

```pip install -e PATH_OF_PLUGIN_FOLDER```

(NOTE: PATH_OF_PLUGIN_FOLDER would be the DESIRED_PATH + 'cea-plugin-bia')


In the CEA console, enter the following command to enable the plugin in CEA:

```cea-config write --general:plugins bia.bia_dli.BiaDliPlugin```

Now you should be able to enter the following command to run the plugin:

```cea bia-dli```

NOTE: When creating your own plugin based on this template, you'll have to adjust the installation instructions above to match.
NOTE: When installing multiple plugins, add them as a comma separated list in the `cea-config write --general:plugins ...` command.
