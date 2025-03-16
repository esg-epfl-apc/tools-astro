#!/usr/bin/env python
# coding: utf-8

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

# flake8: noqa

import json
import os

# import warnings
# warnings.filterwarnings('ignore')

name_clust = "Abell 133"

# # ----------------------------------------------------------------------------------
# # Plotting of the enhanced image, the target location and the corresponding sources
# # ----------------------------------------------------------------------------------
# freq_index = 0
# plot_enhanced_image(clust_name, freq_index)
# plt.show()

# output gathering
_galaxy_meta_data = {}

with open(os.path.join(_galaxy_wd, "galaxy.json"), "w") as fd:
    json.dump(_galaxy_meta_data, fd)
print("*** Job finished successfully ***")
