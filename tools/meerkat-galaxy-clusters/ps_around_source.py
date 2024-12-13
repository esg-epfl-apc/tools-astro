#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import time
import warnings

warnings.filterwarnings("ignore")

from catalogs_methods import *
from data_treatment import *
from plots import *

ra = 83.7  # http://odahub.io/ontology#AngleDegrees
dec = -11  # http://odahub.io/ontology#AngleDegrees
thresh_arcmin = 1  # http://odahub.io/ontology#arcmin

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in ["ra", "dec", "thresh_arcmin"]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

sources = find_sources_around_coordinates(ra, dec, thresh_arcmin)
msg_out = sources["out_msg"]  # http://odahub.io/ontology#String
only_spectrum = True
print(msg_out)

t1 = time.time()
plot_power_spectrum(ra, dec, thresh_arcmin, only_spectrum)
t2 = time.time()
print(t2 - t1)

# output gathering
_galaxy_meta_data = {}

with open(os.path.join(_galaxy_wd, "galaxy.json"), "w") as fd:
    json.dump(_galaxy_meta_data, fd)
print("*** Job finished successfully ***")
