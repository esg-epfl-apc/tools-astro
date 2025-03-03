#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os

# import ctadata.direct

src_name = "1ES 1218+304"  # http://odahub.io/ontology#AstrophysicalObject

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in ["src_name"]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

RA = 185.3414279  # http://odahub.io/ontology#PointOfInterestRA
DEC = 30.17698848  # http://odahub.io/ontology#PointOfInterestDEC
T1 = "2020-01-01T00:00:00.0"  # http://odahub.io/ontology#StartTime
T2 = "2022-04-01T00:00:00.0"  # http://odahub.io/ontology#EndTime
Radius = 15.0  # http://odahub.io/ontology#AngleDegrees ; oda:label "Search cone radius"

# output gathering
_galaxy_meta_data = {}

with open(os.path.join(_galaxy_wd, "galaxy.json"), "w") as fd:
    json.dump(_galaxy_meta_data, fd)
print("*** Job finished successfully ***")
