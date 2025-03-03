#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

# import ctadata.direct

src_name = "1ES 1218+304"  # http://odahub.io/ontology#AstrophysicalObject
RA = 185.3414279  # http://odahub.io/ontology#PointOfInterestRA
DEC = 30.17698848  # http://odahub.io/ontology#PointOfInterestDEC
T1 = "2020-01-01T00:00:00.0"  # http://odahub.io/ontology#StartTime
T2 = "2022-04-01T00:00:00.0"  # http://odahub.io/ontology#EndTime
Radius = 15.0  # http://odahub.io/ontology#AngleDegrees ; oda:label "Search cone radius"

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in ["src_name", "RA", "DEC", "T1", "T2", "Radius"]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

run_list = [1, 2, 3]

run_list = run_list  # http://odahub.io/ontology#ODAAstropyTable

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    ("out_Find_data_run_list", "run_list_galaxy.output", run_list)
)

for _outn, _outfn, _outv in _oda_outs:
    _galaxy_outfile_name = os.path.join(_galaxy_wd, _outfn)
    if isinstance(_outv, str) and os.path.isfile(_outv):
        shutil.move(_outv, _galaxy_outfile_name)
        _galaxy_meta_data[_outn] = {"ext": "_sniff_"}
    elif getattr(_outv, "write_fits_file", None):
        _outv.write_fits_file(_galaxy_outfile_name)
        _galaxy_meta_data[_outn] = {"ext": "fits"}
    elif getattr(_outv, "write_file", None):
        _outv.write_file(_galaxy_outfile_name)
        _galaxy_meta_data[_outn] = {"ext": "_sniff_"}
    else:
        with open(_galaxy_outfile_name, "w") as fd:
            json.dump(_outv, fd, cls=CustomJSONEncoder)
        _galaxy_meta_data[_outn] = {"ext": "json"}

with open(os.path.join(_galaxy_wd, "galaxy.json"), "w") as fd:
    json.dump(_galaxy_meta_data, fd)
print("*** Job finished successfully ***")
