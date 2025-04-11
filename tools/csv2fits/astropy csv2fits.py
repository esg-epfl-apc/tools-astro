#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python

# This script is generated with nb2galaxy

# flake8: noqa

import json
import os
import shutil
import subprocess

import pandas as pd
from astropy.table import Table
from oda_api.json import CustomJSONEncoder

input_file = "https://www.astro.unige.ch/~tucci/Phosphoros/MultiBands_Catalog_1k.fits"  # http://odahub.io/ontology#FileReference

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "C_data_product_" in inp_dic.keys():
    inp_pdic = inp_dic["C_data_product_"]
else:
    inp_pdic = inp_dic
input_file = str(inp_pdic["input_file"])

workdir = os.getcwd()
path_tmp = workdir + "/tmp/"
os.makedirs(path_tmp, exist_ok=True)

# get the input catalog and save it into tmp/ directory as Input_Catalog.fits
csv_file = path_tmp + "file.csv"

read_from_url = False
try:
    desi_output = subprocess.check_output(
        "cp " + input_file + " " + csv_file, shell=True
    ).decode()
except:
    raise RuntimeError("NOT a file")

df = pd.read_csv(csv_file, delimiter="\s+")

t = Table.from_pandas(df)

from oda_api.data_products import ODAAstropyTable

cat = ODAAstropyTable(t)

catalog_table = cat  # http://odahub.io/ontology#ODAAstropyTable

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_astropy_csv2fits_catalog_table",
        "catalog_table_galaxy.output",
        catalog_table,
    )
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
