#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python

# This script is generated with nb2galaxy

# flake8: noqa

import json
import os
import shutil

import tifffile
from astropy.io import fits

file_input = "./input.fits"  # oda:POSIXPath; oda:label "Input file"

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "C_data_product_" in inp_dic.keys():
    inp_pdic = inp_dic["C_data_product_"]
else:
    inp_pdic = inp_dic
file_input = str(inp_pdic["file_input"])

try:
    hdul = fits.open(file_input)
    data = hdul[0].data
    header = hdul[0].header
    data = data.astype(data.dtype.newbyteorder("="))
except:
    raise RuntimeError("The input file should have the FITS format.")

image_out_path = "./output.tiff"
tifffile.imwrite(image_out_path, data)
dict_json = dict(header)

file_output = image_out_path
header_json = dict_json

# output gathering
_galaxy_meta_data = {}
_simple_outs = []
_simple_outs.append(
    ("out_fits2tiff_file_output", "file_output_galaxy.output", file_output)
)
_simple_outs.append(
    ("out_fits2tiff_header_json", "header_json_galaxy.output", header_json)
)

try:
    import numpy as np  # noqa: E402

    _numpy_available = True
except ImportError:
    _numpy_available = False

for _outn, _outfn, _outv in _simple_outs:
    _galaxy_outfile_name = os.path.join(_galaxy_wd, _outfn)
    if isinstance(_outv, str) and os.path.isfile(_outv):
        shutil.move(_outv, _galaxy_outfile_name)
        _galaxy_meta_data[_outn] = {"ext": "_sniff_"}
    elif _numpy_available and isinstance(_outv, np.ndarray):
        with open(_galaxy_outfile_name, "wb") as fd:
            np.savez(fd, _outv)
        _galaxy_meta_data[_outn] = {"ext": "npz"}
    else:
        with open(_galaxy_outfile_name, "w") as fd:
            json.dump(_outv, fd)
        _galaxy_meta_data[_outn] = {"ext": "expression.json"}

with open(os.path.join(_galaxy_wd, "galaxy.json"), "w") as fd:
    json.dump(_galaxy_meta_data, fd)
print("*** Job finished successfully ***")
