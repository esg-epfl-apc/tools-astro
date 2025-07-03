#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python

# This script is generated with nb2galaxy

# flake8: noqa

import json
import os

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
    data = data.astype(data.dtype.newbyteorder("="))
except:
    raise RuntimeError("The input file should have the FITS format.")

image_out_path = "./output.tiff"
tifffile.imwrite(image_out_path, data)

file_output = image_out_path

# output gathering
_galaxy_meta_data = {}

with open(os.path.join(_galaxy_wd, "galaxy.json"), "w") as fd:
    json.dump(_galaxy_meta_data, fd)
print("*** Job finished successfully ***")
