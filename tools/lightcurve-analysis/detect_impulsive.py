#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

light_curve = "spiacs_lc_query.fits"  # http://odahub.io/ontology#POSIXPath

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for vn, vv in inp_pdic.items():
    if vn != "_selector":
        globals()[vn] = type(globals()[vn])(vv)

import numpy as np
from astropy.io import fits
from oda_api.data_products import PictureProduct

lc = fits.open(light_curve)[1].data

lc

lc.columns.names

t = lc["TIME"]

for rate_name in ["RATE", "FLUX"]:
    if rate_name in lc.columns.names:
        rate = lc[rate_name]

rate_err = lc["ERROR"]

bkg = np.mean(rate)

peak_i = rate.argmax()

from matplotlib import pylab as plt

plt.figure()

plt.errorbar(t, rate - bkg, rate_err)

plt.axvline(t[peak_i], c="r")

plt.xlabel("seconds since T$_0$")
plt.ylabel("counts/s")

plt.savefig("lc.png")

lc = PictureProduct.from_file("lc.png")

detection_image = lc  # oda:ODAPictureProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_detect_impulsive_detection_image",
        "detection_image_galaxy.output",
        detection_image,
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
