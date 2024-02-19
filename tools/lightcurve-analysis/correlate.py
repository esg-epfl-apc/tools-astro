#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

light_curve_1 = "spiacs_lc_query.fits"  # http://odahub.io/ontology#POSIXPath
light_curve_2 = "spiacs_lc_query.fits"  # http://odahub.io/ontology#POSIXPath
mode = (
    "Background"  # oda:String; oda:allowed_value 'Background','Peak','Focus'.
)

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

from astropy.io import fits
from astropy.time import Time
from matplotlib import pylab as plt
from oda_api.data_products import LightCurveDataProduct, PictureProduct

lc1 = fits.open(light_curve_1)[1].data
lc2 = fits.open(light_curve_2)[1].data

plt.figure()

plt.plot(lc1["TIME"], lc1["RATE"])
plt.plot(lc2["TIME"], lc1["RATE"])

outlcpng = PictureProduct.from_file("lc.png")
outlc = LightCurveDataProduct.from_arrays(
    times=Time(lc1["TIME"], format="mjd"),
    rates=lc1["RATE"] - lc2["RATE"],
    errors=lc2["ERROR"],
)

detection_image = outlcpng  # oda:ODAPictureProduct
processed_lc = outlc  # oda:LightCurve

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_correlate_detection_image",
        "detection_image_galaxy.output",
        detection_image,
    )
)
_oda_outs.append(
    ("out_correlate_processed_lc", "processed_lc_galaxy.output", processed_lc)
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
