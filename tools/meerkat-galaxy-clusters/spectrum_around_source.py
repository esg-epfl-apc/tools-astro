#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
from catalogs_methods import *
from data_treatment import *
from download import *
from oda_api.data_products import PictureProduct
from oda_api.json import CustomJSONEncoder
from parameters import *
from plots import *

# import warnings
# warnings.filterwarnings('ignore')

# ------------------------------------------------------------------------------------
# ra/dec coordinates of target point of the sky, and the aperture to find sources in
# ------------------------------------------------------------------------------------
ra = 260  # http://odahub.io/ontology#PointOfInterestRA
dec = -82  # http://odahub.io/ontology#PointOfInterestDEC
thresh_arcmin = (
    5  # http://odahub.io/ontology#arcmin, oda:label "Radius [arcmin]"
)

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in ["ra", "dec", "thresh_arcmin"]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

# Handling errors and exceptions
# -----------------------------------------------------
if ra < 0 or ra > 360:
    raise ValueError(
        "Wrong value: the right ascension (ra) must be between 0 and 360 degrees!"
    )
if dec < -90 or dec > 90:
    raise ValueError(
        "Wrong value: the declination (dec) must be between -90 and 90 degrees!"
    )

# Test if the ra/dec coordinates lie in a MGCLS field (cluster)
# -----------------------------------------------------------------------------
clust_name = test_if_in_clusters(ra, dec)

if clust_name == 0:

    raise Exception(
        "Error: the entered coordinates do not lie within any MGCLS field (cluster)"
    )

else:

    sources = find_sources_around_coordinates(ra, dec, thresh_arcmin)
    msg_out = sources["out_msg"]

    if sources["found_source"] == True:

        plot_spectrum(ra, dec, thresh_arcmin)
        plt.savefig("spectrum.png", format="png", bbox_inches="tight")
        bin_image = PictureProduct.from_file("spectrum.png")

test_out = msg_out  # http://odahub.io/ontology#String ; oda:label "Message"
png = bin_image  # http://odahub.io/ontology#ODAPictureProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(("out_spectrum_around_source_png", "png_galaxy.output", png))

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
_simple_outs = []
_simple_outs.append(
    ("out_spectrum_around_source_test_out", "test_out_galaxy.output", test_out)
)
_numpy_available = True

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
