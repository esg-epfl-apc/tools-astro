#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os

from catalogs_methods import *
from data_treatment import *
from download import *
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

# -----------------------------------------------------
# Handling errors and exceptions
# -----------------------------------------------------
# if ra<0 or ra>360:
#     raise ValueError("Wrong value: the right ascension (ra) must be between 0 and 360 degrees!")
# if dec<-90 or dec>90:
#     raise ValueError("Wrong value: the declination (dec) must be between -90 and 90 degrees!")

# # Test if the ra/dec coordinates lie in a MGCLS field (cluster)
# # -----------------------------------------------------------------------------
# clust_name = test_if_in_clusters(ra,dec)

# if clust_name == 0:

#     raise Exception("Error: the entered coordinates do not lie within any MGCLS field (cluster)")

# ----------------------------------------------------------------------------------
# Determination of the sources found (or not) in the given aperture
# ----------------------------------------------------------------------------------
sources = find_sources_around_coordinates(ra, dec, thresh_arcmin)
msg_out = sources["out_msg"]

# ----------------------------------------------------------------------------------
# Computation of the power spectrum of the sources within the aperture
# ----------------------------------------------------------------------------------
# if sources['found_source'] == True:
#     plot_spectrum(ra, dec, thresh_arcmin)
#     plt.savefig('spectrum.png',format='png',bbox_inches='tight')

# if sources['found_source'] == True:
#     bin_image = PictureProduct.from_file('spectrum.png')
# if sources['found_source'] == True:
#     picture = bin_image # http://odahub.io/ontology#ODAPictureProduct

test_out = msg_out  # http://odahub.io/ontology#String

# output gathering
_galaxy_meta_data = {}

with open(os.path.join(_galaxy_wd, "galaxy.json"), "w") as fd:
    json.dump(_galaxy_meta_data, fd)
print("*** Job finished successfully ***")
