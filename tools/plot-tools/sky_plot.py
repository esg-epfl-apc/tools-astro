#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

ra_col = 0  # http://odahub.io/ontology#Integer
dec_col = 1  # http://odahub.io/ontology#Integer
weight_col = 2  # http://odahub.io/ontology#Integer
data_file = "data.tsv"  # oda:POSIXPath
binsz = 0.02  # http://odahub.io/ontology#Float
window_size_RA = 2.0  # http://odahub.io/ontology#Degree
window_size_DEC = 2.0  # http://odahub.io/ontology#Degree
skiprows = 0  # http://odahub.io/ontology#Integer

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

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from gammapy.maps import Map
from oda_api.data_products import ImageDataProduct

def create_test_data():
    Nevents = 1000
    ra = 45 + np.random.randn(Nevents)
    dec = 45 + np.random.randn(Nevents)
    w = np.abs(np.random.randn(Nevents)) + 0.1
    data = np.vstack([ra, dec, w]).transpose()
    np.savetxt(data_file, data)

# create_test_data()

data = np.loadtxt(data_file, skiprows=skiprows)
w = data[:, weight_col]
ra = data[:, ra_col]
dec = data[:, dec_col]

source = SkyCoord(ra=np.mean(ra) * u.deg, dec=np.mean(dec) * u.deg)

map = Map.create(
    binsz=binsz,
    width=(window_size_RA * u.deg, window_size_DEC * u.deg),
    frame="icrs",
    axes=[],
    skydir=SkyCoord(source),
)

map.fill_by_coord({"lat": dec * u.deg, "lon": ra * u.deg}, weights=w)

map.plot()
plt.savefig("map.png")

map.write("map.fits", overwrite=True)
fits_image = ImageDataProduct.from_fits_file("map.fits")

# plot = PictureProduct.from_file('map.png')
# fits_image=fits_image #  # http://odahub.io/ontology#Image

plot = plot  # http://odahub.io/ontology#ODAPictureProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(("out_sky_plot_plot", "plot_galaxy.output", plot))

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
