#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python

# This script is generated with nb2galaxy

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

fn = "data.tsv"  # oda:POSIXPath
skiprows = 0  # http://odahub.io/ontology#Integer
sep = "whitespace"  # http://odahub.io/ontology#String ; oda:allowed_value "auto", "comma", "tab", "whitespace", "semicolon"

ra_col = "c3"  # http://odahub.io/ontology#String
dec_col = "c4"  # http://odahub.io/ontology#String
weight_col = ""  # http://odahub.io/ontology#String
binsz = 0.02  # http://odahub.io/ontology#Float
window_size_RA = 2.0  # http://odahub.io/ontology#Degree
window_size_DEC = 2.0  # http://odahub.io/ontology#Degree

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "C_data_product_" in inp_dic.keys():
    inp_pdic = inp_dic["C_data_product_"]
else:
    inp_pdic = inp_dic
fn = str(inp_pdic["fn"])
skiprows = int(inp_pdic["skiprows"])
sep = str(inp_pdic["sep"])
ra_col = str(inp_pdic["ra_col"])
dec_col = str(inp_pdic["dec_col"])
weight_col = str(inp_pdic["weight_col"])
binsz = float(inp_pdic["binsz"])
window_size_RA = float(inp_pdic["window_size_RA"])
window_size_DEC = float(inp_pdic["window_size_DEC"])

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from gammapy.maps import Map
from oda_api.data_products import ImageDataProduct, PictureProduct

separators = {
    "tab": "\t",
    "comma": ",",
    "semicolon": ";",
    "whitespace": "\s+",
    "space": " ",
}

df = None

if sep == "auto":
    for name, s in separators.items():
        try:
            df = pd.read_csv(fn, sep=s, index_col=False, skiprows=skiprows)
            if len(df.columns) > 2:
                sep = s
                print("Detected separator: ", name)
                break
        except Exception as e:
            print("Separator ", s, " failed", e)
    assert sep != "auto", "Failed to find valid separator"

if df is None:
    df = pd.read_csv(fn, sep=separators[sep], index_col=False)

df.columns

def read_data(df, colname, optional=False):
    for i, c in enumerate(df.columns):
        if colname == f"c{i+1}":
            print(colname, c)
            return df[c].values
        elif colname == c:
            print(colname, c)
            return df[c].values

    assert optional, colname + " column not found"
    return None

ra = read_data(df, ra_col)
dec = read_data(df, dec_col)
w = read_data(df, weight_col, optional=True)
if w is None:
    w = np.ones_like(ra)

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

plot = PictureProduct.from_file("map.png")

plot = plot  # http://odahub.io/ontology#ODAPictureProduct
fits_image = fits_image  # http://odahub.io/ontology#Image

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(("out_sky_plot_plot", "plot_galaxy.output", plot))
_oda_outs.append(
    ("out_sky_plot_fits_image", "fits_image_galaxy.output", fits_image)
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
