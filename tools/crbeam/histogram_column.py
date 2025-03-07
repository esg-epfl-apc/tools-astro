#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

fn = "testfile.tsv"  # oda:POSIXPath
sep = "comma"
column = "c0"
weights = "c1"
binning = "logarithmic"  # http://odahub.io/ontology#String ; oda:allowed_value "linear","logarithmic"
minval = 1e10
maxval = 1e13
nbins = 15
xlabel = "Energy, eV"
ylabel = "Ncounts"

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
import pandas as pd

if sep == "tab":
    sep = "\t"
elif sep == "comma":
    sep = ","

df = pd.read_csv(fn, sep=sep, index_col=False)

from numpy import log10

if binning == "linear":
    bins = np.linspace(minval, maxval, nbins + 1)
else:
    bins = np.logspace(log10(minval), log10(maxval), nbins + 1)
bins

import matplotlib.pyplot as plt

for i, c in enumerate(df.columns):
    if column == f"c{i}":
        colname = c
    elif weights == f"c{i}":
        weightname = c
print(colname, weightname)

plt.figure()
h = plt.hist(df[colname], weights=df[weightname], bins=bins)

if binning == "logarithmic":
    plt.xscale("log")
    plt.yscale("log")
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.savefig("Histogram.png", format="png", dpi=150)
hist_counts = h[0]
hist_bins = h[1]
hist_mins = hist_bins[:-1]
hist_maxs = hist_bins[1:]

from astropy.table import Table
from oda_api.data_products import ODAAstropyTable, PictureProduct

names = ("bins_min", "bins_max", "counts")
res = ODAAstropyTable(Table([hist_mins, hist_maxs, hist_counts], names=names))

plot = PictureProduct.from_file("Histogram.png")

histogram_data = res  # http://odahub.io/ontology#ODAAstropyTable
histogram_picture = plot  # http://odahub.io/ontology#ODAPictureProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_histogram_column_histogram_data",
        "histogram_data_galaxy.output",
        histogram_data,
    )
)
_oda_outs.append(
    (
        "out_histogram_column_histogram_picture",
        "histogram_picture_galaxy.output",
        histogram_picture,
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
