#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

fn = "testfile_l.tsv"  # oda:POSIXPath
sep = "whitespace"  # http://odahub.io/ontology#String ; oda:allowed_value "auto", "comma", "tab", "whitespace", "semicolon"
column = "c0"  # http://odahub.io/ontology#String
weights = "c1"  # http://odahub.io/ontology#String
binning = "logarithmic"  # http://odahub.io/ontology#String ; oda:allowed_value "linear","logarithmic"
minval = 1e10  # http://odahub.io/ontology#Float
maxval = 1e13  # http://odahub.io/ontology#Float
nbins = 15  # http://odahub.io/ontology#Integer
xlabel = "Energy, [eV]"  # http://odahub.io/ontology#String
ylabel = "Flux E^2, [eV]"  # http://odahub.io/ontology#String
spec_power = 2.0  # http://odahub.io/ontology#Float

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in [
    "fn",
    "sep",
    "column",
    "weights",
    "binning",
    "minval",
    "maxval",
    "nbins",
    "xlabel",
    "ylabel",
    "spec_power",
]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
            df = pd.read_csv(fn, sep=s, index_col=False)
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

weightname = ""
colname = ""

for i, c in enumerate(df.columns):

    if column == f"c{i}":
        colname = c
    elif weights == f"c{i}":
        weightname = c
    elif column == c:
        colname = c
    elif weights == c:
        weightname = c

assert len(weightname) > 0 or len(weights) == 0, "weight column not found"
assert len(colname) > 0, "value column not found"

print(colname, weightname)

values = df[colname].values

if len(weightname) > 0:
    weights = df[weightname].values
else:
    weights = np.ones_like(values)

values, weights

from numpy import log10

if binning == "linear":
    bins = np.linspace(minval, maxval, nbins + 1)
else:
    bins = np.logspace(log10(minval), log10(maxval), nbins + 1)
bins

# # generate test E^-alpha spec
# alpha = -1
# n_part = 100000
# log10_E = np.log10(bins[0]) + np.random.rand(n_part) * (np.log10(bins[-1]) - np.log10(bins[0]))
# values = np.power(10, log10_E)
# weights = np.power(values, alpha + 1)
# np.min(values)/minval, np.max(values)/maxval, len(values)

bin_val, _ = np.histogram(values, weights=weights, bins=bins)
len(bin_val), len(bins)
bin_width = bins[1:] - bins[:-1]
flux = bin_val / bin_width
if binning == "linear":
    spec_point = 0.5 * (bins[1:] + bins[:-1])
else:
    spec_point = np.sqrt(bins[1:] * bins[:-1])

plt.figure()
h = plt.plot(spec_point, flux * spec_point**spec_power)

if binning == "logarithmic":
    plt.xscale("log")
    plt.yscale("log")

plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.savefig("spectrum.png", format="png", dpi=150)

from astropy.table import Table
from oda_api.data_products import ODAAstropyTable, PictureProduct

names = ("bins_min", "bins_max", "flux")

res = ODAAstropyTable(Table([bins[:-1], bins[1:], flux], names=names))

plot = PictureProduct.from_file("spectrum.png")

histogram_data = res  # http://odahub.io/ontology#ODAAstropyTable
histogram_picture = plot  # http://odahub.io/ontology#ODAPictureProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_spectrum_histogram_data",
        "histogram_data_galaxy.output",
        histogram_data,
    )
)
_oda_outs.append(
    (
        "out_spectrum_histogram_picture",
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
