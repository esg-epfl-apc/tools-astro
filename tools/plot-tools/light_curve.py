#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

fn = "testfile_l.tsv"  # oda:POSIXPath
sep = "whitespace"  # http://odahub.io/ontology#String ; oda:allowed_value "comma", "tab", "space", "whitespace", "semicolon"
column = "c5"  # http://odahub.io/ontology#String
weights_column = ""  # http://odahub.io/ontology#String
binning = "logarithmic"  # http://odahub.io/ontology#String ; oda:allowed_value "linear","logarithmic"
minval = 0  # http://odahub.io/ontology#Float
maxval = 0  # http://odahub.io/ontology#Float
use_quantile_values = False  # https://odahub.io/ontology/#Boolean
nbins = 15  # http://odahub.io/ontology#Integer
xlabel = "time, s"  # http://odahub.io/ontology#String
ylabel = "Ncounts"  # http://odahub.io/ontology#String

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
    "weights_column",
    "binning",
    "minval",
    "maxval",
    "use_quantile_values",
    "nbins",
    "xlabel",
    "ylabel",
]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

assert minval >= 0 or not use_quantile_values
assert maxval >= 0 or not use_quantile_values
assert minval <= 1 or not use_quantile_values
assert maxval <= 1 or not use_quantile_values
assert minval < maxval or minval == 0 or maxval == 0

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

def weighted_quantile(
    values, quantiles, sample_weight=None, values_sorted=False, old_style=False
):
    """Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of initial array
    :param old_style: if True, will correct output to be consistent with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(
        quantiles <= 1
    ), "quantiles should be in [0, 1]"

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with np.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)

weightname = ""
colname = ""

for i, c in enumerate(df.columns):

    if column == f"c{i}":
        colname = c
    elif weights_column == f"c{i}":
        weightname = c
    elif column == c:
        colname = c
    elif weights_column == c:
        weightname = c

assert (
    len(weightname) > 0 or len(weights_column) == 0
), "weight column not found"
assert len(colname) > 0, "value column not found"

print(colname, weightname)

delays = df[colname].values

if len(weightname) > 0:
    weights = df[weightname].values
else:
    weights = np.ones_like(delays)

if binning != "linear":
    min_positive_val = np.min(delays[delays > 0])
    delays[delays <= 0] = (
        min_positive_val  # replace zero delays with minimal positive value
    )

if use_quantile_values:
    minval, maxval = weighted_quantile(
        delays, [minval, maxval], sample_weight=weights
    )
    if minval == maxval:
        print("ignoreing minval and maxval (empty range)")
        minval = np.min(delays)
        maxval = np.max(delays)
else:
    if minval == 0:
        minval = np.min(delays)
    if maxval == 0:
        maxval = np.max(delays)

if minval == maxval:
    print("correcting minval and maxval (empty range)")
    maxval = minval * 1.1 if minval > 0 else 1e-100

from numpy import log10

if binning == "linear":
    bins = np.linspace(minval, maxval, nbins + 1)
else:
    bins = np.logspace(log10(minval), log10(maxval), nbins + 1)
bins

plt.figure()
h = plt.hist(delays, weights=weights, bins=bins)

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
        "out_light_curve_histogram_data",
        "histogram_data_galaxy.output",
        histogram_data,
    )
)
_oda_outs.append(
    (
        "out_light_curve_histogram_picture",
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
