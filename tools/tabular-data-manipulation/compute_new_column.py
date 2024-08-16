#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

fn = "testfile.tsv"  # oda:POSIXPath
new_column = "sum"
expression = "c1 + c2"
sep = "comma"

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

def filter_df(row):
    for i, c in enumerate(df.columns):
        globals()[f"c{i}"] = row[c]

    return eval(expression)

df[new_column] = df.apply(filter_df, axis=1)

df

df.to_csv("outfile.tsv", sep=sep, index=False)

# from oda_api.data_products import BinaryProduct, PictureProduct
# bin_data = BinaryProduct.from_file("outfile.tsv")

outputfile = "outfile.tsv"  # http://odahub.io/ontology#test

# output gathering
_galaxy_meta_data = {}
_simple_outs = []
_simple_outs.append(
    (
        "out_compute_new_column_outputfile",
        "outputfile_galaxy.output",
        outputfile,
    )
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
