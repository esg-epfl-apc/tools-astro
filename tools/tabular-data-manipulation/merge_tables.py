#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

fn1 = "testfile.tsv"  # oda:POSIXPath
fn2 = "testfile.tsv"  # oda:POSIXPath
new_column = "sum"
sep = "auto"  # oda:allowed_value "auto", "comma", "tab"

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
elif sep == "auto":
    for s in [",", "\t"]:
        try:
            df = pd.read_csv(fn1, sep=s, index_col=False)
            if len(df.columns) > 2:
                sep = s
                print("Detected separator: ", sep)
                break
        except Exception as e:
            print("Separator ", s, " failed", e)
    pd.read_csv(fn1, sep=sep, index_col=False)

    if sep == "auto":
        raise Exception("Separator not detected")

df1 = pd.read_csv(fn1, sep=sep, index_col=False)
df2 = pd.read_csv(fn2, sep=sep, index_col=False)

if len(df1) != len(df2):
    raise Exception("Dataframes have different lengths")

df = pd.concat([df1, df2], axis=1)

df.to_csv("outfile.tsv", sep="\t", index=False)

# from oda_api.data_products import BinaryProduct, PictureProduct
# bin_data = BinaryProduct.from_file("outfile.tsv")

outputfile = "outfile.tsv"  # http://odahub.io/ontology#test

# output gathering
_galaxy_meta_data = {}
_simple_outs = []
_simple_outs.append(
    ("out_merge_tables_outputfile", "outputfile_galaxy.output", outputfile)
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
