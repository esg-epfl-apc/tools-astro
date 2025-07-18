#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

fn = "testfile.tsv"  # oda:POSIXPath
new_column = "sum"
expression = "c1 __gt__ (c2 + v0 - v1)"
variables = "1 2"
# expression = "c1 + c2"
sep = "auto"  # oda:allowed_value "auto", "comma", "tab"
action = "add"  # oda:allowed_value "add", "filter"

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

# this is a patch due to some anomaly in the ODA bot
for k, v in [
    ("X", "&"),
    ("__gt__", ">"),
    ("__lt__", "<"),
    ("__tc__", " "),
]:
    new_expression = expression.replace(k, v)
    if new_expression != expression:
        print(f"Replaced {k} with {v}")
        expression = new_expression

    new_variables = variables.replace(k, v)
    if new_variables != variables:
        print(f"Replaced {k} with {v}")
        variables = new_variables

import numpy as np
import pandas as pd

if sep == "tab":
    sep = "\t"
elif sep == "comma":
    sep = ","
elif sep == "auto":
    for s in [",", "\t"]:
        try:
            df = pd.read_csv(fn, sep=s, index_col=False)
            if len(df.columns) > 2:
                sep = s
                print("Detected separator: ", sep)
                break
        except Exception as e:
            print("Separator ", s, " failed", e)
    pd.read_csv(fn, sep=sep, index_col=False)

    if sep == "auto":
        raise Exception("Separator not detected")

df = pd.read_csv(fn, sep=sep, index_col=False)

for i, v in enumerate(variables.split()):
    globals()[f"v{i}"] = float(v)

def filter_df(row):
    for i, c in enumerate(df.columns):
        globals()[f"c{i}"] = row[c]

    return eval(expression)

v = df.apply(filter_df, axis=1)

if action == "add":
    df[new_column] = v
elif action == "filter":
    df = df[v.astype(bool)]

df

df.to_csv("outfile.tsv", sep="\t", index=False)

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
