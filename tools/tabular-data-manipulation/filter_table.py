#!/usr/bin/env python
# coding: utf-8

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

# flake8: noqa

import json
import os

fn = "testfile.tsv"  # oda:POSIXPath
new_column = "sum"
expression = "c1 > 1e3"
sep = "comma"

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

df = df[filter_df(df)]

df.to_csv("outfile.tsv", sep=sep, index=False)

outputfile = "outfile.tsv"  # http://odahub.io/ontology#test

# output gathering
_galaxy_meta_data = {}

with open(os.path.join(_galaxy_wd, "galaxy.json"), "w") as fd:
    json.dump(_galaxy_meta_data, fd)
print("*** Job finished successfully ***")
