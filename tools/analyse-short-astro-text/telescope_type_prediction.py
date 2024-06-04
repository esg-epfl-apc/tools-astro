#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

from astropy.table import Table
from oda_api.json import CustomJSONEncoder
from pipeline_predict_sensitivity import predict_sensitivity
from pipeline_source_classes import detect_source_classes
from pipeline_sources import rule_based_source_detector
from pipeline_telescope import rule_based_telescope_detector

text = "The TeV-detected BL Lac object 1ES 1218+304 (z0.184) Swift"  # http://odahub.io/ontology#LongString
number = "ATEL16"  # http://odahub.io/ontology#LongString

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

atel_text = text
atel_ = number
### Settings
data_path = f"data/"
model_file = (
    f"{data_path}/probabilistic_models/first_single_second_single.dat.npy"
)

### Required files
telescope_ontology = f"{data_path}/telescope_observatory_survey.ttl"
file_dict_uri_int = f"{data_path}/dictionary_telescope_uri_int_id.json"

dict_path = f"{data_path}/dictionary_names_otypes.json"
simbad_node_file = f"{data_path}/simbad_otypes_nodes.csv"

### Run pipeline
df_tel = rule_based_telescope_detector(
    atel_, atel_text, telescope_ontology, file_dict_uri_int
)
df_sor = rule_based_source_detector(
    atel_, atel_text.lower(), data_path, dict_path
)
df_cla = detect_source_classes(
    atel_, atel_text.lower(), df_sor, simbad_node_file
)
df_in, df_pred = predict_sensitivity(model_file, df_tel, first_type="single")

t_tel = Table.from_pandas(df_tel.map(str))
t_sor = Table.from_pandas(df_sor.map(str))
t_cla = Table.from_pandas(df_cla.map(str))
t_in = Table.from_pandas(df_in.map(str))
t_pred = Table.from_pandas(df_pred.map(str))

telescope_astropy_table = t_tel  # http://odahub.io/ontology#ODAAstropyTable
source_astropy_table = t_sor  # http://odahub.io/ontology#ODAAstropyTable
source_class_astropy_table = t_cla  # http://odahub.io/ontology#ODAAstropyTable
input_sensitivity_astropy_table = (
    t_in  # http://odahub.io/ontology#ODAAstropyTable
)
pred_sensitivity_astropy_table = (
    t_pred  # http://odahub.io/ontology#ODAAstropyTable
)

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_telescope_type_prediction_telescope_astropy_table",
        "telescope_astropy_table_galaxy.output",
        telescope_astropy_table,
    )
)
_oda_outs.append(
    (
        "out_telescope_type_prediction_source_astropy_table",
        "source_astropy_table_galaxy.output",
        source_astropy_table,
    )
)
_oda_outs.append(
    (
        "out_telescope_type_prediction_source_class_astropy_table",
        "source_class_astropy_table_galaxy.output",
        source_class_astropy_table,
    )
)
_oda_outs.append(
    (
        "out_telescope_type_prediction_input_sensitivity_astropy_table",
        "input_sensitivity_astropy_table_galaxy.output",
        input_sensitivity_astropy_table,
    )
)
_oda_outs.append(
    (
        "out_telescope_type_prediction_pred_sensitivity_astropy_table",
        "pred_sensitivity_astropy_table_galaxy.output",
        pred_sensitivity_astropy_table,
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
