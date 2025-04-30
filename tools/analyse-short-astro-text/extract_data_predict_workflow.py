#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python

# This script is generated with nb2galaxy

# flake8: noqa

import json
import os
import shutil

import pandas as pd
from astropy.table import Table
from fetch_atel import fetch_atel
from fetch_gcn import fetch_gcn
from oda_api.data_products import ODAAstropyTable
from oda_api.json import CustomJSONEncoder
from pipeline_astrobert import get_astroBERT_cleaned_result
from pipeline_create_url_vector import create_url_vector
from pipeline_ra_dec import rule_based_ra_dec_detector
from pipeline_source_classes import detect_source_classes
from pipeline_sources import query_info_sources, rule_based_source_detector
from pipeline_telescope import rule_based_telescope_detector
from pipeline_vectorize_text import vectorize_text
from predict_vectorised_text import predict_vector

text = ""  # http://odahub.io/ontology#LongString ; oda:label "Text (optional)"
number = 16672  # http://odahub.io/ontology#Integer ; oda:label "Number of ATEL or GCN"
origin_type = "ATEL"  # oda:String ; oda:allowed_value "atel","gcn" ; oda:label "Select ATEL or GCN"

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "C_data_product_" in inp_dic.keys():
    inp_pdic = inp_dic["C_data_product_"]
else:
    inp_pdic = inp_dic
text = str(inp_pdic["text"])
number = int(inp_pdic["number"])
origin_type = str(inp_pdic["origin_type"])

if origin_type == "ATEL":
    text_id = "ATel #{}".format(number)

    if text == "":
        text_id_text = fetch_atel(number)
        if text_id_text is None:
            raise ValueError("Not possible to fetch ATel number #{number}.")
    else:
        text_id_text = text
else:
    text_id = "GCN #{}".format(number)

    if text == "":
        text_id_text = fetch_gcn(number)
        if text_id_text is None:
            raise ValueError("Not possible to fetch GCN number #{number}.")
    else:
        text_id_text = text

### Settings
try:
    data_path = os.path.dirname(__file__) + "/data/"
except:
    data_path = "data/"
print(data_path)

### Required files
telescope_ontology = f"{data_path}/telescope_observatory_survey.ttl"
file_dict_sens_inst = (
    f"{data_path}/dictionary_telescope_type_2_instrument.json"
)

simbad_node_file = f"{data_path}/simbad_otypes_nodes.csv"

### Run pipeline
df_tel = rule_based_telescope_detector(
    text_id, text_id_text, telescope_ontology
)
df_astrobert = get_astroBERT_cleaned_result(text_id, text_id_text)

astrobert_sources = list(
    df_astrobert[df_astrobert["entity_group"] == "CelestialObject"]
    .dropna()
    .word.values
)
regex_sources = rule_based_source_detector(text_id, text_id_text)

df_sor, df_unk_sor = query_info_sources(
    text_id, list(set(regex_sources + astrobert_sources))
)
df_pos_0 = rule_based_ra_dec_detector(text_id, text)
df_pos = (
    pd.concat(
        [df_pos_0, df_sor],
        axis=0,
        names=["Main ID Name", "RA", "Dec"],
        ignore_index=True,
    )
    .drop(labels=["Raw Source Name", "OTYPE"], axis=1)
    .drop_duplicates(subset=["RA", "Dec"])
)

df_cla = detect_source_classes(
    text_id, text_id_text.lower(), df_sor, simbad_node_file
)

df_vec = vectorize_text(text_id, data_path, df_tel, df_sor, df_astrobert)

df_vec_init_pred = predict_vector(text_id, data_path, df_vec)
df_vec_url, df_url_scores = create_url_vector(
    text_id, data_path, file_dict_sens_inst, df_vec_init_pred, df_sor
)

t_tel = ODAAstropyTable(
    Table.from_pandas(df_tel.map(str).map(lambda x: x.encode("utf-8").strip()))
)
t_sor = ODAAstropyTable(
    Table.from_pandas(df_sor.map(str).map(lambda x: x.encode("utf-8").strip()))
)
t_pos = ODAAstropyTable(
    Table.from_pandas(df_pos.map(str).map(lambda x: x.encode("utf-8").strip()))
)
t_unk_sor = ODAAstropyTable(
    Table.from_pandas(
        df_unk_sor.map(str).map(lambda x: x.encode("utf-8").strip())
    )
)
t_cla = ODAAstropyTable(
    Table.from_pandas(df_cla.map(str).map(lambda x: x.encode("utf-8").strip()))
)
t_astrobert = ODAAstropyTable(
    Table.from_pandas(
        df_astrobert.map(str).map(lambda x: x.encode("utf-8").strip())
    )
)
t_vec_init_pred = ODAAstropyTable(
    Table.from_pandas(
        df_vec_init_pred.map(str).map(lambda x: x.encode("utf-8").strip())
    )
)
t_vec_url = ODAAstropyTable(
    Table.from_pandas(
        df_vec_url.map(str).map(lambda x: x.encode("utf-8").strip())
    )
)
t_url_scores = ODAAstropyTable(
    Table.from_pandas(
        df_url_scores.map(str).map(lambda x: x.encode("utf-8").strip())
    )
)

print(df_sor)
print(df_unk_sor)
print(df_pos_0)
print(df_pos)
print(df_tel)
print(df_cla)
print(df_astrobert)
print(df_vec_init_pred)

table_telescopes = t_tel  # http://odahub.io/ontology#ODAAstropyTable
table_sources = t_sor  # http://odahub.io/ontology#ODAAstropyTable
table_source_positions = t_pos  # http://odahub.io/ontology#ODAAstropyTable
table_unknown_sources = t_unk_sor  # http://odahub.io/ontology#ODAAstropyTable
table_source_classes = t_cla  # http://odahub.io/ontology#ODAAstropyTable
table_astrobert_results = (
    t_astrobert  # http://odahub.io/ontology#ODAAstropyTable
)
table_vectorized_text = (
    t_vec_init_pred  # http://odahub.io/ontology#ODAAstropyTable
)
table_vectorized_url = t_vec_url  # http://odahub.io/ontology#ODAAstropyTable
table_vectorized_url_scores = (
    t_url_scores  # http://odahub.io/ontology#ODAAstropyTable
)

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_extract_data_predict_workflow_table_telescopes",
        "table_telescopes_galaxy.output",
        table_telescopes,
    )
)
_oda_outs.append(
    (
        "out_extract_data_predict_workflow_table_sources",
        "table_sources_galaxy.output",
        table_sources,
    )
)
_oda_outs.append(
    (
        "out_extract_data_predict_workflow_table_source_positions",
        "table_source_positions_galaxy.output",
        table_source_positions,
    )
)
_oda_outs.append(
    (
        "out_extract_data_predict_workflow_table_unknown_sources",
        "table_unknown_sources_galaxy.output",
        table_unknown_sources,
    )
)
_oda_outs.append(
    (
        "out_extract_data_predict_workflow_table_source_classes",
        "table_source_classes_galaxy.output",
        table_source_classes,
    )
)
_oda_outs.append(
    (
        "out_extract_data_predict_workflow_table_astrobert_results",
        "table_astrobert_results_galaxy.output",
        table_astrobert_results,
    )
)
_oda_outs.append(
    (
        "out_extract_data_predict_workflow_table_vectorized_text",
        "table_vectorized_text_galaxy.output",
        table_vectorized_text,
    )
)
_oda_outs.append(
    (
        "out_extract_data_predict_workflow_table_vectorized_url",
        "table_vectorized_url_galaxy.output",
        table_vectorized_url,
    )
)
_oda_outs.append(
    (
        "out_extract_data_predict_workflow_table_vectorized_url_scores",
        "table_vectorized_url_scores_galaxy.output",
        table_vectorized_url_scores,
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
