#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

from astropy.table import Table
from oda_api.data_products import ODAAstropyTable
from oda_api.json import CustomJSONEncoder
from pipeline_astrobert import get_astroBERT_cleaned_result
from pipeline_create_url_vector import create_url_vector
from pipeline_source_classes import detect_source_classes
from pipeline_sources import query_info_sources, rule_based_source_detector
from pipeline_telescope import rule_based_telescope_detector
from pipeline_vectorize_text import vectorize_text
from predict_vectorised_text import predict_vector

text = "We received an alert about a candidate fast transient from the FINK broker (Peloton, Ishida & Moller et al.) on 2024-06-22 09:58:42.001 UTC, named ZTF24aasjjkf (https://fink-portal.org/ZTF24aasjjkf). This transient was classified by FINK on the same day with a 17 percent probability of being a fast transient. We triggered the TAROT-TCA and TRT-SRO telescopes for further investigation, observing in the R, r', and i' bands, and collecting data from 2024-06-23 around 04:55:14 UTC to 2024-06-24 01:00:00 UTC. Our results, along with public ZTF data, can be found here: https://skyportal-icare.ijclab.in2p3.fr/public/sources/ZTF24aasjjkf/version/07cf393bb53513abe80d3b7863938d8f TAROT-TCA telescope observations can be also labeled as generic in the plot. We utilized STDWeb (Karpov et al.) to perform our photometry, using the Pan-STARRS DR1 (PS1) catalog and subtracting the constant flux of background with a reference PS1 image. We note that in the TAROT-TCA data, we have slight contamination from other stars. However, we mitigated this effect by subtracting the PS1 reference image. According to the ZTF data in the r' band, both before and after our observations, as well as ATLAS data (in the orange band) and the Gaia Alert (2024LXA), the source might be 10 days old or more. Its multi-band light curve is not consistent with a fast transient like a kilonova (a decay rate of 0.15 per day in r-band) but resembles a supernova (also classified by FINK on 2024-06-24 at 08:00 UTC). The source has been independently classified by Jianlin Xu et al. as a CV. We thank the FINK team for their valuable collaboration with GRANDMA. GRANDMA is a worldwide telescope network (https://grandma.lal.in2p3.fr/) devoted to the observation of transients in the context of multi-messenger astrophysics (Antier et al. 2020 MNRAS 497, 5518). Kilonova-Catcher (KNC) is the citizen science program of GRANDMA (http://kilonovacatcher.in2p3.fr/)."  # http://odahub.io/ontology#LongString
number = "ATel #16672"  # http://odahub.io/ontology#String

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in ["text", "number"]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

text_id_text = text
text_id = number
### Settings
try:
    data_path = os.path.dirname(__file__) + "/data/"
except:
    data_path = "data/"
print(data_path)

### Required files
telescope_ontology = f"{data_path}/telescope_observatory_survey.ttl"
file_dict_uri_int = f"{data_path}/dictionary_telescope_uri_int_id.json"
file_dict_sens_inst = (
    f"{data_path}/dictionary_telescope_type_2_instrument.json"
)

simbad_node_file = f"{data_path}/simbad_otypes_nodes.csv"

### Run pipeline
df_tel = rule_based_telescope_detector(
    text_id, text_id_text, telescope_ontology, file_dict_uri_int
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

df_cla = detect_source_classes(
    text_id, text_id_text.lower(), df_sor, simbad_node_file
)

df_vec = vectorize_text(text_id, data_path, df_tel, df_sor, df_astrobert)

df_vec_init_pred = predict_vector(text_id, data_path, df_vec)
df_vec_url, df_url_scores = create_url_vector(
    text_id, data_path, file_dict_sens_inst, df_vec_init_pred, df_sor
)

t_tel = ODAAstropyTable(Table.from_pandas(df_tel.map(str)))
t_sor = ODAAstropyTable(Table.from_pandas(df_sor.map(str)))
t_unk_sor = ODAAstropyTable(Table.from_pandas(df_unk_sor.map(str)))
t_cla = ODAAstropyTable(Table.from_pandas(df_cla.map(str)))
t_astrobert = ODAAstropyTable(Table.from_pandas(df_astrobert.map(str)))
t_vec_init_pred = ODAAstropyTable(Table.from_pandas(df_vec_init_pred.map(str)))
t_vec_url = ODAAstropyTable(Table.from_pandas(df_vec_url.map(str)))
t_url_scores = ODAAstropyTable(Table.from_pandas(df_url_scores.map(str)))

print(df_sor)
print(df_unk_sor)
print(df_tel)
print(df_cla)
print(df_astrobert)
print(df_vec_init_pred)

table_telescopes = t_tel  # http://odahub.io/ontology#ODAAstropyTable
table_sources = t_sor  # http://odahub.io/ontology#ODAAstropyTable
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
