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
from pipeline_predict_sensitivity import predict_sensitivity
from pipeline_provide_workflows import provide_workflows
from pipeline_source_classes import detect_source_classes
from pipeline_sources import rule_based_source_detector
from pipeline_telescope import rule_based_telescope_detector

text = "We received an alert about a candidate fast transient from the FINK broker (Peloton, Ishida & Moller et al.) on 2024-06-22 09:58:42.001 UTC, named ZTF24aasjjkf (https://fink-portal.org/ZTF24aasjjkf). This transient was classified by FINK on the same day with a 17 percent probability of being a fast transient. We triggered the TAROT-TCA and TRT-SRO telescopes for further investigation, observing in the R, r', and i' bands, and collecting data from 2024-06-23 around 04:55:14 UTC to 2024-06-24 01:00:00 UTC. Our results, along with public ZTF data, can be found here: https://skyportal-icare.ijclab.in2p3.fr/public/sources/ZTF24aasjjkf/version/07cf393bb53513abe80d3b7863938d8f TAROT-TCA telescope observations can be also labeled as generic in the plot. We utilized STDWeb (Karpov et al.) to perform our photometry, using the Pan-STARRS DR1 (PS1) catalog and subtracting the constant flux of background with a reference PS1 image. We note that in the TAROT-TCA data, we have slight contamination from other stars. However, we mitigated this effect by subtracting the PS1 reference image. According to the ZTF data in the r' band, both before and after our observations, as well as ATLAS data (in the orange band) and the Gaia Alert (2024LXA), the source might be 10 days old or more. Its multi-band light curve is not consistent with a fast transient like a kilonova (a decay rate of 0.15 per day in r-band) but resembles a supernova (also classified by FINK on 2024-06-24 at 08:00 UTC). The source has been independently classified by Jianlin Xu et al. as a CV. We thank the FINK team for their valuable collaboration with GRANDMA. GRANDMA is a worldwide telescope network (https://grandma.lal.in2p3.fr/) devoted to the observation of transients in the context of multi-messenger astrophysics (Antier et al. 2020 MNRAS 497, 5518). Kilonova-Catcher (KNC) is the citizen science program of GRANDMA (http://kilonovacatcher.in2p3.fr/)."  # http://odahub.io/ontology#LongString
number = "ATel #16672"  # http://odahub.io/ontology#LongString

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
model_file = f"{data_path}/probabilistic_models/first_single_second_single_with_tel_repetition.dat.npy"

### Required files
telescope_ontology = f"{data_path}/telescope_observatory_survey.ttl"
file_dict_uri_int = f"{data_path}/dictionary_telescope_uri_int_id.json"
file_dict_sens_inst = (
    f"{data_path}/dictionary_telescope_type_2_instrument.json"
)

dict_path = f"{data_path}/dictionary_names_otypes.json"
simbad_node_file = f"{data_path}/simbad_otypes_nodes.csv"

### Run pipeline
df_tel = rule_based_telescope_detector(
    atel_, atel_text, telescope_ontology, file_dict_uri_int
)
df_sor = rule_based_source_detector(atel_, atel_text, data_path, dict_path)
df_cla = detect_source_classes(
    atel_, atel_text.lower(), df_sor, simbad_node_file
)
df_astrobert = get_astroBERT_cleaned_result(atel_, atel_text)

df_in, df_pred = predict_sensitivity(model_file, df_tel, first_type="single")
df_workflows = provide_workflows(df_in, df_pred, df_sor, file_dict_sens_inst)

t_tel = ODAAstropyTable(Table.from_pandas(df_tel.map(str)))
t_sor = ODAAstropyTable(Table.from_pandas(df_sor.map(str)))
t_cla = ODAAstropyTable(Table.from_pandas(df_cla.map(str)))
t_astrobert = ODAAstropyTable(Table.from_pandas(df_astrobert.map(str)))
t_in = ODAAstropyTable(Table.from_pandas(df_in.map(str)))
t_pred = ODAAstropyTable(Table.from_pandas(df_pred.map(str)))
t_workflows = ODAAstropyTable(Table.from_pandas(df_workflows.map(str)))

print(t_sor)
print(t_tel)
print(t_cla)
print(t_in)
print(t_pred)
print(t_workflows)

table_telescopes = t_tel  # http://odahub.io/ontology#ODAAstropyTable
table_sources = t_sor  # http://odahub.io/ontology#ODAAstropyTable
table_source_classes = t_cla  # http://odahub.io/ontology#ODAAstropyTable
table_astrobert_results = (
    t_astrobert  # http://odahub.io/ontology#ODAAstropyTable
)
table_atel_sensitivity = t_in  # http://odahub.io/ontology#ODAAstropyTable
table_predicted_sensitivity = (
    t_pred  # http://odahub.io/ontology#ODAAstropyTable
)
table_suggested_workflows = (
    t_workflows  # http://odahub.io/ontology#ODAAstropyTable
)

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_telescope_type_prediction_table_telescopes",
        "table_telescopes_galaxy.output",
        table_telescopes,
    )
)
_oda_outs.append(
    (
        "out_telescope_type_prediction_table_sources",
        "table_sources_galaxy.output",
        table_sources,
    )
)
_oda_outs.append(
    (
        "out_telescope_type_prediction_table_source_classes",
        "table_source_classes_galaxy.output",
        table_source_classes,
    )
)
_oda_outs.append(
    (
        "out_telescope_type_prediction_table_astrobert_results",
        "table_astrobert_results_galaxy.output",
        table_astrobert_results,
    )
)
_oda_outs.append(
    (
        "out_telescope_type_prediction_table_atel_sensitivity",
        "table_atel_sensitivity_galaxy.output",
        table_atel_sensitivity,
    )
)
_oda_outs.append(
    (
        "out_telescope_type_prediction_table_predicted_sensitivity",
        "table_predicted_sensitivity_galaxy.output",
        table_predicted_sensitivity,
    )
)
_oda_outs.append(
    (
        "out_telescope_type_prediction_table_suggested_workflows",
        "table_suggested_workflows_galaxy.output",
        table_suggested_workflows,
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
