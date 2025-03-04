import pandas as pd
import numpy as np
import json
import re
import sys

from rdflib import Graph, Literal, URIRef, Namespace
from rdflib.namespace import SKOS, DCTERMS, RDF, RDFS

from aux_functions import list_tel, compute_sensitivity, compute_sensitivity_int

ODA = Namespace("https://odahub.io/ontology#")
g_label_site = len("http://www.w3.org/2000/01/rdf-schema#label")


def find_entity(g, class_, text_id_text, text_id_text_upper):
    dict_ = {"label":{"val": [], "URI": [], "Sensitivity": []}, "altLabel":{"val": [], "URI": [], "Sensitivity": []}, "hiddenLabel":{"val": [], "URI": [], "Sensitivity": []}}
    for i, [u_telescope, p, o] in enumerate(g.triples((None, RDF.type, class_))):
            
        exists_label_telescope = 0
        for s, p, label_telescope in g.triples((u_telescope, RDFS.label, None)):

            result = re.search("\\b(" + label_telescope.lower() + ")([1-2]{0,1})\\b", text_id_text)
            if result:
                exists_label_telescope = 1
                val_ = result.group(0)
                add_ = val_[len(label_telescope):]
                dict_["label"]["val"].append(label_telescope)
                dict_["label"]["URI"].append(u_telescope)
                dict_["label"]["Sensitivity"].append(compute_sensitivity(list_tel(u_telescope, g)))

        exists_altlabel_telescope = 0
        if exists_label_telescope == 0:
            for s, p, altlabel_telescope in g.triples((u_telescope, SKOS.altLabel, None)):

                result = re.search("\\b(" + altlabel_telescope + ")\\b", text_id_text_upper)
                if result:
                    exists_altlabel_telescope = 1
                    val_ = result.group(0)
                    dict_["altLabel"]["val"].append(val_)
                    dict_["altLabel"]["URI"].append(u_telescope)
                    dict_["altLabel"]["Sensitivity"].append(compute_sensitivity(list_tel(u_telescope, g)))

            if exists_altlabel_telescope == 0:
                for s, p, hiddenlabel_telescope in g.triples((u_telescope, SKOS.hiddenLabel, None)):

                    result = re.search("\\b" + hiddenlabel_telescope.lower() + "\\b", text_id_text)
                    if result:
                        val_ = result.group(0)
                        dict_["hiddenLabel"]["val"].append(hiddenlabel_telescope)
                        dict_["hiddenLabel"]["URI"].append(u_telescope)
                        dict_["hiddenLabel"]["Sensitivity"].append(compute_sensitivity(list_tel(u_telescope, g)))
                
    return dict_


def rule_based_telescope_detector(text_id, text_id_text, telescope_ontology, file_dict_uri_int):
    g = Graph()
    g.parse(telescope_ontology, format="n3")
    
    with open(file_dict_uri_int, "r") as fp:
        dict_uri_int = json.load(fp)
    
    
    text_id_text_lower = text_id_text.lower()

    dict_observatory         = find_entity(g, ODA.observatory,         text_id_text_lower, text_id_text)
    dict_survey              = find_entity(g, ODA.survey,              text_id_text_lower, text_id_text)
    dict_telescope           = find_entity(g, ODA.telescope,           text_id_text_lower, text_id_text)
    dict_misctelescope       = find_entity(g, ODA.misctelescope,       text_id_text_lower, text_id_text)
    dict_telescopetype       = find_entity(g, ODA.telescopetype,       text_id_text_lower, text_id_text)

    dict_spacetelescope      = find_entity(g, ODA.spacetelescope,      text_id_text_lower, text_id_text)
    dict_instrument          = find_entity(g, ODA.instrument,          text_id_text_lower, text_id_text)
    dict_institution         = find_entity(g, ODA.institution,         text_id_text_lower, text_id_text)
    dict_radiotelescope      = find_entity(g, ODA.radiotelescope,      text_id_text_lower, text_id_text)

    tel_sur_obs = []
    type_key = []
    uri_list = []
    sens_list = []

    for key in ["label", "altLabel", "hiddenLabel"]:
        list_key = dict_institution[key]["val"] + dict_spacetelescope[key]["val"] + dict_telescope[key]["val"] + dict_survey[key]["val"] + dict_observatory[key]["val"] + dict_radiotelescope[key]["val"] + dict_instrument[key]["val"] + dict_telescopetype[key]["val"] + dict_misctelescope[key]["val"]
        tel_sur_obs += list_key

        list_uri_key = dict_institution[key]["URI"] + dict_spacetelescope[key]["URI"] + dict_telescope[key]["URI"] + dict_survey[key]["URI"] + dict_observatory[key]["URI"] + dict_radiotelescope[key]["URI"] + dict_instrument[key]["URI"] + dict_telescopetype[key]["URI"] + dict_misctelescope[key]["URI"]
        uri_list += list_uri_key
        
        sens_list += dict_institution[key]["Sensitivity"] + dict_spacetelescope[key]["Sensitivity"] + dict_telescope[key]["Sensitivity"] + dict_survey[key]["Sensitivity"] + dict_observatory[key]["Sensitivity"] + dict_radiotelescope[key]["Sensitivity"] + dict_instrument[key]["Sensitivity"] + dict_telescopetype[key]["Sensitivity"] + dict_misctelescope[key]["Sensitivity"]

        type_key += [key]*len(list_key)
    
    dict_data = {"TEXT_ID": [text_id] * len(tel_sur_obs), "Telescope": tel_sur_obs, "LabelType": type_key, "URI": uri_list, "Sensitivity": sens_list, "Total Sensitivity": [compute_sensitivity_int(sens_list)] * len(tel_sur_obs)}
    
    df_data = pd.DataFrame(dict_data)
    df_data.drop_duplicates(subset=['URI'], inplace=True)

    return df_data
