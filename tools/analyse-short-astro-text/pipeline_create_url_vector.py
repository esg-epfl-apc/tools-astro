import pandas as pd
import numpy as np
import json
import re
import sys
import urllib

from rdflib import Graph, Literal, URIRef, Namespace
from rdflib.namespace import SKOS, DCTERMS, RDF, RDFS

from aux_functions import list_tel, compute_sensitivity, compute_sensitivity_int
from pipeline_vectorize_text import otype_to_index


def get_link(ra, dec, T1, T2, instrument, src_name):
    params = {
        # "T1": T1,
        # "T2": T2,
        # "T_format": "isot",
        "instrument": instrument
    }

    if src_name != None and src_name != "":
        params["src_name"] = src_name

    if ra != None:
        params["RA"] = ra
    
    if dec != None:
        params["DEC"] = dec

    rest_url = urllib.parse.urlencode(params)
    return f"https://www.astro.unige.ch/mmoda/?{rest_url}"
    


def source_to_source_type_index(df_dict, df_astrobert_source_types, df_sor):
    df_sor_n = pd.DataFrame()
    if not df_astrobert_source_types.empty and not df_sor_n.empty:
        df_sor_n = pd.concat((df_sor[["Main ID Name", "OTYPE", "RA", "Dec"]], df_astrobert_source_types[["Main ID Name", "OTYPE", "RA", "Dec"]]))
    
    elif not df_astrobert_source_types.empty:
        df_sor_n = df_astrobert_source_types
        
    elif not df_sor.empty:
        df_sor_n = df_sor
    
    if df_sor_n.empty:
        return {"": {"Indices": [], "RA": None, "Dec": None}}
    
    else:
        dict_out = {}
        for main_id, input_otypes, ra_, dec_ in zip(df_sor_n["Main ID Name"].values, df_sor_n["OTYPE"].values, df_sor_n["RA"].values, df_sor_n["Dec"].values):
            output_otype = []
            if main_id != "NotKnown" and input_otypes != "NotKnown":
                for input_otype in set(input_otypes.split("|")):
                    output_otype.append(input_otype)

                output_otype = set(output_otype)
                dict_out[main_id] = {"Indices": otype_to_index(output_otype, df_dict), "RA": ra_, "Dec": dec_}

        if len(dict_out.keys()) != 0:
            return dict_out
        else:
            return {"": {"Indices": [], "RA": None, "Dec": None}}


def create_url_vector(atel_, data_path, file_dict_sens_inst, df_vec_init_pred, df_astrobert_source_types, df_sor):
    otype_label = f"{data_path}/dict_source_otypes_considered_for_prediction.csv"
    df_dict = pd.read_csv(otype_label)
    
    dict_source_to_type_indx = source_to_source_type_index(df_dict, df_astrobert_source_types, df_sor)
    inst_2_inst_name = {
        "spi_acs": ["INTEGRAL", "SPI-ACS"],
        "cta"    : ["CTA/CTAO"],
        "hess"   : ["HESS"],
        "isgri"  : ["INTEGRAL", "ISGRI"],
        "jemx"   : ["INTEGRAL", "JEM-X"],
        "icecube": ["IceCube"],
        "antares": ["ANTARES"],
        "gw"     : ["LIGO/VIRGO"]
    }
    
    with open(file_dict_sens_inst, "r") as fp:
        dict_sens_inst = json.load(fp)
    pred = df_vec_init_pred["Follow-up Vector Prediction"].values
    legend = df_vec_init_pred["Legend"].values
    
    dict_ = {}
    dict_["Legend"] = legend
    dict_[atel_] = df_vec_init_pred[atel_].values
    dict_["Follow-up"] = pred

    pred_norm_vector = pred / np.sum(pred**2)
    dict_url_scores = {"URL Name": [], "Scores": [], "URL": []}

    for source_name in dict_source_to_type_indx.keys():
        counter = 0
        for i, (value, name) in enumerate(zip(pred, legend)):
            if i >= 41 and i < 50 and name in dict_sens_inst.keys():
                for inst_ in dict_sens_inst[legend[i]]:
                    if inst_ in inst_2_inst_name.keys():
                        url_vec_telescope_telescope_type = np.zeros(59)
                        url_vec_telescope_telescope_type[i] = 1
                        
                        for indx_ in dict_source_to_type_indx[source_name]["Indices"]:
                            url_vec_telescope_telescope_type[indx_] = 1

                        for inst_name in inst_2_inst_name[inst_]:
                            inst_indx = np.where(legend==inst_name)[0].squeeze()
                            url_vec_telescope_telescope_type[inst_indx] = 1

                        dict_[f"URL_{counter}{source_name}"] = url_vec_telescope_telescope_type
                        
                        url_norm_vector = (url_vec_telescope_telescope_type)/np.sum(url_vec_telescope_telescope_type**2)
                        score =  np.dot(url_norm_vector, pred_norm_vector)                        
                        dict_url_scores["URL Name"].append(f"URL_{counter}{source_name}")
                        dict_url_scores["Scores"].append(score)
                        dict_url_scores["URL"].append(get_link(dict_source_to_type_indx[source_name]["RA"], dict_source_to_type_indx[source_name]["Dec"], "AAA", "BBB", inst_, source_name))
                     
                        counter += 1

        ### add only the telescope type in order to represent the instruments like Polar, CTA/CTAO, and workflows LegacySurvey DESI, SGWB that are not part of the size 59 vector
        for indx_tel_type in [41, 44, 45, 48]:
            if pred[indx_tel_type] != 0:
                url_vec_telescope_telescope_type = np.zeros(59)
                for indx_ in dict_source_to_type_indx[source_name]["Indices"]:
                    url_vec_telescope_telescope_type[indx_] = 1

                url_vec_telescope_telescope_type[indx_tel_type] = 1

                dict_[f"URL_{counter}{source_name}"] = url_vec_telescope_telescope_type
                
                url_norm_vector = (url_vec_telescope_telescope_type)/np.sum(url_vec_telescope_telescope_type**2)
                score =  np.dot(url_norm_vector, pred_norm_vector)                        
                
                
                for inst_ in dict_sens_inst[legend[indx_tel_type]]:
                    if not inst_ in inst_2_inst_name.keys(): 
                        dict_url_scores["URL Name"].append(f"URL_{counter}{source_name}")
                        dict_url_scores["Scores"].append(score)
                        dict_url_scores["URL"].append(get_link(dict_source_to_type_indx[source_name]["RA"], dict_source_to_type_indx[source_name]["Dec"], "AAA", "BBB", inst_, source_name))
                        
                counter += 1

    df_out = pd.DataFrame(dict_)

    return df_out, pd.DataFrame(dict_url_scores).sort_values("Scores", ascending=False)
