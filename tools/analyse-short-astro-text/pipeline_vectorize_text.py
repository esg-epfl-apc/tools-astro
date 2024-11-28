import pandas as pd
import numpy as np
import rdflib
import os
import json
from aux_functions import bits, compute_sensitivity, find_sensitivity, list_tel, find_name_workflow
import time

from pipeline_sources import query_simbad


import warnings
warnings.filterwarnings("ignore", 'This pattern is interpreted as a regular expression, and has match groups')



######## Vectorize
# From 0  to 40 (41 indicies) represent the source types

# From 41 (41 + log2(bit)) to 49 represent the telescope type

# From 50 to 58 represent the workflow (telescope/instrument)


def otype_to_index(list_otype, df_dict):
    out_indices = []

    for t_o in list_otype:
        df_tmp = df_dict[df_dict["All types"]==t_o]
        if len(df_tmp) == 0:
            if t_o == "SN I" or t_o == "SN Ia" or t_o == "SN Ib" or t_o == "SN Ib/c" or t_o == "SN Ic" or t_o == "SN Ic-BL" or t_o == "SN II" or t_o == "SN IIP" or t_o == "SN Ia-91T-like" or t_o == "SN IIn" or t_o == "SN IIb" or t_o == "SN Ia-91bg-like":
                df_tmp = df_dict[df_dict["All types"]=="SN*"]
            elif t_o == "Other":
                df_tmp = df_dict[df_dict["All types"]=="?"]
            elif t_o == "Varstar":
                df_tmp = df_dict[df_dict["All types"]=="V*"]
            elif t_o == "CV":
                df_tmp = df_dict[df_dict["All types"]=="CV*"]
            elif t_o == "M dwarf" or t_o == "C?*":
                df_tmp = df_dict[df_dict["All types"]=="*"]
            elif t_o == "Galaxy":
                df_tmp = df_dict[df_dict["All types"]=="G"]
                
            
        if len(df_tmp) != 0:
            out_indices.append(df_tmp["Index in vector"].values[0])
        else:
            print(t_o)
        
    return out_indices


def get_source_indices(df_celestialobj, df_dict):
    output_otype = []
    for input_otypes in df_celestialobj["OTYPE"].values:
        if not pd.isnull(input_otypes):
            for input_otype in set(input_otypes.split("|")):
                output_otype.append(input_otype)
                
    output_otype = set(output_otype)
    return otype_to_index(output_otype, df_dict)


def vectorize_text(text_id, data_path, df_regex_telescopes, df_celestialobj, df_astrobert):

    otype_label = f"{data_path}/dict_source_otypes_considered_for_prediction.csv"
    df_dict = pd.read_csv(otype_label)
    dict_out = {}

    dict_out["Legend"] = list(df_dict[["Representatives", "Index in vector"]].drop_duplicates().sort_values('Index in vector').Representatives.values) + ["gamma-ray", "x-ray", "ultraviolet", "optical", "infrared", "radio", "cosmic-ray", "gravitational-wave", "neutrino"] + ["INTEGRAL", "ISGRI", "JEM-X", "SPI-ACS", "ANTARES", "LIGO/VIRGO", "IceCube", "HESS", "CTA/CTAO"]
    data_vector = np.zeros(59, dtype='int32')


    ##### REGEX
    ## Telescope Type
    if not df_regex_telescopes.empty:
        _, list_bits_, _ = find_sensitivity(df_regex_telescopes["Total Sensitivity"].values[0])
        for bit_ in list_bits_:
            indx_2 = int(np.log2(bit_)) + 41
            data_vector[indx_2] += 1

        ## Workflow
        for indx_3 in find_name_workflow(df_regex_telescopes.URI.values):
            data_vector[indx_3] += 1
   
    ### astroBERT + REGEX Source Type
    if not df_celestialobj.empty:
        out_indices_ = get_source_indices(df_celestialobj, df_dict)
        for indx_1 in out_indices_:
            data_vector[indx_1] += 1


    ##### astroBERT
    if not df_astrobert.empty:

        ## Telescope Type
        df_tmp = df_astrobert[df_astrobert["entity_group"] == "Wavelength"].dropna()
        if not df_tmp.empty:
            patterns = ["(x(-| )ray)", "(gamma(-| )ray)", "(optical)", "(radio)", "(ultra(-|)violet)", "(infra(-|)red)"]
            names    = ["x-ray",       "gamma-ray",       "optical",   "radio",   "ultraviolet",       "infrared"]
            for pattern_, name_ in zip(patterns, names):
                indx_2 = int(np.log2(bits(ask=name_))) + 41
                data_vector[indx_2] += len(df_tmp[df_tmp.word.str.contains(pattern_, case=False, regex=True)])

            patterns = ["IR",       "UV"]
            names    = ["infrared", "ultraviolet"]
            for pattern_, name_ in zip(patterns, names):
                indx_2 = int(np.log2(bits(ask=name_))) + 41
                data_vector[indx_2] += len(df_tmp[df_tmp.word.str.contains(pattern_, case=True, regex=True)])


    dict_out[f"{text_id}"] = data_vector

    return pd.DataFrame(dict_out)