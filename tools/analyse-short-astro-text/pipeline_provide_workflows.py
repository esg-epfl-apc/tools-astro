import numpy as np
import pandas as pd
import os
import json
import glob
import urllib.parse

from aux_functions import find_sensitivity, bits2sens


def get_src_names(df_sor):
    dict_out = {}
    
    if not df_sor.empty:
        for raw_name in df_sor["Raw Source Name"]:
            main_name = df_sor["Main ID Name"][df_sor["Raw Source Name"] == raw_name].values[0]
            ra  = df_sor["RA"][df_sor["Raw Source Name"] == raw_name].values[0]
            dec = df_sor["Dec"][df_sor["Raw Source Name"] == raw_name].values[0]

            if main_name == "NotKnown":
                src_name = raw_name
            else:
                src_name = main_name

            if ra == "NotKnown":
                ra = None

            if dec == "NotKnown":
                dec = None

            dict_out[src_name] =  {"RA": ra, "Dec": dec}

    if len(dict_out.keys()) == 0:
        dict_out[None] = {"RA": None, "Dec": None}

    return dict_out


def provide_workflows(df_in, df_pred, df_sor, file_dict_sens_inst):
    with open(file_dict_sens_inst, "r") as fp:
        dict_sens_inst = json.load(fp)
    sens_ = list(set(list(df_in["Initial Sensitivity"].values) + list(df_pred["Predicted Sensitivity"].values)))
    
    dict_out = {"Sensitivity" :[], "Instrument": [], "Workflow": []}
    
    T1  = "CCC"
    T2  = "DDD"

    dict_src_names = get_src_names(df_sor)

    for s_ in sens_:
        if s_ in dict_sens_inst.keys():
            for inst_ in dict_sens_inst[s_]:
                for src_name in dict_src_names.keys():
                    dict_out["Sensitivity"].append(s_)
                    dict_out["Instrument"].append(inst_)
                    dict_out["Workflow"].append(get_link(dict_src_names[src_name]["RA"], dict_src_names[src_name]["Dec"], T1, T2, inst_, src_name))
        
    return pd.DataFrame(dict_out)


def get_link(ra, dec, T1, T2, instrument, src_name):
    params = {
        # "T1": T1,
        # "T2": T2,
        # "T_format": "isot",
        "instrument": instrument
    }

    if src_name != None:
        params["src_name"] = src_name

    if ra != None:
        params["RA"] = ra
    
    if dec != None:
        params["DEC"] = dec

    rest_url = urllib.parse.urlencode(params)
    return f"https://www.astro.unige.ch/mmoda/?{rest_url}"
    