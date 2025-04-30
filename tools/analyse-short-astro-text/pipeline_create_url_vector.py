import pandas as pd
import numpy as np
import json
import urllib

from pipeline_vectorize_text import otype_to_index
from aux_functions import get_dict_instruments_URL_MMODA


def get_link(ra, dec, T1, T2, instrument, src_name):
    params = {
        # "T1": T1,
        # "T2": T2,
        # "T_format": "isot",
        "instrument": instrument
    }

    if src_name is not None and src_name != "":
        params["src_name"] = src_name

    if ra is not None:
        params["RA"] = ra

    if dec is not None:
        params["DEC"] = dec

    rest_url = urllib.parse.urlencode(params)
    return f"https://www.astro.unige.ch/mmoda/?{rest_url}"


def source_to_source_type_index(df_dict, df_sor):
    if df_sor.empty:
        return {"": {"Indices": [], "RA": None, "Dec": None}}

    else:
        dict_out = {}
        for main_id, input_otypes, ra_, dec_ in zip(df_sor["Main ID Name"].values, df_sor["OTYPE"].values, df_sor["RA"].values, df_sor["Dec"].values):
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


def create_url_vector(text_id, data_path, file_dict_sens_inst, df_vec_init_pred, df_sor):
    otype_label = f"{data_path}/dict_source_otypes_considered_for_prediction.csv"
    df_dict = pd.read_csv(otype_label)

    dict_source_to_type_indx = source_to_source_type_index(df_dict, df_sor)
    inst_2_inst_name = get_dict_instruments_URL_MMODA()

    with open(file_dict_sens_inst, "r") as fp:
        dict_sens_inst = json.load(fp)
    pred = df_vec_init_pred["Follow-up Vector Prediction"].values
    legend = df_vec_init_pred["Legend"].values

    dict_ = {}
    dict_["Legend"] = legend
    dict_[text_id] = df_vec_init_pred[text_id].values
    dict_["Follow-up"] = pred

    pred_norm_vector = pred / np.linalg.norm(pred)
    dict_url_scores = {"URL Name": [], "Scores": [], "URL": []}

    for source_name in dict_source_to_type_indx.keys():
        counter = 0

        # Loop through the telescope types: gamma-ray, x-ray, etc.
        for i, (value, name) in enumerate(zip(pred, legend)):
            if i >= 41 and i < 50 and name in dict_sens_inst.keys():

                # Loop through the instruments corresponding to a telescope type: e.g.for gamma-ray, there are spi_acs, polar, grb_detection, magic, etc (see the file file_dict_sens_inst)
                for inst_ in dict_sens_inst[legend[i]]:

                    # if the instrument is explicitly part of the input/output URL vector
                    if inst_ in inst_2_inst_name.keys():
                        url_vec_telescope_telescope_type = np.zeros(len(legend))
                        url_vec_telescope_telescope_type[i] = 1

                        for indx_ in dict_source_to_type_indx[source_name]["Indices"]:
                            url_vec_telescope_telescope_type[indx_] = 1

                        for inst_name in inst_2_inst_name[inst_]:
                            inst_indx = np.where(legend == inst_name)[0].squeeze()
                            url_vec_telescope_telescope_type[inst_indx] = 1

                        dict_[f"URL_{counter}{source_name}"] = url_vec_telescope_telescope_type

                        url_norm_vector = (url_vec_telescope_telescope_type)/np.sum(url_vec_telescope_telescope_type**2)
                        score = np.round(np.dot(url_norm_vector, pred_norm_vector), decimals=5)
                        dict_url_scores["URL Name"].append(f"URL_{counter}{source_name}")
                        dict_url_scores["Scores"].append(score)
                        dict_url_scores["URL"].append(get_link(dict_source_to_type_indx[source_name]["RA"], dict_source_to_type_indx[source_name]["Dec"], "AAA", "BBB", inst_, source_name))

                        counter += 1

        # add only the telescope type in order to represent the instruments like Polar, CTA/CTAO, and workflows LegacySurvey DESI, SGWB that are not part of the size 59 vector
        for indx_tel_type in [41, 44, 45, 48]:
            if pred[indx_tel_type] != 0:
                url_vec_telescope_telescope_type = np.zeros(len(legend))
                for indx_ in dict_source_to_type_indx[source_name]["Indices"]:
                    url_vec_telescope_telescope_type[indx_] = 1

                url_vec_telescope_telescope_type[indx_tel_type] = 1

                dict_[f"URL_{counter}{source_name}"] = url_vec_telescope_telescope_type

                url_norm_vector = (url_vec_telescope_telescope_type)/np.sum(url_vec_telescope_telescope_type**2)
                score = np.round(np.dot(url_norm_vector, pred_norm_vector), decimals=5)

                for inst_ in dict_sens_inst[legend[indx_tel_type]]:
                    if inst_ not in inst_2_inst_name.keys():
                        dict_url_scores["URL Name"].append(f"URL_{counter}{source_name}")
                        dict_url_scores["Scores"].append(score)
                        dict_url_scores["URL"].append(get_link(dict_source_to_type_indx[source_name]["RA"], dict_source_to_type_indx[source_name]["Dec"], "AAA", "BBB", inst_, source_name))

                counter += 1

    df_out = pd.DataFrame(dict_)

    return df_out, pd.DataFrame(dict_url_scores).sort_values("Scores", ascending=False)
