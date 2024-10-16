import numpy as np
import pandas as pd
import os
import glob

from aux_functions import find_sensitivity, bits2sens


def prepare_sensitivity_vector(model, sensitivity, type_=None):    
    sens_vec = np.zeros(model.shape[0])
    
    if type_ == "all":
        sens_vec[sensitivity] = 1
    
    elif type_ == "single":
        _, list_b, _ = find_sensitivity(sensitivity)
        for bit_ in list_b:
            sens_vec[bit_] = 1

    elif type_ == "pairoftwo":
        _, _, list_b_12 = find_sensitivity(sensitivity)
        for bit_ in list_b_12:
            sens_vec[bit_] = 1
    # else:
        # print(f"You chose {type_}. This does not exist!")
        # os._exit(0)
    
    sens_vec[0] = 0
    return sens_vec


def predict_sensitivity(model_file, df_tmp, first_type=None, norm_="max"):
    predictions_ = []
    scores_      = []
    input_ = []
    
    if len(df_tmp) > 0:

        first_sens = df_tmp["ATEL Sensitivity"].values[0]

        model = np.load(model_file)
        pred_sens_vec = np.zeros(model.shape[0])

        first_sens_vec  = prepare_sensitivity_vector(model, first_sens, type_=first_type)

        bits_sens_vec = np.argwhere(first_sens_vec)

        if bits_sens_vec.shape[0] == 1:
            bits_sens_vec = [bits_sens_vec.squeeze()]
        else:
            bits_sens_vec = bits_sens_vec.squeeze()

        for bit_ in bits_sens_vec:
            y__ = model[bit_, :]
            y__[0] = 0

            if norm_ == "max":
                norm_fac = np.max(y__)
            elif norm_ == "sum":
                norm_fac = np.sum(y__)
            else:
                # print("wrong norm fac")
                norm_fac = 1

            if norm_fac == 0:
                norm_fac = 1

            pred_sens_vec += (y__ / norm_fac)

        if norm_ == "max":
            norm_fac = np.max(pred_sens_vec)
        elif norm_ == "sum":
            norm_fac = np.sum(pred_sens_vec)
        else:
            # print("wrong norm fac")
            norm_fac = 1

        if norm_fac == 0:
            norm_fac = 1

        normalized_pred_sens_vec = pred_sens_vec / norm_fac

        ### Choose the first 3 predictions and shorten the vector
        arg_pred_ = np.argsort(normalized_pred_sens_vec, axis=0)
        np.put_along_axis(normalized_pred_sens_vec, arg_pred_[:-3], 0, axis=0)

        indx_ = np.argwhere(normalized_pred_sens_vec)
        scores_ = normalized_pred_sens_vec[indx_].squeeze() / np.max(normalized_pred_sens_vec[indx_])

        for bit_ in indx_:
            predictions_.append(bits2sens(bit_))

        ### Clean the input vector
        indx_i = np.argwhere(first_sens_vec)
        for bit_ in indx_i:
            input_.append(bits2sens(bit_))
        
    # input_, 
    dict_in = {"Initial Sensitivity":[]}
    for i, in_ in enumerate(input_):
        dict_in["Initial Sensitivity"].append(in_)
        
    
    dict_out = {"Predicted Sensitivity":[], "Score":[]}
    for i, pred_ in enumerate(predictions_):
        dict_out["Predicted Sensitivity"].append(pred_)
        dict_out["Score"].append(scores_[i])
        
    return pd.DataFrame(dict_in), pd.DataFrame(dict_out)
