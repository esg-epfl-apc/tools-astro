import re
import pandas as pd
import numpy as np
from transformers import AutoModelForTokenClassification, AutoTokenizer
from transformers import TokenClassificationPipeline


def split_text_in_phrases(text_id, text_):
    list_proto_phrases= re.split(r"(\. [A-Z])", text_)
    for i in range(1, len(list_proto_phrases) - 1, 2):
        back_  = list_proto_phrases[i][0]
        front_ = list_proto_phrases[i][-1]
        list_proto_phrases[i+1] = front_ + list_proto_phrases[i+1]
        list_proto_phrases[i-1] = list_proto_phrases[i-1] + back_

    list_phrases = []
    for i in range(0, len(list_proto_phrases), 2):
        list_phrases.append(list_proto_phrases[i])

    text_check = " ".join(list_phrases)
    if text_check != text_:
        print(text_id)
    return list_phrases


def apply_astroBERT(text_id, body_text_0):

    # load astroBERT for NER-DEAL
    remote_model_path = 'adsabs/astroBERT'
    # you need to load the astroBERT trained for NER-DEAL, which is on a seperate branch
    revision = 'NER-DEAL'

    astroBERT_NER_DEAL = AutoModelForTokenClassification.from_pretrained(
        pretrained_model_name_or_path=remote_model_path,
        revision=revision,
    )

    astroBERT_tokenizer = AutoTokenizer.from_pretrained(
        pretrained_model_name_or_path=remote_model_path,
        add_special_tokens=True,
        do_lower_case=False,
        model_max_length = 512,
    )

    # use the Hugginface Pipeline class
    NER_pipeline=TokenClassificationPipeline(
        model=astroBERT_NER_DEAL,
        tokenizer=astroBERT_tokenizer,
        task='astroBERT NER_DEAL',
        aggregation_strategy='average',
        ignore_labels=['O']
    )

    dict_out = {"TEXT_ID": [], "word":[], "start":[], "end":[], "score":[], "entity_group":[], "Phrase":[]}

    text = " ".join(body_text_0.split()).replace("Â°", "o").replace("Âº", "o").replace("−","-").replace('°', "o")
    list_phrases = split_text_in_phrases(text_id, text)

    for phrase_ in list_phrases:
        result = NER_pipeline(phrase_)

        for u in result:
            ent_ = u["entity_group"]
            if ent_ in ["Instrument", "Telescope", "Wavelength", "CelestialObject", "CelestialRegion", "EntityOfFutureInterest", "Mission", "Observatory", "Survey"]:
                dict_out["TEXT_ID"].append(text_id)
                dict_out["Phrase"].append(phrase_)

                dict_out["word"].append(u["word"])
                dict_out["score"].append(u["score"])
                dict_out["start"].append(u["start"])
                dict_out["end"].append(u["end"])
                dict_out["entity_group"].append(ent_)
    
    return pd.DataFrame(dict_out)


def get_astroBERT_cleaned_result(text_id, body_text_0):
    list_entities = ["Instrument", "Telescope", "Wavelength", "CelestialObject", "CelestialRegion", "EntityOfFutureInterest", "Mission", "Observatory", "Survey"]

    df_raw = apply_astroBERT(text_id, body_text_0)
    dict_out = {"TEXT_ID": [], "word":[], "start":[], "end":[], "Score":[], "Phrase":[], "entity_group":[]}

    for entity_to_study in list_entities:
        df_tmp0 = df_raw[df_raw["entity_group"] == entity_to_study] 
        phrases_ = np.unique(df_tmp0["Phrase"])

        for phrase_ in phrases_:
            df_tmp1 = df_tmp0[df_tmp0["Phrase"]==phrase_]
            if len(df_tmp1) == 1:
                dict_out["TEXT_ID"].append(text_id)
                dict_out["Phrase"].append(df_tmp1.Phrase.values[0])
                dict_out["word"].append(df_tmp1.word.values[0])
                dict_out["start"].append(df_tmp1.start.values[0])
                dict_out["end"].append(df_tmp1.end.values[0])
                dict_out["Score"].append(df_tmp1.score.values[0])
                dict_out["entity_group"].append(entity_to_study)

            else:
                df_tmp1.sort_values(by=['start'])
                for s_i, (s_, e_, sc_) in enumerate(zip(df_tmp1.start.values, df_tmp1.end.values, df_tmp1.score.values)):
                    if s_i == 0:
                        s_o = s_
                        e_o = e_
                        sc_o = sc_
                        sc_s = sc_
                        word_size = 1
                    else:

                        if s_ <= e_o + 1:
                            e_o = e_
                            sc_o = sc_
                            sc_s += sc_
                            word_size += 1

                        else:
                            dict_out["TEXT_ID"].append(text_id)
                            dict_out["Phrase"].append(phrase_)
                            dict_out["word"].append(phrase_[s_o: e_o])
                            dict_out["start"].append(s_o)
                            dict_out["end"].append(e_o)
                            dict_out["Score"].append(sc_s / word_size)
                            dict_out["entity_group"].append(entity_to_study)
                            
                            s_o = s_
                            e_o = e_
                            sc_o = sc_
                            sc_s = sc_
                            word_size = 1

                    if s_i == len(df_tmp1) - 1:
                        dict_out["TEXT_ID"].append(text_id)
                        dict_out["Phrase"].append(phrase_)
                        dict_out["word"].append(phrase_[s_o: e_o])
                        dict_out["start"].append(s_o)
                        dict_out["end"].append(e_o)
                        dict_out["Score"].append(sc_s / word_size)
                        dict_out["entity_group"].append(entity_to_study)

    return pd.DataFrame(dict_out)


def clean_ra_dec