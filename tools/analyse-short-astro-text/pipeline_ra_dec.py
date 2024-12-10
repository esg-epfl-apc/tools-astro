from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd
import numpy as np
import json
import re
import sys


def split_text_in_phrases(atel_, text_):
    list_proto_phrases= re.split(r"(\. [a-z])", text_)
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
        print(atel_)

    return list_phrases


def create_pattern_list():
    pattern_list_ra     = []
    pattern_list_dec    = []
    pattern_list_ra_dec = []
    pattern_list_table  = []
    ra_text  = "r(\\.|)a(\\.|\\:|)"
    dec_text = "dec(l|)(\\.|)"


    units_minutes = "(\\'|m|\\:|\\' |m |\\: |)"
    units_minutes_mix_all = "(\\'|m|\\:| )"
    units_seconds = '(\\"|s|\\:|\\" |s |\\: |)'
    units_seconds_mix_all = '(\\"|s|\\:| |)'

    degree_ = "((deg)|d|)"
    assignment_char = "(( |)(=|\\:)( |))"

    units_dec_min = "\\'m"
    units_dec_sec = '\\"s'
    units_dec_deg = "dego"
    ra_value  = f"([0-9\\.\\:\\s{units_dec_min}{units_dec_sec}hdeg]{{2,}})"
    dec_value = f"(\\+|-|)([0-9\\.\\:\\s{units_dec_min}{units_dec_sec}{units_dec_deg}]{{2,}})"    
    ra_value_deg    = "(\\d{1,3}(\\.|)\\d{0,})"
    dec_value_deg   = "(\\+|-|)(\\d{1,2}(\\.|)\\d{0,})"

    ra_dec = f"({ra_text},( |){dec_text})"

    ### GOOD
    pattern_list_ra  += [f"\\b({ra_text}{assignment_char}{ra_value})\\b"]
    pattern_list_dec += [f"\\b({dec_text}{assignment_char}{dec_value})\\b"]
        
    pattern_J2000 = "( |)(-|)((\((((j|)2000(.0|))|(deg))\))|(2000))"
    pattern_list_ra  += [f"({ra_text}{pattern_J2000}{assignment_char}{ra_value})"]
    pattern_list_dec += [f"({dec_text}{pattern_J2000}{assignment_char}{dec_value})"]


    ### TABLES
    pattern_list_table += [f"([0-9]{{1,2}}(\\:)[0-9]{{1,2}}(\\:)[0-9]{{1,2}}(\\.|)[0-9]{{0,}})((( |)(,|\\|)( |))|( ))(\\+|-|)([0-9]{{1,2}}(\\:)[0-9]{{1,2}}(\\:)[0-9]{{1,2}}(\\.|)[0-9]{{0,}})"]
    pattern_list_table += [f"([0-9]{{1,2}}( )[0-9]{{1,2}}( )[0-9]{{1,2}}(\\.|)[0-9]{{0,}})((( |)(,|\\|)( |))|( ))(\\+|-|)([0-9]{{1,2}}( )[0-9]{{1,2}}( )[0-9]{{1,2}}(\\.|)[0-9]{{0,}})"]
    pattern_list_table += [f"([0-9]{{1,2}}(h)[0-9]{{1,2}}{units_minutes}[0-9]{{1,2}}(\\.|)[0-9]{{0,}}{units_seconds})((( |)(,|\\|)( |))|( ))(\\+|-|)([0-9]{{1,2}}(d|(deg))[0-9]{{1,2}}{units_minutes}[0-9]{{1,2}}(\\.|)[0-9]{{0,}}{units_seconds})"]


    ### PAIRS
    pattern_list_ra_dec += [f"\(j2000 {ra_dec}\){assignment_char}\({ra_value_deg}( |)(,|)( |){dec_value_deg}\)( |){degree_}"]
    pattern_list_ra_dec += [f"\({ra_dec}( |)(j|)2000(.0|)\){assignment_char}(\(|){ra_value}( |)(,|)( |){dec_value}(\)|)"]
    pattern_list_ra_dec += [f"\({ra_dec} {ra_value}( |)(,)( |){dec_value}\)"]
    pattern_list_ra_dec += [f"({ra_text}( |)\((j|)2000(.0|)\) {ra_value}), ({dec_text}( |)\((j|)2000(.0|)\) {dec_value})"]

    pattern_list_ra_dec += [f"\\b({ra_text} {ra_value}(, | |; |,|;){dec_text} {dec_value})\\b"]
    pattern_list_ra_dec  += [f"\\b({ra_text}{assignment_char}{ra_value})( )({dec_text}{assignment_char}{dec_value})\\b"]
    
    pattern_list_ra_dec += [f"\\b{ra_dec}{assignment_char}{ra_value}( |)(,|)( |){dec_value}\\b"]
    pattern_list_ra_dec += [f"\({ra_dec}\){assignment_char}(\(|){ra_value}( |)(,|)( |){dec_value}(\)|)"]

    pattern_list_ra_dec += [f"({ra_text}(\\/|,|, ){dec_text}( |)(\((j|)2000(.0|)\)|){assignment_char}{ra_value}( |,|, ){dec_value})"]
    
    pattern_list_ra_dec += [f"({ra_text}( and ){dec_text} {ra_value}( and ){dec_value})"]

    return pattern_list_ra_dec, pattern_list_ra, pattern_list_dec, pattern_list_table


def ra_dec_detector(text_id, text_id_text):
    pattern_list_ra_dec, pattern_list_ra, pattern_list_dec, pattern_list_table = create_pattern_list()

    counter = 0

    text_id_text = " ".join(text_id_text.split()).replace("Â°", "o").replace("Âº", "o").replace("−","-").replace('°', "o")
    list_phrases = split_text_in_phrases(text_id, text_id_text.lower())

    dict_data = {"TEXT_ID": [], "Positions": [], "Start": [], "End": [], "Phrase": []}

    for phrase_ in list_phrases:
        for pattern_ in pattern_list_ra_dec + pattern_list_ra + pattern_list_dec + pattern_list_table:
            for m in re.finditer(pattern_, phrase_.lower()):
                pos_ = m.group(0)
                start, end = m.span()

                dict_data["TEXT_ID"].append(text_id)
                dict_data["Start"].append(start)
                dict_data["End"].append(end)
                dict_data["Positions"].append(pos_)
                dict_data["Phrase"].append(phrase_)


    df_data = pd.DataFrame(dict_data)
    return df_data


def merge_ra_dec(text_id, df_init):
    df_init.drop_duplicates(inplace=True)
    dict_data = {"TEXT_ID": [], "Positions": [], "Start": [], "End": [], "Phrase": []}
    phrases_, c = np.unique(df_init["Phrase"].values, return_counts=True)

    for p_n, phrase_ in enumerate(phrases_):
        df_tmp0 = df_init[df_init["Phrase"] == phrase_]
        if len(df_tmp0) > 1:
            df_tmp = df_tmp0.sort_values("Start")
            start_ = df_tmp["Start"].values
            end_   = df_tmp["End"].values
            for i in range(1, len(start_)):
                if start_[i] <= end_[i-1]:
                    start_[i] = start_[i-1]
                    max_ = max(end_[i-1], end_[i])
                    end_[i-1], end_[i] = max_, max_
                    end_[i-1] = -1
                    start_[i-1] = -1

            for s_i, e_i in zip(start_, end_):
                if s_i != -1:
                    dict_data["TEXT_ID"]   += [text_id]
                    dict_data["Start"]     += [s_i]
                    dict_data["End"]       += [e_i]
                    dict_data["Positions"] += [phrase_[s_i: e_i]]
                    dict_data["Phrase"]    += [phrase_]

    df_data = pd.DataFrame(dict_data)
    df_data.drop_duplicates(inplace=True)
    return df_data


def clean_ra(ra, ra_text, pattern_J2000):    
    ra_new = " ".join(str(ra).split()).replace("±", "+/-").replace("—", "-").replace("−", "-").replace("−","-")
    
    ra_new = re.sub(f"{ra_text}{pattern_J2000}", "" , ra_new)
    ra_new = re.sub(f"{ra_text}", "" , ra_new)
    
    ra_new = re.sub("[^0-9+-\.deg]",":", ra_new)
     
    while len(ra_new) > 1 and (ra_new[-1] in [":", "."]):
        ra_new = ra_new[:-1]
    
    while len(ra_new) > 1 and (ra_new[0] in [":", "."]):
        ra_new = ra_new[1:]
    
    result = re.match(f"(\\+|)[0-9]{{1,2}}[:]{{1,2}}[0-9]{{1,2}}[:]{{1,2}}[0-9]{{1,2}}(:\\.|\\.|)([0-9]){{0,}}", ra_new)
    if result:
        if result.group(0) == ra_new:
            ra_new = ra_new.replace("::", ":")
            ra_new = ra_new.replace(":.", ".")
    
    ### Remove some incorect pos
    result = re.match(f"(\\+|)[0-9]{{4,}}(\\.|)([0-9]){{0,}}", ra_new)
    if result:
        if result.group(0) == ra_new:
            ra_new = ":"
        
    ra_new = ra_new.replace(":deg", " deg")
    
    return ra_new


def clean_dec(dec, dec_text, pattern_J2000):    
    dec_new = " ".join(str(dec).split()).replace("±", "+/-").replace("—", "-").replace("−", "-").replace("−","-").replace("--", "-")

    dec_new = re.sub(f"{dec_text}{pattern_J2000}", "" , dec_new)
    dec_new = re.sub(f"{dec_text}", "" , dec_new)

    dec_new = re.sub("[^0-9+-\.deg]",":", dec_new)
    
    while len(dec_new) != 1 and (dec_new[-1] in [":", "."]):
        dec_new = dec_new[:-1]
    
    while len(dec_new) != 1 and (dec_new[0] in [":", "."]):
        dec_new = dec_new[1:]
    
    result = re.match(f"(\\+|\\-|)[0-9]{{1,2}}(deg|d|:)[:]{{0,1}}[0-9]{{1,2}}[:]{{1,2}}[0-9]{{1,2}}(:\\.|\\.|)([0-9]){{0,}}", dec_new)
    if result:
        if result.group(0) == dec_new:
            dec_new = dec_new.replace("deg:", ":")
            dec_new = dec_new.replace("deg", ":")
            dec_new = dec_new.replace("d:", ":")
            dec_new = dec_new.replace("d", ":")
            dec_new = dec_new.replace("::", ":")
            dec_new = dec_new.replace(":.", ".")
    
    dec_new = dec_new.replace(":deg", " deg")
    
    ### Remove some incorect pos
    result = re.match(f"(\\+|\\-|)[0-9]{{4,}}(\\.|)([0-9]){{0,}}", dec_new)
    if result:
        if result.group(0) == dec_new:
            dec_new = ":"
    
    return dec_new


def clean_ra_dec(ra_dec, ra_text, dec_text, pattern_J2000):
    ra_dec_n = re.sub(f"{pattern_J2000}",   "",  ra_dec)
    ra_dec_n = re.sub(f"{ra_text}",         "",  ra_dec_n)
    ra_dec_n = re.sub(f"{dec_text}",        "",  ra_dec_n)
    if ra_dec_n[-1] == "o":
        ra_dec_n = ra_dec_n[:-1]
    ra_dec_n = re.sub("(o)", "d", ra_dec_n)
    ra_dec_n = re.sub("('')", "", ra_dec_n)
    ra_dec_n = re.sub("(')", "m", ra_dec_n)
    ra_dec_n = re.sub("[^0-9+-\.hmd\s:]", "", ra_dec_n)

    ra_dec_n = re.sub("[,]", "", ra_dec_n)
    while len(ra_dec_n) != 1 and (ra_dec_n[-1] in [":", ".", " "]):
        ra_dec_n = ra_dec_n[:-1]

    while len(ra_dec_n) != 1 and (ra_dec_n[0] in [":", ".", " "]):
        ra_dec_n = ra_dec_n[1:]
        
    return ra_dec_n


def astropy_test(df_init):
    ra_text  = "(r(\\.|)a(\\.|\\:|))"
    dec_text = "(dec(l|)(\\.|\\:|))"

    ra_dec_pattern = f"({ra_text},( |){dec_text})"
    pattern_J2000 = "( |)(-|)((\((((j|)2000(.0|))|(deg))\))|(2000))"

    rest_ra_dec  = []
    counter_rest = 0
    counter_ra_dec_try = 0
    good_ra_dec = []

    for text_ in list(set(df_init.Phrase)):
        df_tmp1 = df_init[df_init.Phrase == text_]
        ra_values  = []
        dec_values = []

        ra_start  = []
        dec_start = []
        ra_end    = []
        dec_end   = []

        df_tmp1.sort_values("Start")

        for ra_dec, s_, e_ in zip(df_tmp1.Positions, df_tmp1.Start, df_tmp1.End):    
            try:
                ra_dec_n = ra_dec.replace("|", " ")
                ra_dec_n = ra_dec_n.replace(",", " ")
                cc = SkyCoord(ra_dec_n, unit=(u.hourangle, u.deg))
                good_ra_dec.append(cc)
            except:
                ra_s = re.findall(ra_text, ra_dec)
                dec_s = re.findall(dec_text, ra_dec)
                if len(ra_s) > 0 and len(dec_s) > 0:
                    counter_ra_dec_try += 1
                    ra_dec_n = clean_ra_dec(ra_dec, ra_text, dec_text, pattern_J2000)
                    try:
                        cc = SkyCoord(ra_dec_n, unit=(u.hourangle, u.deg))
                        good_ra_dec.append(cc)
                    except:
                        rest_ra_dec.append(ra_dec)

                elif len(ra_s) > 0:
                    ra_values.append(ra_dec)
                    ra_start.append(s_)
                    ra_end.append(e_)

                elif len(dec_s) > 0:
                    dec_values.append(ra_dec)
                    dec_start.append(s_)
                    dec_end.append(e_)

                else:
                    counter_rest += 1

        if len(ra_values) < len(dec_values):

            for ra_, s_ra, e_ra in zip(ra_values, ra_start, ra_end):
                min_diff = 1000
                dec_pair = ""
                for dec_, s_dec, e_dec in zip(dec_values, dec_start, dec_end):
                    diff_ = s_dec - e_ra
                    if diff_ < min_diff:
                        min_diff = diff_
                        dec_pair = dec_

                c_ra  = clean_ra(ra_, ra_text, pattern_J2000)
                c_dec = clean_dec(dec_pair, dec_text, pattern_J2000)
                try:
                    cc = SkyCoord(ra=c_ra, dec=c_dec, unit=(u.hourangle, u.deg))
                    good_ra_dec.append(cc)
                except:
                    try:
                        cc = SkyCoord(ra=c_ra, dec=c_dec, unit=(u.deg, u.deg))
                        good_ra_dec.append(cc)
                    except:
                        rest_ra_dec.append(f"{ra_}|{dec_pair}")
        else:

            for dec_, s_dec, e_dec in zip(dec_values, dec_start, dec_end):
                min_diff = 1000
                ra_pair = ""
                for ra_, s_ra, e_ra in zip(ra_values, ra_start, ra_end):
                    diff_ = s_dec - e_ra
                    if diff_ < min_diff:
                        min_diff = diff_
                        ra_pair = ra_

                c_ra  = clean_ra(ra_pair, ra_text, pattern_J2000)
                c_dec = clean_dec(dec_, dec_text, pattern_J2000)
                try:
                    cc = SkyCoord(ra=c_ra, dec=c_dec, unit=(u.hourangle, u.deg))
                    good_ra_dec.append(cc)
                except:
                    try:
                        cc = SkyCoord(ra=c_ra, dec=c_dec, unit=(u.deg, u.deg))
                        good_ra_dec.append(cc)
                    except:
                        rest_ra_dec.append(f"{ra_pair}|{dec_}")
                        
    return good_ra_dec


def rule_based_ra_dec_detector(text_id, text_id_text):
    df_init  = ra_dec_detector(text_id, text_id_text)
    df_final = merge_ra_dec(text_id, df_init)
    good_ra_dec = astropy_test(df_final)
    print(good_ra_dec)
    dict_out = {"TEXT_ID":[], "RA":[], "DEC":[]}
    for ra_dec in good_ra_dec:
        dict_out["TEXT_ID"].append(text_id)
        dict_out["RA"].append(ra_dec.ra.deg)
        dict_out["DEC"].append(ra_dec.dec.deg)

    return pd.DataFrame(dict_out)
