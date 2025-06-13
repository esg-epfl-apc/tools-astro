import pandas as pd
import numpy as np
import re
import time
import io
import requests
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy import units as u

from query_tns_aux import parse_data, query_tns_main_name, query_tns_survey_name


def query_simbad(name):
    Simbad.add_votable_fields("otypes")

    table_ids = Simbad.query_objectids(name)
    table_obj = Simbad.query_object(name)

    dict_data = {}

    list_ids = None
    main_id = None
    otype_ = None
    ra = None
    dec = None
    if table_ids:
        df_ids = table_ids.to_pandas()
        df_ids.columns = df_ids.columns.str.upper()
        list_ids = list(df_ids["ID"].str.decode("utf-8").values)
    if table_obj:
        df_obj = table_obj.to_pandas()
        df_obj.columns = df_obj.columns.str.upper()
        main_id = df_obj["MAIN_ID"].values[0]
        otype_ = '|'.join(list(set(df_obj["OTYPES.OTYPE"].values)))
        str_ra_dec = str(df_obj["RA"].values[0]) + " " + str(df_obj["DEC"].values[0])
        ra_dec = SkyCoord(str_ra_dec, unit=(u.hourangle, u.deg))
        ra = ra_dec.ra.value
        dec = ra_dec.dec.value

    dict_data[name] = {"IDs": list_ids, "MAIN_ID": main_id, "OTYPES": otype_, "RA": ra, "DEC": dec, "DISCOVERY_TIME": None}

    return dict_data


def query_tns(name):
    dict_data = {}

    if ((name[0:3] == "at2") or (name[0:2] == "sn")):
        my_text = query_tns_main_name(name[2:])
    else:
        my_text = query_tns_survey_name(name)

    time.sleep(5)

    main_id, list_ids, list_otypes, ra, dec, discovery_time = parse_data(my_text)
    if ra and dec:
        str_ra_dec = ra + ":" + dec
        ra_dec = SkyCoord(str_ra_dec.replace(":", " "), unit=(u.hourangle, u.deg))
        ra = ra_dec.ra.value
        dec = ra_dec.dec.value

    if list_otypes is None:
        otype_ = None
    elif len(list_otypes) == 1:
        otype_ = list_otypes[0]
    else:
        otype_ = '|'.join(list_otypes)

    dict_data[name] = {"IDs": list_ids, "MAIN_ID": main_id, "OTYPES": otype_, "RA": ra, "DEC": dec, "DISCOVERY_TIME": discovery_time}

    return dict_data


def query_fink(name):
    dict_data = {}
    list_ids = None
    main_id = None
    otype_ = None
    ra = None
    dec = None
    discovery_time = None

    query_name = name.replace("ztf", "ZTF")
    r = requests.post('https://api.fink-portal.org/api/v1/objects', json={'objectId': query_name, 'output-format': 'json', 'columns': 'i:ra,i:dec,i:objectId'})
    # Format output in a DataFrame
    df = pd.read_json(io.BytesIO(r.content))
    if not df.empty:
        object_id = df["i:objectId"].values[0]
        main_id = object_id
        ra_tmp = df["i:ra"][df["i:objectId"] == object_id].values
        dec_tmp = df["i:dec"][df["i:objectId"] == object_id].values
        ra = np.mean(np.float32(ra_tmp))
        dec = np.mean(np.float32(dec_tmp))
        list_ids = '|'.join(list(set(df["i:objectId"].values)))

    dict_data[name] = {"IDs": list_ids, "MAIN_ID": main_id, "OTYPES": otype_, "RA": ra, "DEC": dec, "DISCOVERY_TIME": discovery_time}

    return dict_data


def create_pattern_list():
    pattern_list_low = []
    pattern_list_low += ["\\b(ngc) *?([0-9]{1,4})\\b", "\\b(m) *?([0-9]{1,3})\\b"]
    pattern_list_low += ["\\b(ugc) *?([0-9]{1,5})\\b"]

    pattern_list_low += ["\\b(icecube|grb|frb|pks|mrk|hawc|maxi|gw)([ -]?)([0-9\\.\\-\\+]{2,}[a-z]?)\\b"]
    pattern_list_low += ["\\b(ic)([ -]?)([0-9]{1,4})\\b"]
    pattern_list_low += ["\\b(ztf) *?([0-9]{2}[a-z]{7})\\b"]

    pattern_list_low += ["\\b(at|sn) *?([1-2]{1}[0-9]{3}[a-z]{1,4})\\b"]
    pattern_list_low += ["\\b(asas)([ -]?)(sn)([ -]?)([0-9]{2}[a-z]{2,3})\\b"]

    pattern_list_low += ["\\b(ps|gaia) *?([1-2]{1}[0-9]{1}[a-z]{1,4})\\b"]

    pattern_list_low += ["\\b(m31n) *?([0-9]{4})-([0-9]{2}[a-z]{1})\\b"]
    pattern_list_low += ["\\b(ptf|atlas) *?([0-9]{2}[a-z]{1,4})\\b"]

    pattern_list_low += ["\\b(4c) *?((\\+|-)[0-9]{2}\\.[0-9]{2}\\.[0-9]{1})\\b"]
    pattern_list_low += ["\\b(4c) *?((\\+|-)[0-9]{2}\\.[0-9]{2})\\b"]

    # new from lm-astronomy
    pattern_list_low += ["\\b(lsq12) *?([a-z]{3})\\b"]
    pattern_list_low += ["\\b(des14) *?([a-z]{1}[0-9]{1})\\b"]

    pattern_list = []

    HH = "(2[0-4]|[0-1][0-9])"
    MM = "(90|[0-8][0-9])"
    MMm = "(90|[0-8][0-9])([0-9])"
    SS = "(90|[0-8][0-9])"
    DD = "(90|[0-8][0-9])"
    pDDd = "(\\+|-)(90|[0-8][0-9])([0-9])"

    LLL = "([0-9][0-9][0-9])"
    BB = "([0-9][0-9])"
    VVV = "([0-9][0-9][0-9])"

    NAAA = "(([0-9][A-Z]{1,4})|([A-Z]{2,4}))"

    pattern_list += [f"\\b{NAAA} *?([0-9]{{6}})(?!-)\\b"]

    pattern_list += [f"\\b{NAAA} *?(J[0-9]{{6}}(\\+|-)[0-9]{{5}})(?!-)\\b"]
    pattern_list += [f"\\b{NAAA} *?(J[0-9]{{6}}(\\.?)[0-9]{{1}}(\\+|-)[0-9]{{6}})(?!-)\\b"]

    pattern_list += [f"\\b{NAAA} *?(J[0-9]{{4}}(\\.?)[0-9]{{1}}(\\+|-)[0-9]{{4}})[a-z]{{0,1}}(?!-)\\b"]
    pattern_list += [f"\\b{NAAA} *?{HH}{MM}(\\+|-){DD}(?!-)\\b"]
    pattern_list += [f"\\b{NAAA} *?(J?){HH}{MM}{pDDd}(?!-)\\b"]
    pattern_list += [f"\\b{NAAA} *?B{HH}{MM}(\\+|-){DD}(?!-)\\b"]
    pattern_list += [f"\\b{NAAA} *?B{HH}{MM}{pDDd}(?!-)\\b"]
    pattern_list += [f"\\b{NAAA} *?J{HH}{MM}(\\+|-){DD}{MM}(?!-)\\b"]
    pattern_list += [f"\\b{NAAA} *?J{HH}{MM}([0-9])(\\+|-){DD}{MM}(?!-)\\b"]

    pattern_list += [f"\\b{NAAA} *?{HH}{MMm}(\\+|-){DD}{MM}(?!-)\\b"]
    pattern_list += [f"\\b{NAAA} *?{HH}{MM}(\\.?)([0-9])(\\+|-){DD}{MM}(?!-)\\b"]
    pattern_list += [f"\\b{NAAA} *?{HH}{MM}{SS}(\\+|-){DD}{MM}{SS}(?!-)\\b"]
    pattern_list += [f"\\b{NAAA} *?{HH}{MM}{SS}(\\.?)([0-9])(\\+|-){DD}{MM}{SS}(?!-)\\b"]
    pattern_list += [f"\\b{NAAA} *?B{HH}{MM}{SS}(\\.?)([0-9])(\\+|-){DD}{MM}{SS}(?!-)\\b"]
    pattern_list += [f"\\b{NAAA} *?J{HH}{MM}{SS}(\\.?)([0-9][0-9])(\\+|-){DD}{MM}{SS}(\\.?)([0-9])(?!-)\\b"]

    pattern_list += [f"\\b{NAAA} *?{LLL}(\\.?)([0-9])(\\+|-){BB}(\\.?)([0-9])(\\+|-){VVV}(?!-)\\b"]
    pattern_list += [f"\\b{NAAA} *?G{LLL}(\\.?)([0-9])(\\+|-){BB}(\\.?)([0-9])(\\+|-){VVV}(?!-)\\b"]

    FFF = "([0-9][0-9][0-9])"
    FFFF = "([0-9][0-9][0-9][0-9])"
    NNNN = FFFF
    NNNNN = "([0-9][0-9][0-9][0-9][0-9])"

    pattern_list += [f"\\b{NAAA} *?{FFF}-{NNNN}(?!-)\\b"]
    pattern_list += [f"\\b{NAAA} *?{FFFF}-{NNNNN}(?!-)\\b"]

    YY = "[0-9][0-9]"
    MM = "[0-1][0-9]"
    DD = "[0-3][0-9]"
    pattern_list += [f"\\b{NAAA} *?{YY}{MM}{DD}A{{0,1}}(?!-)\\b"]

    return pattern_list, pattern_list_low


def rule_based_source_detector(text_id, text_id_text):
    pattern_list, pattern_list_low = create_pattern_list()

    regex_sources = []

    for pattern in pattern_list_low:
        for m in re.finditer(pattern, text_id_text.lower()):
            source_ = m.group(0).replace(" ", "")

            if "asas" in source_:
                source_ = source_.replace("-", "")
                source_ = source_.replace("sn", "sn-")

            regex_sources.append(source_)

    for pattern in pattern_list:
        for m in re.finditer(pattern, text_id_text):
            source_ = m.group(0)
            if source_.replace(" ", "").lower() not in regex_sources:
                regex_sources.append(source_)

    return list(set(regex_sources))


def query_info_sources(text_id, sources):
    if len(sources) != 0:
        source_list = []
        otype_list = []
        mainid_list = []
        ra_list = []
        dec_list = []
        pattern_string = re.compile("[a-z]")
        dict_unknown = {}
        dict_unknown = {"Raw Source Name": []}

        for source_name_0 in sources:
            # TODO:  remove all non a-z0-9 characters from the beginning and ending of the source name
            # Now: remove "," from the beginning and ending of the source name
            if source_name_0[0] == ",":
                source_name = source_name_0[1:]
            elif source_name_0[-1] == ",":
                source_name = source_name_0[:-1]
            else:
                source_name = source_name_0

            if pattern_string.findall(source_name.lower()):
                dict_otype = query_simbad(source_name)
                if dict_otype[source_name]["MAIN_ID"] is None:
                    dict_otype = query_tns(source_name)
                if dict_otype[source_name]["MAIN_ID"] is None:
                    dict_otype = query_fink(source_name)

                if dict_otype[source_name]["MAIN_ID"] is None:
                    dict_unknown["Raw Source Name"].append(source_name)
                else:
                    mainid_list.append(dict_otype[source_name]["MAIN_ID"])
                    otype_list.append(dict_otype[source_name]["OTYPES"])
                    ra_list.append(dict_otype[source_name]["RA"])
                    dec_list.append(dict_otype[source_name]["DEC"])
                    source_list.append(source_name)
            else:
                dict_unknown["Raw Source Name"].append(source_name_0)

        dict_data = {"TEXT_ID": [text_id] * len(source_list), "Raw Source Name": source_list, "Main ID Name": mainid_list, "OTYPE": otype_list, "RA": ra_list, "Dec": dec_list}
        df_save = pd.DataFrame(dict_data)
        df_save.replace({None: "NotKnown"}, inplace=True)
        return df_save.drop_duplicates(subset=['Main ID Name']), pd.DataFrame(dict_unknown)

    return pd.DataFrame(), pd.DataFrame()
