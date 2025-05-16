import pandas as pd
import re


def rule_based_class_detector(simbad_node_file, text_id_text):
    df = pd.read_csv(simbad_node_file)
    pattern_list = list(df["Description"].values)

    classes = []

    for pattern in pattern_list:
        for m in re.finditer(f"\\b{pattern.lower()}\\b", text_id_text):
            source_ = m.group(0)
            classes.append(source_)

    return classes


def source_class(df_in, simbad_node_file):
    out_class_list = []
    if len(df_in) > 0:
        df_dict = pd.read_csv(simbad_node_file)

        class_list = []

        otypes_ = df_in["OTYPE"].values
        for otypes in otypes_:
            if otypes is not None:
                for otype in set(otypes.split("|")):
                    class_list.append(otype)

        for otype in set(class_list):
            if "?" in otype:
                out_class_list.append(otype)
            classes = df_dict["Description"][df_dict["Id"] == otype].values
            if len(classes) != 0:
                out_class_list.append(classes[0])

    return out_class_list


def detect_source_classes(text_id, text_id_text, df_sources, simbad_node_file):
    classes_1 = rule_based_class_detector(simbad_node_file, text_id_text.lower())
    classes_2 = source_class(df_sources, simbad_node_file)
    classes = classes_1 + classes_2

    if len(classes) != 0:
        out_classes = list(set(classes))

        dict_data = {"TEXT_ID": [text_id] * len(out_classes), "Source Classes": out_classes}
        df_data = pd.DataFrame(dict_data)
        return df_data

    return pd.DataFrame()
