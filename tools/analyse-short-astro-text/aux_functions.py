import numpy as np


def list_tel(key, G):
    cols = []

    for e in G.query(f'SELECT ?a ?b ?c WHERE {{<{key}> <https://odahub.io/ontology#sensitivity> ?c}} LIMIT 500'):
        cols.append(str(e[2]))

    if len(cols) == 0:
        for e in G.query(f'SELECT ?a ?b ?c WHERE {{?a <http://purl.org/dc/terms/isRequiredBy> <{key}>}} LIMIT 500'):
            for f in G.query(f'SELECT ?a ?b ?c WHERE {{<{e[0]}> <https://odahub.io/ontology#sensitivity> ?c}} LIMIT 500'):
                cols.append(str(f[2]))

    # if len(cols) == 0:
    #     aux_ = []
    #     for key in dict_.keys():
    #         if dict_[key] == -i:
    #             for e in G.query(f'SELECT ?a ?b ?c WHERE {{?a <http://purl.org/dc/terms/isPartOf> <{key}>}} LIMIT 500'):
    #                 # for f in G.query(f'SELECT ?a ?b ?c WHERE {{<{e[0]}> <https://odahub.io/ontology#sensitivity> ?c}} LIMIT 500'):
    #                 aux_.append(e[0])
    #             if len(aux_) == 1:
    #                 print(aux_)

    return list(set(cols))


def provide_telescope_sensitivities(dict_, G):
    colors = []
    N_TEL = len(dict_.keys())
    for key in dict_.keys():
        for e in G.query(f'SELECT ?a ?b ?c WHERE {{<{key}> <https://odahub.io/ontology#sensitivity> ?c}} LIMIT 500'):
            colors.append(str(e[2]))

    return list(set(colors)), N_TEL


def bits(ask="try"):
    if ask == "None":
        return 0    # (0 0 0 0 0 0 0 0 0)
    if ask == "gamma-ray":
        return 1    # (0 0 0 0 0 0 0 0 1)
    if ask == "x-ray":
        return 2    # (0 0 0 0 0 0 0 1 0)
    if ask == "ultraviolet":
        return 4    # (0 0 0 0 0 0 1 0 0)
    if ask == "optical":
        return 8    # (0 0 0 0 0 1 0 0 0)
    if ask == "infrared":
        return 16   # (0 0 0 0 1 0 0 0 0)
    if ask == "radio":
        return 32   # (0 0 0 1 0 0 0 0 0)
    if ask == "cosmic-ray":
        return 64   # (0 0 1 0 0 0 0 0 0)
    if ask == "gravitational-wave":
        return 128  # (0 1 0 0 0 0 0 0 0)
    if ask == "neutrino":
        return 256  # (1 0 0 0 0 0 0 0 0)
    # print(f"You have asked for {ask}. Which does not exist. Please check bits() function.")
    # os._exit(1)


def bits2sens(ask=None):
    if ask == 0:
        return "None"
    if ask == 1:
        return "gamma-ray"
    if ask == 2:
        return "x-ray"
    if ask == 4:
        return "ultraviolet"
    if ask == 8:
        return "optical"
    if ask == 16:
        return "infrared"
    if ask == 32:
        return "radio"
    if ask == 64:
        return "cosmic-ray"
    if ask == 128:
        return "gravitational-wave"
    if ask == 256:
        return "neutrino"
    # print(f"You have asked for {ask}. Which does not exist. Please check bits() function.")
    # os._exit(1)


def compute_sensitivity(colors):
    final_ = 0
    for color in colors:
        final_ = np.bitwise_or(bits(ask=color), final_)

    return final_


def compute_sensitivity_int(colors_int):
    final_ = 0
    for color in colors_int:
        final_ = np.bitwise_or(color, final_)

    return final_


def find_sensitivity(bit_integer):
    list_o = []
    list_b = []
    list_b_12 = []

    colors = ['gravitational-wave', 'radio', 'ultraviolet', 'optical', 'gamma-ray', 'infrared', 'x-ray', 'cosmic-ray', 'neutrino']

    for sens in colors:
        bit_ = bits(ask=sens)
        if np.bitwise_and(bit_, bit_integer) == bit_:
            list_o.append(sens)
            list_b.append(bit_)

    for i_1 in range(len(colors)):
        s_1 = colors[i_1]
        b_1 = bits(ask=s_1)
        for i_2 in range(i_1, len(colors), 1):
            s_2 = colors[i_2]
            b_2 = bits(ask=s_2)

            if i_1 != i_2:
                b_12 = np.bitwise_or(b_1, b_2)

                if np.bitwise_and(b_12, bit_integer) == b_12:
                    list_b_12.append(b_12)
            else:
                if b_2 == bit_integer:
                    list_b_12.append(b_2)

    return list_o, list_b, list_b_12


def find_name_workflow(list_telescopes):
    indices = []
    if "https://odahub.io/ontology#international-gamma-ray-astrophysics-laboratory" in list_telescopes:
        indices.append(50)
    if "https://odahub.io/ontology#imager-on-board-the-integral-satellite" in list_telescopes:
        indices.append(51)
    if "https://odahub.io/ontology#joint-european-x-ray-monitor" in list_telescopes:
        indices.append(52)
    if "https://odahub.io/ontology#spectrometer-on-integral" in list_telescopes:
        indices.append(53)
    if "https://odahub.io/ontology#astronomy-with-a-neutrino-telescope-and-abyss-environmental-research-project" in list_telescopes:
        indices.append(54)
    if "https://odahub.io/ontology#laser-interferometer-gravitational-wave-observatory" in list_telescopes:
        indices.append(55)
    if "https://odahub.io/ontology#virgo" in list_telescopes:
        indices.append(55)
    if "https://odahub.io/ontology#iceCube-neutrino-observatory" in list_telescopes:
        indices.append(56)
    if "https://odahub.io/ontology#high-energy-stereoscopic-system" in list_telescopes:
        indices.append(57)
    if "https://odahub.io/ontology#cherenkov-telescope-array" in list_telescopes:
        indices.append(58)
    if "https://odahub.io/ontology#cherenkov-telescope-array-observatory" in list_telescopes:
        indices.append(58)
    return indices


def get_list_instruments_CNN_MMODA():
    return ["INTEGRAL", "ISGRI", "JEM-X", "SPI-ACS", "ANTARES", "LIGO/VIRGO", "IceCube", "HESS", "CTA/CTAO"]


def get_dict_instruments_URL_MMODA():
    return {
        "spi_acs": ["INTEGRAL", "SPI-ACS"],
        "cta": ["CTA/CTAO"],
        "hess": ["HESS"],
        "isgri": ["INTEGRAL", "ISGRI"],
        "jemx": ["INTEGRAL", "JEM-X"],
        "icecube": ["IceCube"],
        "antares": ["ANTARES"],
        "gw": ["LIGO/VIRGO"]
    }
