import numpy as np
import pandas as pd

from pipeline_telescope import rule_based_telescope_detector
from pipeline_sources import rule_based_source_detector, query_info_sources
from pipeline_source_classes import detect_source_classes
from pipeline_astrobert import get_astroBERT_cleaned_result
from pipeline_vectorize_text import vectorize_text
from predict_vectorised_text import predict_vector
from pipeline_create_url_vector import create_url_vector

if __name__ == "__main__":

    ### Settings
    data_path = "/home/jovyan/work/analyse-short-astro-text/data/"
    model_file = f"{data_path}/probabilistic_models/first_single_second_single_with_tel_repetition.dat.npy"

    ### Input
    # text_id_text=" We report photometry of Type II SN 2024iss (TNSTR-2024-1494, TNSCR-2024-1517) gathered with the OAUNI 51cm telescope (arXiv:1512.03104) at Huancayo Observatory, Peru. CCD imaging in R filter was performed on three consecutive nights under non-photometric conditions with airmasses lower than 1.5. The integration time was 70x20s=1400s in each observation. Our measurements yielded: UCAC4 field stars were used for the zero point calibration. Our observations are ~3.7 weeks after the GOTO detection (TNSTR-2024-1494). The OAUNI project is supported by UNI, TWAS, IGP and ProCiencia-Concytec (Convenio 133-2020 Fondecyt). "
    # text_id = "ATel #16662"
    
    # text_id_text = " On 29 August 2024, we performed follow-up observations of candidate optical counterparts of newly-discovered Einstein Probe X-ray transients, EP240809a (ATel #16765) and EP240709a (ATel #16704), with the Southern African Large Telescope (SALT). The observations were executed using the Robert Stobie Spectrograph (RSS) in low-resolution mode (R ~ 1000). EP240709a was observed with a total exposure time of 2400 secs (starting at UT 21:50:36) with a grating angle of 14.75o, covering a wavelength range of ~4000 - 7000 Angstroms. The extracted 1D spectrum is devoid of any prominent lines. EP240809a was observed with a total exposure time of 3200 secs (starting at UT 20:41:37) and a grating angle of 15.125o, covering a wavelength region of ~4100-7200 Angstroms. The 1D spectrum reveals no significant lines and a steep rising red continuum."
    # text_id = "ATel #16801"
    
    text_id_text = '''During the ongoing All Sky Automated Survey for SuperNovae (ASAS-SN or "Assassin"), using data from the quadruple 14-cm "Brutus" telescope in Haleakala, Hawaii, we discovered a new transient source, most likely a supernova, in the galaxy 2MASX J15591858+1336487. ASASSN-17bd (AT 2017nk) was discovered in images obtained on UT 2017-01-23.61 at V~17.3 mag. We do not detect (V>17.5) the object in images taken on UT 2017-01-16.66 and before. An image obtained on 2017-01-23 by J. Brimacombe confirms the discovery of the transient. This figure shows the archival SDSS g-band image of the host (left) and the J. Brimacombe confirmation image (right). The red circle has a radius of 5" and is centered on the position of the transient in the J. Brimacombe image. The position of ASASSN-17bd is approximately 1.8" North and 2.4" West from the center of the galaxy 2MASX J15591858+1336487 (z=0.034554, d=147 Mpc, via NED), giving an absolute V-band magnitude of approximately -18.6 (m-M=35.77, A_V=0.129). Properties of the new source and photometry are summarized in the tables below: Object       RA (J2000)     DEC (J2000)      Disc. UT Date   Disc. V mag  Approx. Abs. Mag   Offset from Host (")  ASASSN-17bd  15:59:18.432   +13:36:50.89    2017-01-23.61   17.3           -18.6               3.0  Obs. UT Date         V mag  2017-01-16.66        >17.5 2017-01-23.61         17.3  Follow-up observations are encouraged. While we are participating in the TNS system to minimize potential confusion, ASAS-SN will continue using ASASSN-17xx transient names as our primary nomenclature (including supernovae, but also other classes of transients), and we encourage others to do the same. We prefer merging the names as ASASSN-17xx (AT 2017xyz) to preserve, rather than anonymize, the origin of the transient. We thank Las Cumbres Observatory and its staff for their continued support of ASAS-SN. ASAS-SN is funded in part by the Gordon and Betty Moore Foundation through grant GBMF5490 to the Ohio State University, NSF grant AST-1515927, the Center for Cosmology and AstroParticle Physics (CCAPP) at OSU, the Chinese Academy of Sciences South America Center for Astronomy (CASSACA), and the Mt. Cuba Astronomical Foundation. For more information about the ASAS-SN project, see the ASAS-SN Homepage and the list of all ASAS-SN transients.'''
    text_id = "ATel #10000"
    
    ### Required files
    telescope_ontology   = f"{data_path}/telescope_observatory_survey.ttl"
    file_dict_uri_int    = f"{data_path}/dictionary_telescope_uri_int_id.json"
    file_dict_sens_inst  = f"{data_path}/dictionary_telescope_type_2_instrument.json"

    simbad_node_file  = f"{data_path}/simbad_otypes_nodes.csv"

    ### Run pipeline
    df_tel = rule_based_telescope_detector(text_id, text_id_text, telescope_ontology, file_dict_uri_int)
    df_astrobert = get_astroBERT_cleaned_result(text_id, text_id_text)

    astrobert_sources = list(df_astrobert[df_astrobert["entity_group"] == "CelestialObject"].dropna().word.values)
    regex_sources = rule_based_source_detector(text_id, text_id_text)

    df_sor, df_unk_sor = query_info_sources(text_id, list(set(regex_sources + astrobert_sources)))

    df_cla = detect_source_classes(text_id, text_id_text.lower(), df_sor, simbad_node_file)

    df_vec = vectorize_text(text_id, data_path, df_tel, df_sor, df_astrobert)

    df_vec_init_pred = predict_vector(text_id, data_path, df_vec)
    df_vec_url, df_url_scores = create_url_vector(text_id, data_path, file_dict_sens_inst, df_vec_init_pred, df_sor)

    print(df_tel)
    print(df_sor)
    print(df_unk_sor)
    print(df_cla)
    print(df_astrobert)

    print(df_vec_init_pred)
    print(df_vec_url)
    print(df_url_scores)