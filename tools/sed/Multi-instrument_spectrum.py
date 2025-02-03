#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from numpy import sqrt
from oda_api.api import DispatcherAPI, ProgressReporter
from oda_api.json import CustomJSONEncoder
from oda_api.token import discover_token

pr = ProgressReporter()

from oda_api.token import discover_token

token = discover_token(allow_invalid=True)

workdir = os.getcwd()

src_name = "Crab"  # http://odahub.io/ontology#AstrophysicalObject
RA = 83.6325  # http://odahub.io/ontology#PointOfInterestRA
DEC = 22.0175  # http://odahub.io/ontology#PointOfInterestDEC
T1 = "2000-03-06T13:26:48.0"  # http://odahub.io/ontology#StartTime
T2 = "2024-03-06T15:32:27.0"  # http://odahub.io/ontology#EndTime

do_mwa = True  # http://odahub.io/ontology#Boolean ; oda:label "MWA (radio)"
do_jemx = True  # http://odahub.io/ontology#Boolean ; oda:label "INTEGRAL/JEM-X (X-ray)"
do_isgri = True  # http://odahub.io/ontology#Boolean ; oda:label "INTEGRAL/ISGRI (hard X-ray)"
do_fermi = True  # http://odahub.io/ontology#Boolean ; oda:label "Fermi/LAT (gamma-ray)"
do_magic = (
    True  # http://odahub.io/ontology#Boolean ; oda:label "MAGIC (gamma-ray)"
)
do_icecube = (
    True  # http://odahub.io/ontology#Boolean ; oda:label "IceCube (neutrino)"
)
do_auger = (
    True  # http://odahub.io/ontology#Boolean ; oda:label "Auger (UHECR)"
)
do_hess = (
    True  # http://odahub.io/ontology#Boolean ; oda:label "HESS (gamma-ray)"
)
do_gaia = (
    True  # http://odahub.io/ontology#Boolean ; oda:label "GAIA (visible)"
)
do_legacysurvey = True  # http://odahub.io/ontology#Boolean ; oda:label "DESI Legacy Survey (visible/infrared)"

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in [
    "src_name",
    "RA",
    "DEC",
    "T1",
    "T2",
    "do_mwa",
    "do_jemx",
    "do_isgri",
    "do_fermi",
    "do_magic",
    "do_icecube",
    "do_auger",
    "do_hess",
    "do_gaia",
    "do_legacysurvey",
]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

FLAG_mwa = 0
FLAG_jemx = 0
FLAG_isgri = 0
FLAG_fermi = 0
FLAG_magic = 0
FLAG_hess = 0
FLAG_icecube = 0
FLAG_auger = 0
FLAG_gaia = 0
FLAG_legacysurvey = 0

disp_list = []
parameters = []

E0 = 20.0  # keV

def dn_de(E, A, Gam):
    return A * (E / E0) ** (-Gam)

def exp_counts_isgri(A, Gam):
    spectrum_dndE = dn_de(ENERG, A, Gam)
    tmp = np.zeros(NErec)
    for i in range(NErec):
        tmp[i] += sum(spectrum_dndE * dENERG * resp_isgri[:, i])
    return tmp

def exp_counts_jemx(A, Gam):
    spectrum_dndE = dn_de(ENERG, A, Gam)
    tmp = np.zeros(NErec)
    for i in range(NErec):
        # tmp[i]+=sum(spectrum_dndE*dENERG*resp_jemx[:,i])
        tmp[i] += sum(spectrum_dndE * dENERG * resp_jemx[:, i] * aeff)
    return tmp

pr.report_progress(stage="Sending requests for data products", progress=1)
if do_fermi:
    try:
        disp = DispatcherAPI(
            url="https://www.astro.unige.ch/mmoda//dispatch-data",
            instrument="mock",
            wait=False,
        )
        par_dict = {
            "Background_region_radius": 4.0,
            "RA": RA,
            "DEC": DEC,
            "Emax": 1000000.0,
            "Emin": 100.0,
            "NEbins": 8,
            "Radius": 15.0,
            "Source_region_radius": 2.0,
            "T1": "2017-03-06T13:26:48.000",
            "T2": "2017-07-06T15:32:27.000",
            "T_format": "isot",
            "instrument": "fermi_lat",
            "product": "Spectrum_aperture",
            "product_type": "Real",
            "src_name": src_name,
            "token": token,
        }
        data_collection_fermi = disp.get_product(**par_dict)
        disp_list.append(disp)
        parameters.append(par_dict)
    except Exception as e:
        print(e)

pr.report_progress(stage="Sending requests for data products", progress=2)
if do_jemx:
    try:
        disp = DispatcherAPI(
            url="https://www.astro.unige.ch/mmoda//dispatch-data",
            instrument="mock",
            wait=False,
        )
        par_dict = {
            "DEC": DEC,
            "E1_keV": 3.0,
            "E2_keV": 20.0,
            "RA": RA,
            "T1": T1,
            "T2": T2,
            "T_format": "isot",
            "detection_threshold": "7",
            "instrument": "jemx",
            "integral_data_rights": "public",
            "jemx_num": 1,
            "max_pointings": 50,
            "osa_version": "OSA11.2",
            "product": "jemx_spectrum",
            "product_type": "Real",
            "radius": 4.0,
            "src_name": src_name,
            "token": token,
        }
        data_collection_jemx1 = disp.get_product(**par_dict)
        disp_list.append(disp)
        parameters.append(par_dict)
    except Exception as e:
        print(e)

pr.report_progress(stage="Sending requests for data products", progress=3)
if do_jemx:
    try:
        disp = DispatcherAPI(
            url="https://www.astro.unige.ch/mmoda//dispatch-data",
            instrument="mock",
            wait=False,
        )
        par_dict = {
            "DEC": DEC,
            "E1_keV": 3.0,
            "E2_keV": 20.0,
            "RA": RA,
            "T1": T1,
            "T2": T2,
            "T_format": "isot",
            "detection_threshold": "7",
            "instrument": "jemx",
            "integral_data_rights": "public",
            "jemx_num": 2,
            "max_pointings": 50,
            "osa_version": "OSA11.2",
            "product": "jemx_spectrum",
            "product_type": "Real",
            "radius": 4.0,
            "src_name": src_name,
            "token": token,
        }
        data_collection_jemx2 = disp.get_product(**par_dict)
        disp_list.append(disp)
        parameters.append(par_dict)
    except Exception as e:
        print(e)

pr.report_progress(stage="Sending requests for data products", progress=4)
if do_isgri:
    try:
        disp = DispatcherAPI(
            url="https://www.astro.unige.ch/mmoda//dispatch-data",
            instrument="mock",
            wait=False,
        )
        par_dict = {
            "RA": RA,
            "DEC": DEC,
            "E1_keV": 28.0,
            "E2_keV": 40.0,
            "T1": T1,
            "T2": T2,
            "T_format": "isot",
            "detection_threshold": "7",
            "instrument": "isgri",
            "integral_data_rights": "public",
            "max_pointings": 50,
            "osa_version": "OSA11.2",
            "product": "isgri_spectrum",
            "product_type": "Real",
            "radius": 15.0,
            "src_name": src_name,
            "token": token,
        }
        data_collection_isgri = disp.get_product(**par_dict)
        disp_list.append(disp)
        parameters.append(par_dict)
    except Exception as e:
        print(e)

pr.report_progress(stage="Sending requests for data products", progress=5)
if do_legacysurvey:
    try:
        disp = DispatcherAPI(
            url="https://www.astro.unige.ch/mmoda//dispatch-data",
            instrument="mock",
            wait=False,
        )
        par_dict = {
            "DEC": DEC,
            "RA": RA,
            "Radius": 0.1,
            "data_release": 10,
            "instrument": "desi_legacy_survey",
            "product": "Spectrum",
            "product_type": "Real",
            "src_name": src_name,
            "token": token,
        }
        data_collection_legacysurvey = disp.get_product(**par_dict)
        disp_list.append(disp)
        parameters.append(par_dict)
    except Exception as e:
        print(e)

pr.report_progress(stage="Sending requests for data products", progress=6)
if do_gaia:
    try:
        disp = DispatcherAPI(
            url="https://www.astro.unige.ch/mmoda//dispatch-data",
            instrument="mock",
            wait=False,
        )
        par_dict = {
            "DEC": DEC,
            "RA": RA,
            "T1": T1,
            "T2": T2,
            "T_format": "isot",
            "data_release": 3,
            "instrument": "gaia",
            "product": "Spectrum_from_photometry",
            "product_type": "Real",
            "radius_photometry": 60.0,
            "src_name": src_name,
            "token": token,
        }
        data_collection_gaia = disp.get_product(**par_dict)
        disp_list.append(disp)
        parameters.append(par_dict)
    except Exception as e:
        print(e)

pr.report_progress(stage="Sending requests for data products", progress=7)
if do_hess:
    try:
        disp = DispatcherAPI(
            url="https://www.astro.unige.ch/mmoda//dispatch-data",
            instrument="mock",
            wait=False,
        )
        par_dict = {
            "RA": RA,
            "DEC": DEC,
            "Efit_max": 10.0,
            "Efit_min": 0.2,
            "Emax": 100.0,
            "Emin": 0.1,
            "NEbins": 15,
            "R_s": 0.2,
            "Radius": 1.0,
            "T1": T1,
            "T2": T2,
            "T_format": "isot",
            "instrument": "hess",
            "product": "Spectrum_public_dl3",
            "product_type": "Real",
            "src_name": "Crab",
            "token": token,
        }
        data_collection_hess = disp.get_product(**par_dict)
        disp_list.append(disp)
        parameters.append(par_dict)
    except Exception as e:
        print(e)

pr.report_progress(stage="Sending requests for data products", progress=8)
if do_magic:
    try:
        disp = DispatcherAPI(
            url="https://www.astro.unige.ch/mmoda//dispatch-data",
            instrument="mock",
            wait=False,
        )
        par_dict = {
            "RA": RA,
            "DEC": DEC,
            "T1": T1,
            "T2": T2,
            "Efit_max": 10.0,
            "Efit_min": 0.2,
            "Emax": 20.0,
            "Emin": 0.1,
            "NEbins": 30,
            "NSB": 0,
            "Offset": "0.4 deg",
            "R_s": 0.2,
            "Radius_search": 2.0,
            "T_format": "isot",
            "instrument": "magic",
            "product": "Spectrum_public_dl3",
            "product_type": "Real",
            "src_name": src_name,
            "token": token,
        }
        data_collection_magic = disp.get_product(**par_dict)
        disp_list.append(disp)
        parameters.append(par_dict)
    except Exception as e:
        print(e)

pr.report_progress(stage="Sending requests for data products", progress=9)
try:
    if do_mwa:
        disp = DispatcherAPI(
            url="https://www.astro.unige.ch/mmoda//dispatch-data",
            instrument="mock",
            wait=False,
        )
        par_dict = {
            "DEC": DEC,
            "RA": RA,
            "Radius": 0.1,
            "T1": T1,
            "T2": T2,
            "T_format": "isot",
            "instrument": "mwa",
            "product": "Spectrum",
            "product_type": "Real",
            "src_name": src_name,
            "token": token,
        }
        job = disp.get_product(**par_dict)
        disp_list.append(disp)
        parameters.append(par_dict)
except Exception as e:
    print(e)

pr.report_progress(stage="Sending requests for data products", progress=10)
if do_icecube:
    try:
        disp = DispatcherAPI(
            url="https://www.astro.unige.ch/mmoda//dispatch-data",
            instrument="mock",
            wait=False,
        )
        par_dict = {
            "RA": RA,
            "DEC": DEC,
            "IC40": False,
            "IC59": False,
            "IC79": False,
            "IC86_I": True,
            "IC86_II_VII": True,
            "Slope": 3.0,
            "Spectrum_type": "Free_slope",
            "instrument": "icecube",
            "product": "Spectrum",
            "product_type": "Real",
            "src_name": src_name,
            "token": token,
        }
        data_collection_icecube = disp.get_product(**par_dict)
        disp_list.append(disp)
        parameters.append(par_dict)
    except Exception as e:
        print(e)

pr.report_progress(stage="Sending requests for data products", progress=11)
if do_auger:
    try:
        disp = DispatcherAPI(
            url="https://www.astro.unige.ch/mmoda//dispatch-data",
            instrument="mock",
            wait=False,
        )
        par_dict = {
            "RA": RA,
            "DEC": DEC,
            "T1": T1,
            "T2": T2,
            "Emax": 3.162e20,
            "Emin": 3.162e19,
            "NEbins": 2,
            "Source_region_radius": 2.0,
            "T_format": "isot",
            "instrument": "auger",
            "product": "Spectrum",
            "product_type": "Real",
            "src_name": src_name,
            "token": token,
        }
        data_collection_auger = disp.get_product(**par_dict)
        disp_list.append(disp)
        parameters.append(par_dict)
    except Exception as e:
        print(e)

state_list = np.zeros(len(disp_list)).astype(bool)
for i, disp in enumerate(disp_list):
    print(i, disp.instrument, disp.query_status)
    if disp.query_status == "done":
        state_list[i] = True
state_list
# disp.poll()

pr.report_progress(stage="Producing instrument products", progress=12)

while True:
    instruments = ""
    for i, disp in enumerate(disp_list):
        if state_list[i] == False:
            disp.poll()
            state_list[i] = disp.is_complete
            instruments += disp.instrument + ", "
    pr.report_progress(
        stage="Working on " + instruments,
        progress=12 + int(sum(state_list) / len(state_list) * 76),
    )
    if np.all(state_list):
        break

Energies = []
instruments = []
fluxes = []
fluxes_max = []
fluxes_min = []
ULs = []
for i, disp in enumerate(disp_list):
    print(disp.instrument)
    pr.report_progress(
        stage="Collecting the results " + disp.instrument,
        progress=88 + int(10 * (i / len(disp_list))),
    )
    try:
        par_dict = parameters[i]
        data_collection = disp_list[i].get_product(**par_dict)
        if par_dict["instrument"] == "icecube":
            try:
                dic = data_collection.as_list()
                for i in range(len(dic)):
                    if dic[i]["prod_name"] == "table_0":
                        FLAG_icecube = 1
                        tab = data_collection.table_0.table
                        E_icecube = tab["Energy[TeV]"]
                        UL_icecube = tab["F_max_90[TeV/cm2s]"]
                        F_ic = tab["F_best[TeV/cm2s]"]
                        F_ic_min = tab["F_min_68[TeV/cm2s]"]
                        F_ic_max = tab["F_max_68[TeV/cm2s]"]
                        Energies.append(E_icecube)
                        fluxes.append(F_ic)
                        fluxes_max.append(F_ic_max)
                        fluxes_min.append(F_ic_min)
                        ULs.append(UL_icecube)
                        for k in range(len(E_icecube)):
                            instruments.append("IceCube")
            except Exception as e:
                print(e)

        if par_dict["instrument"] == "fermi_lat":
            try:
                dic = data_collection.as_list()
                for i in range(len(dic)):
                    if dic[i]["prod_name"] == "spectrum_astropy_table_0":
                        FLAG_fermi = 1
                        tab = data_collection.spectrum_astropy_table_0.table
                        E_fermi = tab["Emean[MeV]"] * 1e-6
                        Emin_fermi = tab["Emin[MeV]"] * 1e-6
                        Emax_fermi = tab["Emax[MeV]"] * 1e-6
                        F_fermi = tab["Flux[MeV/cm2s]"] * 1e-6
                        F_err_fermi = tab["Flux_error[MeV/cm2s]"] * 1e-6
                        UL_fermi = F_err_fermi > F_fermi / 2.0
                        Energies.append(E_fermi)
                        fluxes.append(F_fermi)
                        fluxes_max.append(F_fermi + F_err_fermi)
                        fluxes_min.append(F_fermi - F_err_fermi)
                        ULs.append(UL_fermi)
                        for k in range(len(E_fermi)):
                            instruments.append("Fermi-LAT")
                print(E_fermi)
            except Exception as e:
                print(e)

        if par_dict["instrument"] == "magic":
            try:
                dic = data_collection.as_list()
                for i in range(len(dic)):
                    if dic[i]["prod_name"] == "table_spectrum_1":
                        FLAG_magic = 1
                        tab = data_collection.table_spectrum_1.table
                        Emean_magic = tab["Emean[TeV]"]
                        Emin_magic = tab["Emin[TeV]"]
                        Emax_magic = tab["Emax[TeV]"]
                        Flux_magic = tab["Flux[TeV/cm2s]"]
                        Flux_err_magic = tab["Flux_error[TeV/cm2s]"]
                        UL_magic = Flux_err_magic > Flux_magic / 2.0
                        Energies.append(Emean_magic)
                        fluxes.append(Flux_magic)
                        fluxes_max.append(Flux_magic + Flux_err_magic)
                        fluxes_min.append(Flux_magic - Flux_err_magic)
                        ULs.append(UL_magic)
                        for k in range(len(Emean_magic)):
                            instruments.append("MAGIC")
            except Exception as e:
                print(e)

        if par_dict["instrument"] == "hess":
            try:
                dic = data_collection.as_list()
                for i in range(len(dic)):
                    if dic[i]["prod_name"] == "table_spectrum_1":
                        FLAG_hess = 1
                        tab = data_collection.table_spectrum_1.table
                        E_hess = tab["Emean[TeV]"]
                        Emin_hess = tab["Emin[TeV]"]
                        Emax_hess = tab["Emax[TeV]"]
                        F_hess = tab["Flux[TeV/cm2s]"]
                        F_err_hess = tab["Flux_error[TeV/cm2s]"]
                        UL_hess = F_err_hess > F_hess / 2.0
                        Energies.append(E_hess)
                        fluxes.append(F_hess)
                        fluxes_max.append(F_hess + F_err_hess)
                        fluxes_min.append(F_hess - F_err_hess)
                        ULs.append(UL_hess)
                        for k in range(len(E_hess)):
                            instruments.append("HESS")
            except Exception as e:
                print(e)

        if par_dict["instrument"] == "mwa":
            try:
                tab = data_collection.spectrum_image_table_1.table
                E_mwa = tab["E[eV]"] * 1e-12
                F_mwa = tab["Flux[erg/cm2s]"] / 1.6
                F_mwa_err = tab["Flux_error[erg/cm2s]"] / 1.6
                UL_mwa = F_mwa_err > F_mwa / 2.0
                FLAG_mwa = 1
                Energies.append(E_mwa)
                fluxes.append(F_mwa)
                fluxes_max.append(F_mwa + F_mwa_err)
                fluxes_min.append(F_mwa - F_mwa_err)
                ULs.append(UL_mwa)
                for k in range(len(E_mwa)):
                    instruments.append("MWA")
            except Exception as e:
                print(e)

        if par_dict["instrument"] == "desi_legacy_survey":
            try:
                tab = data_collection.spectrum_table_0.table
                E_desi = tab["Energy[eV]"] / 1e12
                F_desi = tab["Flux[erg/cm2s]"] / 1.6
                Ferr_desi = tab["Flux_err[erg/cm2s]"] / 1.6
                UL_desi = Ferr_desi > F_desi / 2.0
                FLAG_legacysurvey = 1
                Energies.append(E_desi)
                fluxes.append(F_desi)
                fluxes_max.append(F_desi + Ferr_desi)
                fluxes_min.append(F_desi - Ferr_desi)
                ULs.append(UL_desi)
                for k in range(len(E_desi)):
                    instruments.append("DESI-LegacySurvey")
            except Exception as e:
                print(e)

        pr.report_progress(stage="Collecting the results", progress=95)
        if par_dict["instrument"] == "gaia":
            try:
                tab = data_collection.spectrum_astropy_table_0.table
                E_gaia = tab["Emean[eV]"] * 1e-12
                Emin_gaia = tab["Emin[eV]"] * 1e-12
                Emax_gaia = tab["Emax[eV]"] * 1e-12
                F_gaia = tab["Flux[erg/cm2s]"] / 1.6
                Ferr_gaia = tab["Flux_error[erg/cm2s]"] / 1.6
                UL_gaia = Ferr_gaia > F_gaia / 2.0
                FLAG_gaia = 1
                Energies.append(E_gaia)
                fluxes.append(F_gaia)
                fluxes_max.append(F_gaia + Ferr_gaia)
                fluxes_min.append(F_gaia - Ferr_gaia)
                ULs.append(UL_gaia)
                for k in range(len(E_gaia)):
                    instruments.append("GAIA")
            except Exception as e:
                print(e)

        pr.report_progress(stage="Collecting the results", progress=96)
        if par_dict["instrument"] == "auger":
            try:
                dic = data_collection.as_list()
                for i in range(len(dic)):
                    if dic[i]["prod_name"] == "spectrum_astropy_table_0":
                        FLAG_auger = 1
                        tab = data_collection.spectrum_astropy_table_0.table
                        E_auger = tab["Emean[eV]"]
                        Emin_auger = tab["Emin[eV]"]
                        Emax_auger = tab["Emax[eV]"]
                        F_auger = tab["Flux[erg/cm2s]"]
                        F_err_lo_auger = tab["Flux_error_lo[erg/cm2s]"]
                        F_err_hi_auger = tab["Flux_error_hi[erg/cm2s]"]
                        UL_auger = tab["UpLim"]
                        Energies.append(E_auger)
                        fluxes.append(F_auger)
                        fluxes_max.append(F_auger + F_err_hi_auger)
                        fluxes_min.append(F_auger - F_err_lo_auger)
                        ULs.append(UL_auger)
                        for k in range(len(E_auger)):
                            instruments.append("PierreAugerObservatory")
                FLAG_auger = 1
            except Exception as e:
                print(e)

        if par_dict["instrument"] == "jemx":
            if par_dict["jemx_num"] == 1:
                try:
                    prod_list = data_collection.as_list()
                    data_collection.save_all_data()
                    for prod in prod_list:
                        if prod["meta_data:"]["src_name"] == src_name:
                            if (
                                prod["meta_data:"]["product"]
                                == "jemx_spectrum"
                            ):
                                FLAG_jemx = 1
                                fname = (
                                    workdir + "/" + prod["prod_name"] + ".fits"
                                )
                                hdul = fits.open(fname)
                                spec = hdul["JMX1-PHA1-SPE"].data
                                rate = spec["RATE"]
                                rate_err = spec["STAT_ERR"]
                                rate_sys = spec["SYS_ERR"]
                                rate_qual = spec["QUALITY"]
                            elif prod["meta_data:"]["product"] == "jemx_arf":
                                fname = (
                                    workdir + "/" + prod["prod_name"] + ".fits"
                                )
                                hdul = fits.open(fname)
                                arf = hdul["SPECRESP"].data
                                ENERG_LO = arf["ENERG_LO"]
                                ENERG_HI = arf["ENERG_HI"]
                                ENERG = sqrt(ENERG_LO * ENERG_HI)
                                dENERG = ENERG_HI - ENERG_LO
                                aeff = arf["SPECRESP"]
                            elif prod["meta_data:"]["product"] == "jemx_rmf":
                                fname = (
                                    workdir + "/" + prod["prod_name"] + ".fits"
                                )
                                hdul = fits.open(fname)
                                EBOUNDS = hdul["EBOUNDS"].data
                                rmf = hdul["SPECRESP MATRIX"].data
                                Emins = EBOUNDS["E_MIN"]
                                Emaxs = EBOUNDS["E_MAX"]
                                Emeans = sqrt(Emins * Emaxs)
                                NEtrue = len(rmf["MATRIX"])
                                NErec = len(rmf["MATRIX"][0])
                                resp_jemx = np.zeros((NEtrue, NErec))
                                for i in range(NEtrue):
                                    resp_jemx[i] = rmf["MATRIX"][i]
                                Norm = 1e-2
                                Gamma = 2.1
                                Ebins_jemx = np.logspace(0.5, 1.5, 3)
                                Emins_jemx = Ebins_jemx[:-1]
                                Emaxs_jemx = Ebins_jemx[1:]
                                Emeans_jemx = sqrt(Emins_jemx * Emaxs_jemx)
                                F_jemx = np.zeros(len(Emeans_jemx))
                                F_jemx_err = np.zeros(len(Emeans_jemx))
                                for i in range(len(Emins_jemx)):
                                    m = (Emeans > Emins_jemx[i]) & (
                                        Emeans < Emaxs_jemx[i]
                                    )
                                    model_cts = sum(
                                        m * exp_counts_jemx(Norm, Gamma)
                                    )
                                    real_cts = np.nansum(m * rate)
                                    print(model_cts, real_cts)
                                    real_err = sqrt(np.nansum(m * rate_err**2))
                                    ratio = real_cts / model_cts
                                    F_jemx[i] = (
                                        ratio
                                        * dn_de(Emeans_jemx[i], Norm, Gamma)
                                        * Emeans_jemx[i] ** 2
                                        * 1e-9
                                    )
                                    F_jemx_err[i] = (
                                        F_jemx[i] / real_cts * real_err
                                    )
                except Exception as e:
                    print(e)
            else:
                try:
                    prod_list = data_collection.as_list()
                    data_collection.save_all_data()
                    for prod in prod_list:
                        if prod["meta_data:"]["src_name"] == src_name:
                            FLAG_jemx = 1
                            if (
                                prod["meta_data:"]["product"]
                                == "jemx_spectrum"
                            ):
                                fname = (
                                    workdir + "/" + prod["prod_name"] + ".fits"
                                )
                                hdul = fits.open(fname)
                                spec = hdul["JMX2-PHA1-SPE"].data
                                rate = spec["RATE"]
                                rate_err = spec["STAT_ERR"]
                                rate_sys = spec["SYS_ERR"]
                                rate_qual = spec["QUALITY"]
                            elif prod["meta_data:"]["product"] == "jemx_arf":
                                fname = (
                                    workdir + "/" + prod["prod_name"] + ".fits"
                                )
                                hdul = fits.open(fname)
                                arf = hdul["SPECRESP"].data
                                ENERG_LO = arf["ENERG_LO"]
                                ENERG_HI = arf["ENERG_HI"]
                                ENERG = sqrt(ENERG_LO * ENERG_HI)
                                dENERG = ENERG_HI - ENERG_LO
                                aeff = arf["SPECRESP"]
                            elif prod["meta_data:"]["product"] == "jemx_rmf":
                                fname = (
                                    workdir + "/" + prod["prod_name"] + ".fits"
                                )
                                hdul = fits.open(fname)
                                EBOUNDS = hdul["EBOUNDS"].data
                                rmf = hdul["SPECRESP MATRIX"].data
                                Emins = EBOUNDS["E_MIN"]
                                Emaxs = EBOUNDS["E_MAX"]
                                Emeans = sqrt(Emins * Emaxs)
                                NEtrue = len(rmf["MATRIX"])
                                NErec = len(rmf["MATRIX"][0])
                                resp_jemx = np.zeros((NEtrue, NErec))
                                for i in range(NEtrue):
                                    resp_jemx[i] = rmf["MATRIX"][i]
                                Norm = 1e-2
                                Gamma = 2.1
                                Ebins_jemx = np.logspace(0.5, 1.5, 3)
                                Emins_jemx = Ebins_jemx[:-1]
                                Emaxs_jemx = Ebins_jemx[1:]
                                Emeans_jemx = sqrt(Emins_jemx * Emaxs_jemx)
                                F_jemx2 = np.zeros(len(Emeans_jemx))
                                F_jemx2_err = np.zeros(len(Emeans_jemx))
                                for i in range(len(Emins_jemx)):
                                    m = (Emeans > Emins_jemx[i]) & (
                                        Emeans < Emaxs_jemx[i]
                                    )
                                    model_cts = sum(
                                        m * exp_counts_jemx(Norm, Gamma)
                                    )
                                    real_cts = np.nansum(m * rate)
                                    print(model_cts, real_cts)
                                    real_err = sqrt(np.nansum(m * rate_err**2))
                                    ratio = real_cts / model_cts
                                    F_jemx2[i] = (
                                        ratio
                                        * dn_de(Emeans_jemx[i], Norm, Gamma)
                                        * Emeans_jemx[i] ** 2
                                        * 1e-9
                                    )
                                    F_jemx2_err[i] = (
                                        F_jemx2[i] / real_cts * real_err
                                    )
                except Exception as e:
                    print(e)

        if par_dict["instrument"] == "isgri":
            try:
                prod_list = data_collection.as_list()
                data_collection.save_all_data()
                for prod in prod_list:
                    if prod["meta_data:"]["src_name"] == src_name:
                        if prod["meta_data:"]["product"] == "isgri_spectrum":
                            FLAG_isgri = 1
                            fname = workdir + "/" + prod["prod_name"] + ".fits"
                            hdul = fits.open(fname)
                            spec = hdul["ISGR-EVTS-SPE"].data
                            rate = spec["RATE"]
                            rate_err = spec["STAT_ERR"]
                            rate_sys = spec["SYS_ERR"]
                            rate_qual = spec["QUALITY"]
                        elif prod["meta_data:"]["product"] == "isgri_arf":
                            fname = workdir + "/" + prod["prod_name"] + ".fits"
                            hdul = fits.open(fname)
                            arf = hdul["SPECRESP"].data
                            ENERG_LO = arf["ENERG_LO"]
                            ENERG_HI = arf["ENERG_HI"]
                            ENERG = sqrt(ENERG_LO * ENERG_HI)
                            dENERG = ENERG_HI - ENERG_LO
                        elif prod["meta_data:"]["product"] == "isgri_rmf":
                            fname = workdir + "/" + prod["prod_name"] + ".fits"
                            hdul = fits.open(fname)
                            EBOUNDS = hdul["EBOUNDS"].data
                            rmf = hdul["SPECRESP MATRIX"].data
                            Emins = EBOUNDS["E_MIN"]
                            Emaxs = EBOUNDS["E_MAX"]
                            Emeans = sqrt(Emins * Emaxs)
                            NEtrue = len(rmf["MATRIX"])
                            NErec = len(rmf["MATRIX"][0])
                            resp_isgri = np.zeros((NEtrue, NErec))
                            for i in range(NEtrue):
                                resp_isgri[i] = rmf["MATRIX"][i]
                            Norm = 1e-2
                            Gamma = 2.1
                            Ebins_isgri = np.logspace(1.5, 2.5, 3)
                            Emins_isgri = Ebins_isgri[:-1]
                            Emaxs_isgri = Ebins_isgri[1:]
                            Emeans_isgri = sqrt(Emins_isgri * Emaxs_isgri)
                            F_isgri = np.zeros(len(Emeans_isgri))
                            F_isgri_err = np.zeros(len(Emeans_isgri))
                            for i in range(len(Emins_isgri)):
                                m = (Emeans > Emins_isgri[i]) & (
                                    Emeans < Emaxs_isgri[i]
                                )
                                model_cts = sum(
                                    m * exp_counts_isgri(Norm, Gamma)
                                )
                                real_cts = np.nansum(m * rate)
                                real_err = sqrt(np.nansum(m * rate_err**2))
                                ratio = real_cts / model_cts
                                F_isgri[i] = (
                                    ratio
                                    * dn_de(Emeans_isgri[i], Norm, Gamma)
                                    * Emeans_isgri[i] ** 2
                                    * 1e-9
                                )
                                F_isgri_err[i] = (
                                    F_isgri[i] / real_cts * real_err
                                )

                            E_isgri = Emeans_isgri * 1e-9
                            Emins_isgri = Emins_isgri * 1e-9
                            Emaxs_isgri = Emaxs_isgri * 1e-9
                            UL_isgri = F_isgri_err > F_isgri / 2.0
                            Energies.append(E_isgri)
                            fluxes.append(F_isgri)
                            fluxes_max.append(F_isgri + F_isgri_err)
                            fluxes_min.append(F_isgri - F_isgri_err)
                            ULs.append(UL_isgri)
                            for k in range(len(E_isgri)):
                                instruments.append("INTEGRAL-ISGRI")

            except Exception as e:
                print(e)
    except Exception as e:
        print(e)

if FLAG_jemx == 1:
    tmp = F_jemx2 / F_jemx2_err**2 + F_jemx / F_jemx_err**2
    tmp1 = 1 / F_jemx2_err**2 + 1 / F_jemx_err**2
    F_jemx_all = tmp / tmp1
    F_jemx_all_err = 1 / sqrt(tmp1)
    F_jemx_all, F_jemx_all_err
    UL_jemx = F_jemx_all_err > F_jemx_all / 2.0
    Energies.append(Emeans_jemx)
    fluxes.append(F_jemx_all)
    fluxes_max.append(F_jemx_all + F_jemx_all_err)
    fluxes_min.append(F_jemx_all - F_jemx_all_err)
    ULs.append(UL_jemx)
    for k in range(len(Emeans_jemx)):
        instruments.append("INTEGRAL-JEMX")

E = np.array([])
F = np.array([])
F_max = np.array([])
F_min = np.array([])
UL = np.array([])
for i, EE in enumerate(Energies):
    try:
        E = np.concatenate((E, EE.value))
    except:
        E = np.concatenate((E, EE))
    try:
        F = np.concatenate((F, fluxes[i].value))
    except:
        F = np.concatenate((F, fluxes[i]))
    try:
        F_max = np.concatenate((F_max, fluxes_max[i].value))
    except:
        F_max = np.concatenate((F_max, fluxes_max[i]))
    try:
        F_min = np.concatenate((F_min, fluxes_min[i].value))
    except:
        F_min = np.concatenate((F_min, fluxes_min[i]))
    try:
        UL = np.concatenate((UL, ULs[i].value))
    except:
        UL = np.concatenate((F_min, ULs[i]))

len(instruments), len(F), len(E)
from astropy.table import Table
from oda_api.data_products import ODAAstropyTable

data = [E, F, F_max, F_min, instruments]
names = (
    "Energy[TeV]",
    "F[TeV/cm2s]",
    "F_max[TeV/cm2s]",
    "F_min[TeV/cm2s]",
    "instrument",
)
sed_table = ODAAstropyTable(Table(data, names=names))

if FLAG_mwa > 0:
    plt.errorbar(E_mwa, F_mwa, yerr=F_mwa_err, label="MWA")
    ymin_mwa = min(F_mwa / 3.0)
    ymax_mwa = max(F_mwa * 3.0)

pr.report_progress(stage="SED preparation", progress=99)
ymin_auger = 1e20
ymax_auger = 1e-20
ymin_magic = 1e20
ymax_magic = 1e-20
ymin_hess = 1e20
ymax_hess = 1e-20
ymin_ic = 1e20
ymax_ic = 1e-20
ymin_fermi = 1e20
ymax_fermi = 1e-20
ymin_isgri = 1e20
ymax_isgri = 1e-20
ymin_mwa = 1e20
ymax_mwa = 1e-20

plt.figure(figsize=(12, 7))
if FLAG_legacysurvey == 1:
    plt.errorbar(
        E_desi,
        F_desi,
        yerr=Ferr_desi,
        uplims=UL_desi,
        linestyle="none",
        color="black",
        label="DESI Legacy Survey",
        marker="o",
    )
if FLAG_gaia == 1:
    plt.errorbar(
        E_gaia,
        F_gaia,
        xerr=[E_gaia - Emin_gaia, Emax_gaia - E_gaia],
        yerr=Ferr_gaia,
        uplims=UL_gaia,
        linestyle="none",
        color="magenta",
        label="GAIA",
        marker="*",
    )
if FLAG_isgri > 0:
    m = (E_isgri > 3e-8) & (E_isgri < 2e-7)
    plt.errorbar(
        E_isgri[m],
        F_isgri[m],
        yerr=F_isgri_err[m],
        xerr=[E_isgri - Emins_isgri, Emaxs_isgri - E_isgri],
        uplims=UL_isgri,
        label="INTEGRAL/ISGRI",
    )
    ymin_isgri = min(F_isgri[m]) / 3.0
    ymax_isgri = max(F_isgri[m]) * 3.0
if FLAG_jemx > 0:
    plt.errorbar(
        Emeans_jemx * 1e-9,
        F_jemx,
        yerr=F_jemx_err,
        xerr=[
            (Emeans_jemx - Emins_jemx) * 1e-9,
            (Emaxs_jemx - Emeans_jemx) * 1e-9,
        ],
        label="INTEGRAL/JEMX",
    )
    ymin_isgri = min(F_isgri[m]) / 3.0
    ymax_isgri = max(F_isgri[m]) * 3.0

if FLAG_auger > 0:
    plt.errorbar(
        E_auger,
        F_auger,
        xerr=[E_auger - Emin_auger, Emax_auger - E_auger],
        yerr=[F_err_lo_auger, F_err_hi_auger],
        uplims=UL_auger,
        linestyle="none",
        color="green",
        label="Auger",
    )
    ymin_auger = min(F_auger) / 3.0
    ymax_auger = max(F_auger) * 3.0

if FLAG_icecube > 0:
    if min(F_ic_min) > 1e-14:
        plt.fill_between(
            E_icecube, F_ic_min, F_ic_max, alpha=0.1, color="black"
        )
    plt.plot(E_icecube, UL_icecube, color="black", label="IceCube")
    E_icecube_line = E_icecube[np.argmin(UL_icecube)]
    F_icecube_line = UL_icecube[np.argmin(UL_icecube)]
    plt.annotate(
        "",
        xy=(E_icecube_line, F_icecube_line / 2.0),
        xytext=(E_icecube_line, F_icecube_line),
        arrowprops=dict(facecolor="black", shrink=0.0, width=1, headwidth=4),
    )
    ymin_ic = F_icecube_line / 3.0
    ymax_ic = max(F_ic_max * 3.0)

if FLAG_fermi > 0:
    plt.errorbar(
        E_fermi,
        F_fermi,
        yerr=F_err_fermi,
        xerr=[E_fermi - Emin_fermi, Emax_fermi - E_fermi],
        label="Fermi/LAT",
        uplims=UL_fermi,
    )
    ymin_fermi = min(F_fermi / 3.0)
    ymax_fermi = max(F_fermi * 3.0)

if FLAG_magic > 0:
    plt.errorbar(Emean_magic, Flux_magic, yerr=Flux_err_magic, label="MAGIC")
    ymin_magic = min(Flux_magic / 3.0)
    ymax_magic = max(Flux_magic * 3.0)

if FLAG_hess > 0:
    plt.errorbar(E_hess, F_hess, yerr=F_err_hess, label="HESS")
    ymin_hess = min(F_hess / 3.0)
    ymax_hess = max(F_hess * 3.0)

if FLAG_mwa > 0:
    plt.errorbar(E_mwa, F_mwa, yerr=F_mwa_err, label="MWA", marker="o")
    ymin_mwa = min(F_mwa / 3.0)
    ymax_mwa = max(F_mwa * 3.0)

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Energy [TeV]")
plt.ylabel("Flux [TeV/(cm$^2$s)]")
plt.legend(loc="upper right")
ymin = min([ymin_mwa, ymin_isgri, ymin_ic, ymin_fermi, ymin_magic, ymin_auger])
ymax = max([ymax_mwa, ymax_isgri, ymax_ic, ymax_fermi, ymax_magic, ymax_auger])
plt.ylim(ymin, ymax)
plt.title(src_name)
plt.savefig("SED.png", format="png", bbox_inches="tight")

from oda_api.data_products import PictureProduct

bin_image = PictureProduct.from_file("SED.png")

sed_png = bin_image  # http://odahub.io/ontology#ODAPictureProduct
sed_table = sed_table  # http://odahub.io/ontology#ODAAstropyTable

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    ("out_Multi_instrument_spectrum_sed_png", "sed_png_galaxy.output", sed_png)
)
_oda_outs.append(
    (
        "out_Multi_instrument_spectrum_sed_table",
        "sed_table_galaxy.output",
        sed_table,
    )
)

for _outn, _outfn, _outv in _oda_outs:
    _galaxy_outfile_name = os.path.join(_galaxy_wd, _outfn)
    if isinstance(_outv, str) and os.path.isfile(_outv):
        shutil.move(_outv, _galaxy_outfile_name)
        _galaxy_meta_data[_outn] = {"ext": "_sniff_"}
    elif getattr(_outv, "write_fits_file", None):
        _outv.write_fits_file(_galaxy_outfile_name)
        _galaxy_meta_data[_outn] = {"ext": "fits"}
    elif getattr(_outv, "write_file", None):
        _outv.write_file(_galaxy_outfile_name)
        _galaxy_meta_data[_outn] = {"ext": "_sniff_"}
    else:
        with open(_galaxy_outfile_name, "w") as fd:
            json.dump(_outv, fd, cls=CustomJSONEncoder)
        _galaxy_meta_data[_outn] = {"ext": "json"}

with open(os.path.join(_galaxy_wd, "galaxy.json"), "w") as fd:
    json.dump(_galaxy_meta_data, fd)
print("*** Job finished successfully ***")
