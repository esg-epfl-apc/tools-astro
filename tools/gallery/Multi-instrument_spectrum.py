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

workdir = os.getcwd()

src_name = "Mrk 421"  # http://odahub.io/ontology#AstrophysicalObject
RA = 166.113809  # http://odahub.io/ontology#PointOfInterestRA
DEC = 38.208833  # http://odahub.io/ontology#PointOfInterestDEC
T1 = "2000-10-09T13:16:00.0"  # http://odahub.io/ontology#StartTime
T2 = "2021-10-13T13:16:00.0"  # http://odahub.io/ontology#EndTime

do_isgri = (
    True  # http://odahub.io/ontology#Boolean ; oda:label "INTEGRAL/ISGRI"
)
do_fermi = True  # http://odahub.io/ontology#Boolean ; oda:label "Fermi/LAT"
do_magic = True  # http://odahub.io/ontology#Boolean ; oda:label "MAGIC"
do_icecube = True  # http://odahub.io/ontology#Boolean ; oda:label "IceCube"
do_auger = True  # http://odahub.io/ontology#Boolean ; oda:label "Auger"
do_hess = True  # http://odahub.io/ontology#Boolean ; oda:label "HESS"

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
    "do_isgri",
    "do_fermi",
    "do_magic",
    "do_icecube",
    "do_auger",
    "do_hess",
]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

disp = DispatcherAPI(
    url="https://www.astro.unige.ch/mmoda//dispatch-data", instrument="mock"
)
token = discover_token()

E0 = 20.0  # keV

def dn_de(E, A, Gam):
    return A * (E / E0) ** (-Gam)

def exp_counts(A, Gam):
    spectrum_dndE = dn_de(ENERG, A, Gam)
    tmp = np.zeros(NErec)
    for i in range(NErec):
        tmp[i] += sum(spectrum_dndE * dENERG * resp[:, i])
    return tmp

pr.report_progress(stage="INTEGRAL/ISGRI", progress=10)
FLAG_isgri = 0
if do_isgri:
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

if do_isgri:
    prod_list = data_collection_isgri.as_list()
    for prod in prod_list:
        if prod["meta_data:"]["src_name"] == src_name:
            if prod["meta_data:"]["product"] == "isgri_spectrum":
                FLAG_isgri = 1
                print(prod["meta_data:"]["product"])
                fname = workdir + "/" + prod["prod_name"] + ".fits"
                hdul = fits.open(fname)
                spec = hdul["ISGR-EVTS-SPE"].data
                rate = spec["RATE"]
                rate_err = spec["STAT_ERR"]
                rate_sys = spec["SYS_ERR"]
                rate_qual = spec["QUALITY"]
            elif prod["meta_data:"]["product"] == "isgri_arf":
                print(prod["meta_data:"]["product"])
                fname = workdir + "/" + prod["prod_name"] + ".fits"
                hdul = fits.open(fname)
                arf = hdul["SPECRESP"].data
                ENERG_LO = arf["ENERG_LO"]
                ENERG_HI = arf["ENERG_HI"]
                ENERG = sqrt(ENERG_LO * ENERG_HI)
                dENERG = ENERG_HI - ENERG_LO
            elif prod["meta_data:"]["product"] == "isgri_rmf":
                print(prod["meta_data:"]["product"])
                fname = workdir + "/" + prod["prod_name"] + ".fits"
                hdul = fits.open(fname)
                EBOUNDS = hdul["EBOUNDS"].data
                rmf = hdul["SPECRESP MATRIX"].data
                Emins = EBOUNDS["E_MIN"]
                Emaxs = EBOUNDS["E_MAX"]
                Emeans = sqrt(Emins * Emaxs)
                NEtrue = len(rmf["MATRIX"])
                NErec = len(rmf["MATRIX"][0])
                resp = np.zeros((NEtrue, NErec))
                for i in range(NEtrue):
                    resp[i] = rmf["MATRIX"][i]
                Norm = 1e-2
                Gamma = 2.1
                Ebins_isgri = np.logspace(1.5, 2.5, 3)
                Emins_isgri = Ebins_isgri[:-1]
                Emaxs_isgri = Ebins_isgri[1:]
                Emeans_isgri = sqrt(Emins_isgri * Emaxs_isgri)
                F_isgri = np.zeros(len(Emeans_isgri))
                F_isgri_err = np.zeros(len(Emeans_isgri))
                for i in range(len(Emins_isgri)):
                    m = (Emeans > Emins_isgri[i]) & (Emeans < Emaxs_isgri[i])
                    model_cts = sum(m * exp_counts(Norm, Gamma))
                    real_cts = np.nansum(m * rate)
                    real_err = sqrt(np.nansum(m * rate_err**2))
                    ratio = real_cts / model_cts
                    F_isgri[i] = (
                        ratio
                        * dn_de(Emeans_isgri[i], Norm, Gamma)
                        * Emeans_isgri[i] ** 2
                        * 1e-9
                    )
                    F_isgri_err[i] = F_isgri[i] / real_cts * real_err

                E_isgri = Emeans_isgri * 1e-9
                Emins_isgri = Emins_isgri * 1e-9
                Emaxs_isgri = Emaxs_isgri * 1e-9

pr.report_progress(stage="HESS", progress=30)
FLAG_hess = 0
if do_hess:
    try:
        par_dict = {
            "RA": RA,
            "DEC": DEC,
            "Efit_max": 10.0,
            "Efit_min": 0.2,
            "Emax": 100.0,
            "Emin": 1.0,
            "NEbins": 30,
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
        dic = data_collection_hess.as_list()
        for i in range(len(dic)):
            if dic[i]["prod_name"] == "table_spectrum_1":
                FLAG_hess = 1
                tab = data_collection_hess.table_spectrum_1.table
                E_hess = tab["Emean[TeV]"]
                Emin_hess = tab["Emin[TeV]"]
                Emax_hess = tab["Emax[TeV]"]
                F_hess = tab["Flux[TeV/cm2s]"]
                F_err_hess = tab["Flux_error[TeV/cm2s]"]
    except:
        print("No HESS data")

pr.report_progress(stage="MAGIC", progress=40)
FLAG_magic = 0
if do_magic:
    try:
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
        dic = data_collection_magic.as_list()
        for i in range(len(dic)):
            if dic[i]["prod_name"] == "table_spectrum_1":
                FLAG_magic = 1
                tab = data_collection_magic.table_spectrum_1.table
                Emean_magic = tab["Emean[TeV]"]
                Emin_magic = tab["Emin[TeV]"]
                Emax_magic = tab["Emax[TeV]"]
                Flux_magic = tab["Flux[TeV/cm2s]"]
                Flux_err_magic = tab["Flux_error[TeV/cm2s]"]
    except:
        print("No MAGIC data")

pr.report_progress(stage="Fermi/LAT", progress=50)
FLAG_fermi = 0
if do_fermi:
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
        "src_name": "Crab",
        "token": token,
    }
    data_collection_fermi = disp.get_product(**par_dict)
    dic = data_collection_fermi.as_list()
    for i in range(len(dic)):
        if dic[i]["prod_name"] == "spectrum_astropy_table_0":
            FLAG_fermi = 1
            tab = data_collection_fermi.spectrum_astropy_table_0.table
            E_fermi = tab["Emean[MeV]"] * 1e-6
            Emin_fermi = tab["Emin[MeV]"] * 1e-6
            Emax_fermi = tab["Emax[MeV]"] * 1e-6
            F_fermi = tab["Flux[MeV/cm2s]"] * 1e-6
            F_err_fermi = tab["Flux_error[MeV/cm2s]"] * 1e-6

pr.report_progress(stage="IceCube", progress=80)
FLAG_icecube = 0
if do_icecube:
    try:
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
            "src_name": "Crab",
            "token": token,
        }
        data_collection_icecube = disp.get_product(**par_dict)
        dic = data_collection_icecube.as_list()
        for i in range(len(dic)):
            if dic[i]["prod_name"] == "table_0":
                FLAG_icecube = 1
                tab = data_collection_icecube.table_0.table
                E_icecube = tab["Energy[TeV]"]
                UL_icecube = tab["F_max_90[TeV/cm2s]"]
                F_ic = tab["F_best[TeV/cm2s]"]
                F_ic_min = tab["F_min_68[TeV/cm2s]"]
                F_ic_max = tab["F_max_68[TeV/cm2s]"]
    except:
        print("No IceCube data")

pr.report_progress(stage="Auger", progress=90)
FLAG_auger = 0
if do_auger:
    try:
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
            "src_name": "Crab",
            "token": token,
        }
        data_collection_auger = disp.get_product(**par_dict)
        dic = data_collection_auger.as_list()
        for i in range(len(dic)):
            if dic[i]["prod_name"] == "spectrum_astropy_table_0":
                FLAG_auger = 1
                tab = data_collection_auger.spectrum_astropy_table_0.table
                E_auger = tab["Emean[eV]"]
                Emin_auger = tab["Emin[eV]"]
                Emax_auger = tab["Emax[eV]"]
                F_auger = tab["Flux[erg/cm2s]"]
                F_err_lo_auger = tab["Flux_error_lo[erg/cm2s]"]
                F_err_hi_auger = tab["Flux_error_hi[erg/cm2s]"]
                UL_auger = tab["UpLim"]
    except:
        print("No Auger data")

pr.report_progress(stage="Data product preparation", progress=95)
ymin_auger = 1e20
ymax_auger = 1e-20
ymin_magic = 1e20
ymax_magic = 1e-20
ymin_ic = 1e20
ymax_ic = 1e-20
ymin_fermi = 1e20
ymax_fermi = 1e-20
ymin_isgri = 1e20
ymax_isgri = 1e-20

if FLAG_isgri > 0:
    m = (E_isgri > 3e-8) & (E_isgri < 2e-7)
    plt.errorbar(
        E_isgri[m],
        F_isgri[m],
        yerr=F_isgri_err[m],
        xerr=[E_isgri - Emins_isgri, Emaxs_isgri - E_isgri],
        label="INTEGRAL/ISGRI",
    )
    ymin_isgri = min(F_isgri[m]) / 10.0
    ymax_isgri = max(F_isgri[m]) * 10.0

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
    ymin_auger = min(F_auger) / 10.0
    ymax_auger = max(F_auger) * 10.0

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
    ymin_ic = F_icecube_line / 10.0
    ymax_ic = max(F_ic_max * 10.0)

if FLAG_fermi > 0:
    plt.errorbar(
        E_fermi,
        F_fermi,
        yerr=F_err_fermi,
        xerr=[E_fermi - Emin_fermi, Emax_fermi - E_fermi],
        label="Fermi/LAT",
    )
    ymin_fermi = min(F_fermi / 10.0)
    ymax_fermi = max(F_fermi * 10.0)

if FLAG_magic > 0:
    plt.errorbar(Emean_magic, Flux_magic, yerr=Flux_err_magic, label="MAGIC")
    ymin_magic = min(Flux_magic / 10.0)
    ymax_magic = max(Flux_magic * 10.0)

if FLAG_hess > 0:
    plt.errorbar(E_hess, F_hess, yerr=F_err_hess, label="HESS")
    ymin_magic = min(Flux_magic / 10.0)
    ymax_magic = max(Flux_magic * 10.0)

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Energy [TeV]")
plt.ylabel("Flux [TeV/(cm$^2$s)]")
plt.legend(loc="upper right")
ymin = min([ymin_isgri, ymin_ic, ymin_fermi, ymin_magic, ymin_auger])
ymax = max([ymax_isgri, ymax_ic, ymax_fermi, ymax_magic, ymax_auger])
plt.ylim(ymin, ymax)
plt.title(src_name)
plt.savefig("SED.png", format="png", bbox_inches="tight")

from oda_api.data_products import PictureProduct

bin_image = PictureProduct.from_file("SED.png")

sed_png = bin_image  # http://odahub.io/ontology#ODAPictureProduct

FLAG_isgri = 1

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    ("out_Multi_instrument_spectrum_sed_png", "sed_png_galaxy.output", sed_png)
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
