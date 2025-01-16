#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os

import matplotlib.pyplot as plt
import numpy as np
from oda_api.api import DispatcherAPI
from oda_api.token import discover_token

src_name = "Crab"  # http://odahub.io/ontology#AstrophysicalObject
RA = 83.628700  # http://odahub.io/ontology#PointOfInterestRA
DEC = 22.014700  # http://odahub.io/ontology#PointOfInterestDEC

T1 = "2000-10-09T13:16:00.0"  # http://odahub.io/ontology#StartTime
T2 = "2021-10-13T13:16:00.0"  # http://odahub.io/ontology#EndTime

do_fermi = True  # http://odahub.io/ontology#Boolean ; oda:label "Fermi/LAT"
do_magic = True  # http://odahub.io/ontology#Boolean ; oda:label "MAGIC"
do_icecube = True  # http://odahub.io/ontology#Boolean ; oda:label "IceCube"
do_auger = True  # http://odahub.io/ontology#Boolean ; oda:label "Auger"

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
    "do_fermi",
    "do_magic",
    "do_icecube",
    "do_auger",
]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

disp = DispatcherAPI(
    url="https://www.astro.unige.ch/mmoda//dispatch-data", instrument="mock"
)

token = discover_token()

FLAG_magic = 0
if do_magic:
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
        "T2": "2017-05-06T15:32:27.000",
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

FLAG_icecube = 0
if do_icecube:
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

FLAG_auger = 0
if do_auger:
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

ymin_auger = 1e20
ymax_auger = 1e-20
ymin_magic = 1e20
ymax_magic = 1e-20
ymin_ic = 1e20
ymax_ic = 1e-20
ymin_fermi = 1e20
ymax_fermi = 1e-20

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
    ymax_auger = max(F_auger) / 10.0

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

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Energy [TeV]")
plt.ylabel("Flux [TeV/(cm$^2$s)]")
plt.legend(loc="upper right")
ymin = min([ymin_ic, ymin_fermi, ymin_magic, ymin_auger])
ymax = max([ymax_ic, ymax_fermi, ymax_magic, ymax_auger])
plt.ylim(ymin, ymax)
plt.title(src_name)
plt.savefig("SED.png", format="png", bbox_inches="tight")

from oda_api.data_products import PictureProduct

bin_image = PictureProduct.from_file("SED.png")

sed_png = bin_image  # http://odahub.io/ontology#ODAPictureProduct

# output gathering
_galaxy_meta_data = {}

with open(os.path.join(_galaxy_wd, "galaxy.json"), "w") as fd:
    json.dump(_galaxy_meta_data, fd)
print("*** Job finished successfully ***")
