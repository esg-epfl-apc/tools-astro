#!/usr/bin/env python
# coding: utf-8

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

# flake8: noqa

import json
import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt
from oda_api.api import DispatcherAPI
from oda_api.json import CustomJSONEncoder

disp = DispatcherAPI(
    url="https://www.astro.unige.ch/mmoda//dispatch-data", instrument="mock"
)

redshift = 0.9
token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJhbmRyaWkubmVyb25vdkBnbWFpbC5jb20iLCJlbWFpbCI6ImFuZHJpaS5uZXJvbm92QGdtYWlsLmNvbSIsIm5hbWUiOiJhbmRyaWluZXJvbm92Iiwicm9sZXMiOiJhdXRoZW50aWNhdGVkIHVzZXIsIGFkbWluaXN0cmF0b3IsIHVzZXIgbWFuYWdlciwgZ2VuZXJhbCwgaW50ZWdyYWwtcHJpdmF0ZS1xbGEsIHVuaWdlLWhwYy1mdWxsLCBwdWJsaWMtcG9vbC1ocGMsIGFudGFyZXMsIHNkc3MsIGFwYywgcmVua3UgY29udHJpYnV0b3IsIGdhbGxlcnkgY29udHJpYnV0b3IsIG9kYSB3b3JrZmxvdyBkZXZlbG9wZXIiLCJleHAiOjE3Mzc3MzQ3MjV9.uXgzftZwf55IJGN6VGuO8ztTccFU-XitWtg11fvDocI"

T = 2.73 * (1 + redshift)

B = 30e-9
Ecut = 4e10 / sqrt(2)

par_dict = {
    "B": B,
    "Back_norm": 1.0,
    "Emax": 100000000000000.0,
    "Emin": 1e-10,
    "T": T,
    "Z": 1.4,
    "backgr_dN_dE": "",
    "backgr_file": "",
    "dN_dE": "2.0e-11*pow(E/5e10, -2.)*exp(-E/" + str(Ecut) + ")",
    "electron_file": "",
    "instrument": "synch_ic_brems",
    "n": 0.1,
    "product": "Synchrotron_external-Compton_Bremsstrahlung",
    "product_type": "Real",
    "token": token,
}

data_collection = disp.get_product(**par_dict)

E = data_collection.spectrum_table_0.table["E[eV]"]
synch_30nG = data_collection.spectrum_table_0.table["E2dNdE_synch"]
ics_30nG = data_collection.spectrum_table_0.table["E2dNdE_ics"]

B = 60e-9
Ecut = 4e10 / 2
par_dict = {
    "B": B,
    "Back_norm": 1.0,
    "Emax": 100000000000000.0,
    "Emin": 1e-10,
    "T": T,
    "Z": 1.4,
    "backgr_dN_dE": "",
    "backgr_file": "",
    "dN_dE": "2.0e-11*pow(E/5e10, -2.)*exp(-E/" + str(Ecut) + ")",
    "electron_file": "",
    "instrument": "synch_ic_brems",
    "n": 0.1,
    "product": "Synchrotron_external-Compton_Bremsstrahlung",
    "product_type": "Real",
    "token": token,
}

data_collection1 = disp.get_product(**par_dict)

E = data_collection1.spectrum_table_0.table["E[eV]"]
synch_60nG = data_collection1.spectrum_table_0.table["E2dNdE_synch"]
ics_60nG = data_collection1.spectrum_table_0.table["E2dNdE_ics"]

B = 15e-9
Ecut = 4e10

par_dict = {
    "B": B,
    "Back_norm": 1.0,
    "Emax": 100000000000000.0,
    "Emin": 1e-10,
    "T": T,
    "Z": 1.4,
    "backgr_dN_dE": "",
    "backgr_file": "",
    "dN_dE": "2.0e-11*pow(E/5e10, -2.)*exp(-E/" + str(Ecut) + ")",
    "electron_file": "",
    "instrument": "synch_ic_brems",
    "n": 0.1,
    "product": "Synchrotron_external-Compton_Bremsstrahlung",
    "product_type": "Real",
    "token": token,
}

data_collection2 = disp.get_product(**par_dict)

E = data_collection2.spectrum_table_0.table["E[eV]"]
synch_15nG = data_collection2.spectrum_table_0.table["E2dNdE_synch"]
ics_15nG = data_collection2.spectrum_table_0.table["E2dNdE_ics"]

factor = 3e-11
plt.plot(
    E / (1 + redshift), synch_15nG * factor, color="black", linestyle="dashed"
)
plt.plot(
    E / (1 + redshift),
    ics_15nG * factor,
    color="black",
    linestyle="dashed",
    label=r"$B=15$ nG",
)

factor = 3e-11 / 4.0
plt.plot(E / (1 + redshift), synch_30nG * factor, color="black", linewidth=3)
plt.plot(
    E / (1 + redshift),
    ics_30nG * factor,
    color="black",
    linestyle="solid",
    label=r"$B=30$ nG",
    linewidth=3,
)

factor = 3e-11 / 16.0
# plt.plot(E/(1+redshift),synch_60nG*factor,color='black')
plt.plot(
    E / (1 + redshift), ics_60nG * factor, color="black", label=r"$B=60$ nG"
)

plt.axhline(2.5e-11, color="blue", linestyle="dotted")
plt.text(1e-8, 3e-11, "Eddington limit", color="blue")
plt.text(1e8, 5e-11, "Fermi/LAT", color="black")
plt.text(3e-7, 1.5e-16, "LOFAR", color="red")
plt.text(3e-2, 0.4e-12, "AGN core", color="blue", alpha=0.5)
plt.text(0.7e-8, 1.6e-17, "synchrotron", color="black")
plt.text(1e3, 1.6e-14, "inverse Compton", color="black")

d = np.genfromtxt("data_Porphyrion/porphyrion_LOFAR.csv")
plt.scatter(d[:, 0], d[:, 1], color="red", marker="o")
d = np.genfromtxt("data_Porphyrion/porphyrion_fermi.csv")
plt.plot(d[:, 0], d[:, 1], color="black", linewidth=2)
plt.fill_between(
    d[:, 0], d[:, 1], 0.0 * d[:, 0] + 1, color="black", linewidth=0, alpha=0.25
)
d = np.genfromtxt("data_Porphyrion/porphyrion_core.csv")
plt.plot(d[:, 0], d[:, 1], color="blue", linewidth=4, alpha=0.3)

plt.xscale("log")
plt.yscale("log")

plt.ylim(1.3e-17, 1e-10)
plt.xlim(1e-9, 1e11)
plt.legend(loc="lower right")
plt.xlabel("Energy, eV")
plt.ylabel("Flux, erg/(cm$^2$s)")
plt.savefig("Spectrum.png", format="png", bbox_inches="tight")

from oda_api.data_products import PictureProduct

bin_image = PictureProduct.from_file("Spectrum.png")

png = bin_image  # http://odahub.io/ontology#ODAPictureProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(("out_Porphyrion_png", "png_galaxy.output", png))

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
