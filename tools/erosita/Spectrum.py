#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from oda_api.json import CustomJSONEncoder

get_ipython().run_line_magic("matplotlib", "inline")   # noqa: F821

from oda_api.api import ProgressReporter

pr = ProgressReporter()

import os

workdir = os.getcwd()

class AnalysisError(RuntimeError): ...

src_name = "TXS 0614-083"  # http://odahub.io/ontology#AstrophysicalObject
RA = 94.2635  # http://odahub.io/ontology#PointOfInterestRA
DEC = -8.374  # http://odahub.io/ontology#PointOfInterestDEC
# src_name='Mrk 421' #http://odahub.io/ontology#AstrophysicalObject
# RA = 166.1138083333333  # http://odahub.io/ontology#PointOfInterestRA
# DEC = 38.20883277777778 # http://odahub.io/ontology#PointOfInterestDEC

T1 = "2017-03-06T13:26:48.000"  # http://odahub.io/ontology#StartTime
T2 = "2017-06-06T15:32:27.000"  # http://odahub.io/ontology#EndTime
Radius = 0.1  # http://odahub.io/ontology#AngleDegrees ; oda:label "Radius of the source region"

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in ["src_name", "RA", "DEC", "T1", "T2", "Radius"]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

c = SkyCoord(RA, DEC, unit="degree")

hdul = fits.open(workdir + "/eRASS1_Main.v1.1.fits")
erosita = hdul[1].data
ras_eros = erosita["RA"]
decs_eros = erosita["DEC"]
flux = []
flux_err_down = []
flux_err_up = []
for i in range(5):
    flux.append(erosita["ML_FLUX_P" + str(i + 1)])
    flux_err_down.append(erosita["ML_FLUX_LOWERR_P" + str(i + 1)])
    flux_err_up.append(erosita["ML_FLUX_UPERR_P" + str(i + 1)])
coords_eros = SkyCoord(ras_eros, decs_eros, unit="degree")
seps_eros = c.separation(coords_eros).deg
print(min(seps_eros))
flux = np.array(flux)
flux_err_down = np.array(flux_err_down)
flux_err_up = np.array(flux_err_up)

m = seps_eros < Radius
if sum(m) == 0:
    raise AnalysisError("No data found")

m = seps_eros < Radius
F_eros = np.sum(flux[:, m], axis=1)
F_eros_err_down = np.sqrt(np.sum(flux_err_down[:, m] ** 2, axis=1))
F_eros_err_up = np.sqrt(np.sum(flux_err_up[:, m] ** 2, axis=1))

Emins_eros = np.array([0.2, 0.5, 1.0, 2.0, 5.0])
Emaxs_eros = np.array([0.5, 1.0, 2.0, 5.0, 8.0])
Emids_eros = np.sqrt(Emins_eros * Emaxs_eros)

F_eros = F_eros * Emids_eros / (Emaxs_eros - Emins_eros)
F_eros_dE_err_down = F_eros_err_down * Emids_eros / (Emaxs_eros - Emins_eros)
F_eros_err_up = F_eros_err_up * Emids_eros / (Emaxs_eros - Emins_eros)

UL_eros = F_eros <= F_eros_err_down
F_eros = F_eros * (1 - UL_eros) + UL_eros * 1.5 * F_eros_err_up
F_eros_err_up = F_eros_err_up * (1 - UL_eros) + UL_eros * F_eros_err_up
F_eros_err_down = F_eros_err_down * (1 - UL_eros) + UL_eros * F_eros_err_up
print(F_eros)
# 0.2–0.5 0.5–1.0 1.0–2.0 2.0–5.0 5.0–8.0 4.0–10.0 5.1–6.1 6.2–7.1 7.2–8.2
plt.errorbar(
    Emids_eros,
    F_eros,
    yerr=[F_eros_err_down, F_eros_err_up],
    xerr=[Emids_eros - Emins_eros, Emaxs_eros - Emids_eros],
    uplims=UL_eros,
    linestyle="none",
    marker="o",
    label="eROSITA",
)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("$E$, keV")
plt.ylabel("$E^2 dN/dE$, erg/cm$^2$s")
plt.legend(loc="lower right")
plt.savefig("Spectrum.png", format="png", bbox_inches="tight")

from astropy.table import Table
from oda_api.data_products import ODAAstropyTable, PictureProduct

bin_image = PictureProduct.from_file("Spectrum.png")

data = [
    Emids_eros,
    Emins_eros,
    Emaxs_eros,
    F_eros,
    F_eros_err_down,
    F_eros_err_up,
]
names = (
    "Emean[keV]",
    "Emin[keV]",
    "Emax[keV]",
    "E2dN_dE[erg/cm2s]",
    "E2dN_dE_error_down[erg/cm2s]",
    "E2dN_dE_error_up[erg/cm2s]",
)
spec = ODAAstropyTable(Table(data, names=names))

picture_png = bin_image  # http://odahub.io/ontology#ODAPictureProduct
spectrum_astropy_table = spec  # http://odahub.io/ontology#ODAAstropyTable

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    ("out_Spectrum_picture_png", "picture_png_galaxy.output", picture_png)
)
_oda_outs.append(
    (
        "out_Spectrum_spectrum_astropy_table",
        "spectrum_astropy_table_galaxy.output",
        spectrum_astropy_table,
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
