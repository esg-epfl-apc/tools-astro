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

# src_name='TXS 0614-083' #http://odahub.io/ontology#AstrophysicalObject
# RA = 94.2635  # http://odahub.io/ontology#PointOfInterestRA
# DEC = -8.374 # http://odahub.io/ontology#PointOfInterestDEC

src_name = "Mrk 421"  # http://odahub.io/ontology#AstrophysicalObject
RA = 166.1138083333333  # http://odahub.io/ontology#PointOfInterestRA
DEC = 38.20883277777778  # http://odahub.io/ontology#PointOfInterestDEC

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

Emin = np.array([0.1e3])
Emax = np.array([2.4e3])
Emid = np.sqrt(Emin * Emax)

hdul = fits.open("IX_29_rass_fsc.dat.gz.fits")
fsc = hdul[1].data
ras = fsc["RAdeg"]
decs = fsc["DEdeg"]
rate = fsc["Count"]
rate_err = fsc["e_Count"]
f1 = 10.4e-12 * rate / (Emax - Emin) * Emid
f1_err = 10.4e-12 * rate_err / (Emax - Emin) * Emid
f2 = 5.6e-12 * rate / (Emax - Emin) * Emid
f2_err = 5.6e-12 * rate_err / (Emax - Emin) * Emid
f3 = 3.8e-12 * rate / (Emax - Emin) * Emid
f3_err = 3.8e-12 * rate_err / (Emax - Emin) * Emid
coords = SkyCoord(ras, decs, unit="degree")
seps = c.separation(coords).deg

from astropy.io import fits

hdul = fits.open("cat2rxs.fits.gz")
rosat_cat = hdul[1].data
ras2 = rosat_cat["RA_DEG"]
decs2 = rosat_cat["DEC_DEG"]
coords2 = SkyCoord(ras2, decs2, unit="degree")
seps2 = coords2.separation(c).deg

NH_gal = rosat_cat["NH_gal"]
NORM_p = rosat_cat["NORM_p"]
NORM_err_p = rosat_cat["NORM_err_p"]
GAMMA_p = rosat_cat["GAMMA_p"]
GAMMA_err_p = rosat_cat["GAMMA_err_p"]

mm = seps2 < 1
NH = NH_gal[m]
plt.hist(NH)
np.median(NH)

m = seps < Radius
if sum(m) > 0:
    Emin = np.array([0.1e3])
    Emax = np.array([2.4e3])
    Emid = np.sqrt(Emin * Emax)
    plt.errorbar(Emid, f2[m], xerr=[Emid - Emin, Emax - Emid], yerr=f2_err[m])
    E_rosat = Emid
    Emax_rosat = Emax
    Emin_rosat = Emin
    F_rosat = f2[m]
    F_err_rosat = f2_err[m]

m = seps2 < Radius
if sum(m) > 0:
    N = NORM_p[m][0]
    Nerr = NORM_err_p[m][0]
    G = GAMMA_p[m][0]
    Gerr = GAMMA_err_p[m][0]

    x = np.logspace(-1, 1, 10)
    y = N * x ** (G) * x**2 * 1.6e-9
    y1 = (N + Nerr) * x ** (G + Gerr) * x**2 * 1.6e-9
    ymax = np.maximum(y, y1)
    ymin = np.minimum(y, y1)
    y2 = (N + Nerr) * x ** (G - Gerr) * x**2 * 1.6e-9
    ymax = np.maximum(ymax, y2)
    ymin = np.minimum(ymin, y2)
    y2 = (N - Nerr) * x ** (G + Gerr) * x**2 * 1.6e-9
    ymax = np.maximum(ymax, y2)
    ymin = np.minimum(ymin, y2)
    y2 = (N - Nerr) * x ** (G - Gerr) * x**2 * 1.6e-9
    ymax = np.maximum(ymax, y2)
    ymin = np.minimum(ymin, y2)
plt.plot(x, y, color="black")
plt.fill_between(x, ymin, ymax, color="black", alpha=0.25, linewidth=0)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("$E$, keV")
plt.ylabel("$E^2 dN/dE$, erg/cm$^2$s")
plt.savefig("Spectrum.png", format="png", bbox_inches="tight")

from astropy.table import Table
from oda_api.data_products import ODAAstropyTable, PictureProduct

bin_image = PictureProduct.from_file("Spectrum.png")

data = [x, y, ymin, ymax]
names = (
    "Emean[keV]",
    "E2dN_dE[erg/cm2s]",
    "E2dN_dE_min[erg/cm2s]",
    "E2dN_dE_max[erg/cm2s]",
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
