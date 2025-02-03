#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.constants import h
from astroquery.skyview import SkyView
from numpy import sqrt
from oda_api.data_products import ODAAstropyTable, PictureProduct
from oda_api.json import CustomJSONEncoder
from pyvo import registry  # version >=1.4.1

class AnalysisError(RuntimeError): ...

h_p = h.to(u.eV * u.s).value

src_name = "1ES 0229+200"  # http://odahub.io/ontology#AstrophysicalObject
RA = 38.202562  # http://odahub.io/ontology#PointOfInterestRA
DEC = 20.288191  # http://odahub.io/ontology#PointOfInterestDEC
T1 = "2000-10-09T13:16:00.0"  # http://odahub.io/ontology#StartTime
T2 = "2022-10-10T13:16:00.0"  # http://odahub.io/ontology#EndTime

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in ["src_name", "RA", "DEC", "T1", "T2"]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

Radius = 0.1  # http://odahub.io/ontology#AngleDegrees
Imsize = 0.2  # http://odahub.io/ontology#AngleDegrees
pixsize = 0.01
pixels = int(2 * Imsize / pixsize) + 1
source_pix = int(Imsize / pixsize)
from astropy.coordinates import SkyCoord

coords_s = SkyCoord(RA, DEC, unit="degree")
pos = str(RA) + ", " + str(DEC)
pos

from astropy.io import fits

workdir = os.getcwd()

hdul = fits.open(workdir + "/GLEAM_EGC_v2.fits.gz")
table = hdul[1].data

ras = table["RAJ2000"]
decs = table["DEJ2000"]
# flux=[table['int_flux_084']
columns = table.columns
coords = SkyCoord(ras, decs, unit="degree")
seps = coords.separation(coords_s).deg

freq_list = [
    "GLEAM 72-103 MHz",
    "GLEAM 103-134 MHz",
    "GLEAM 139-170 MHz",
    "GLEAM 170-231 MHz",
]
fmin = np.array([73.0, 103.0, 139.0, 170.0]) * 1e6
fmax = np.array([103.0, 134.0, 170.0, 231.0]) * 1e6
fmid = sqrt(fmin * fmax)
Emid = h_p * fmid
Emax = h_p * fmax
Emin = h_p * fmin

try:
    images = SkyView.get_images(
        position=pos, survey=freq_list, pixels=pixels, radius=Imsize * u.deg
    )
except:
    raise AnalysisError("No data found")
    message = "No data found!"

im = images[-1][0].data
plt.imshow(im)

f = []
ferr = []
for i in range(len(images)):
    im = images[i][0].data
    f.append(im[source_pix, source_pix] * fmid[i] * 1e-23)
    im_flat = im.flatten()
    f1 = max(im_flat)
    f0 = min(im_flat)
    fbins = np.linspace(f0, f1, 1001)
    fvalues = (fbins[:-1] + fbins[1:]) / 2.0
    h = np.histogram(im_flat, bins=fbins)
    profile = h[0]
    profile_cum = np.cumsum(profile)
    profile_cum = profile_cum / profile_cum[-1]
    plt.plot(fvalues, profile_cum)
    m = profile_cum < 0.5 - 0.34
    f1 = fvalues[m][-1]
    m = profile_cum < 0.5 + 0.34
    f2 = fvalues[m][-1]
    ferr.append((f2 - f1) / 2.0 * fmid[i] * 1e-23)
    plt.axvline(np.median(im_flat), color="red")
    print(f[-1], ferr[-1])

m = seps < Radius

if len(table[m]) > 0.0:
    nu = []
    flux = []
    flux_err = []
    for col in columns:
        if (
            (col.name[:8] == "int_flux")
            and (col.name[9] != "w")
            and (col.name[9] != "f")
        ):
            nu.append(int(col.name[9:12]) * 1e6)
            flux.append(sum(table[col.name][m]) * nu[-1] * 1e-23)
            flux_err.append(
                sqrt(sum((table["err_" + col.name][m]) ** 2)) * nu[-1] * 1e-23
            )
    E = h_p * np.array(nu)

if sum(m) > 0:
    plt.errorbar(
        E, flux, flux_err, label="sum of catalog source fluxes within aperture"
    )
plt.errorbar(
    Emid,
    f,
    yerr=ferr,
    xerr=[Emid - Emin, Emax - Emid],
    linestyle="none",
    color="black",
    linewidth=2,
    marker="o",
    label="point source flux at requested position (from images)",
)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("$E$, eV")
plt.ylabel("$E F_E$, erg/cm$^2$s")
plt.legend(loc="lower right")
plt.savefig("Spectrum.png", format="png", bbox_inches="tight")

bin_image = PictureProduct.from_file("Spectrum.png")
from astropy.table import Table

names = ("E[eV]", "frequency[Hz]", "Flux[erg/cm2s]", "Flux_error[erg/cm2s]")
if sum(m) > 0:
    data = [E, nu, flux, flux_err]
    spec = ODAAstropyTable(Table(data, names=names))
else:
    data = [[0], [0], [0], [0]]
    spec = ODAAstropyTable(Table(data, names=names))

names = (
    "E[eV]",
    "Emin[eV]",
    "Emax[eV]",
    "f[Hz]",
    "fmin[Hz]",
    "fmax[Hz]",
    "Flux[erg/cm2s]",
    "Flux_error[erg/cm2s]",
)
data = [Emid, Emin, Emax, fmid, fmin, fmax, f, ferr]
spec_image = ODAAstropyTable(Table(data, names=names))

picture_png = bin_image  # http://odahub.io/ontology#ODAPictureProduct
spectrum_catalog_table = spec  # http://odahub.io/ontology#ODAAstropyTable
spectrum_image_table = spec_image  # http://odahub.io/ontology#ODAAstropyTable

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    ("out_Spectrum_picture_png", "picture_png_galaxy.output", picture_png)
)
_oda_outs.append(
    (
        "out_Spectrum_spectrum_catalog_table",
        "spectrum_catalog_table_galaxy.output",
        spectrum_catalog_table,
    )
)
_oda_outs.append(
    (
        "out_Spectrum_spectrum_image_table",
        "spectrum_image_table_galaxy.output",
        spectrum_image_table,
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
