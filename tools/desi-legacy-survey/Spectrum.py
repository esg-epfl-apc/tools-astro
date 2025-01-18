#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import astropy.units as u
import numpy as np

# Conventional astronomical tools, also to be traced by Renku plugin, there is domain-specific ontology built in
from astropy.coordinates import Angle  # Angles
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.table import Table
from desi import DESILegacySurvey
from numpy import pi, sqrt
from oda_api.json import CustomJSONEncoder

nano_maggies_to_microJy = 3.631  # constant of conversion
nano_maggies_to_Jy = 3.631e-6  # constant of conversion

import matplotlib.pyplot as plt

src_name = "Mrk 421"  # http://odahub.io/ontology#AstrophysicalObject
RA = 166.113808  # http://odahub.io/ontology#PointOfInterestRA
DEC = 38.208833  # http://odahub.io/ontology#PointOfInterestDEC
Radius = 0.1  # http://odahub.io/ontology#AngleMinutes
data_release = (
    10  # http://odahub.io/ontology#Integer ; oda:label "Data Release"
)

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in ["src_name", "RA", "DEC", "Radius", "data_release"]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

def clean_flux(flux, flux_ivar):
    new_flux = flux * nano_maggies_to_microJy
    new_flux_ivar = nano_maggies_to_microJy / np.sqrt(flux_ivar)
    return new_flux, new_flux_ivar

ra_s = RA
dec_s = DEC
dr = data_release
source = SkyCoord(ra_s, dec_s, unit="degree")
Radius = Angle(Radius * u.arcmin)

case_ = 0
if source.galactic.b > 0:
    if dec_s >= 32:
        print("MzLS")
        case_ = 1
    else:
        print("DECALS")
else:
    print("DECALS")

error_message = f"No data found, i.e. (RA, Dec) = ({ra_s}, {dec_s}) is outside the covered region by Data Release {dr}."
try:
    query = DESILegacySurvey.query_region(
        coordinates=source, radius=Radius, data_release=dr
    )
except:
    raise RuntimeError(error_message)

if len(query) == 0:
    raise RuntimeError(error_message)

tap_result = query

wavelength = np.array(
    [4770, 6231, 9134, 3.368e4, 4.618e4, 12.082e4, 22.194e4]
)  # in Angstroem
wavelength = wavelength * 1e-8  # in cm
frequency = 3.0e10 / wavelength  # in Hz
energy = 2 * pi * 6.67e-16 * frequency  # in eV
factor = 1e-23 * frequency  # conversion of Jy to erg/cm2s

filters = {
    "g": 4770,
    "r": 6231,
    "z": 9134,
    "w1": 3.368e4,
    "w2": 4.618e4,
    "w3": 12.082e4,
    "w4": 22.194e4,
}
nsources = len(query["ra"])
ras = np.array(query["ra"])
decs = np.array(query["dec"])

coords = SkyCoord(ra=ras, dec=decs, unit="deg", frame="icrs")
sep = source.separation(coords).arcmin
t = query["type"]
for i in range(nsources):
    if t[i] == "DUP":
        sep[i] = 100.0
index = np.argmin(sep)

label = ["g", "r", "z", "w1", "w2", "w3", "w4"]
flux = np.zeros(len(label))
flux_err = np.zeros(len(label))
for ind in range(nsources):
    if t[ind] != "DUP":
        print(ind)
        y1 = np.array(
            [
                query[ind]["flux_" + i] / query["mw_transmission_" + i][ind]
                for i in label
            ]
        )  # in nanomaggies
        y1_err = np.array(
            [
                1.0
                / sqrt(query[ind]["flux_ivar_" + i])
                / query["mw_transmission_" + i][ind]
                for i in label
            ]
        )  # in nanomaggies
        y = y1 * nano_maggies_to_Jy * factor  # in Jy
        y_err = y1_err * nano_maggies_to_Jy * factor  # in Jy
        flux += y
        flux_err += y_err**2
        plt.errorbar(
            energy,
            y,
            y_err,
            marker="o",
            label="source " + str(query[ind]["objid"]),
            alpha=0.3,
        )
flux_err = sqrt(flux_err)
plt.errorbar(
    energy,
    flux,
    flux_err,
    marker="o",
    color="black",
    linewidth=2,
    label="total",
)
plt.yscale("log")
plt.xscale("log")
plt.legend(loc="lower right")
plt.xlabel("$E$, [eV]", fontsize=16)
plt.ylabel("$E^2dN/dE$, [erg/(cm$^2$s)]", fontsize=16)
plt.savefig("Spectrum.png")

from oda_api.data_products import ODAAstropyTable

data = [energy, flux, flux_err]
names = ("Energy[eV]", "Flux[erg/cm2s]", "Flux_err[erg/cm2s]")
spec = ODAAstropyTable(Table(data, names=names))

from oda_api.data_products import PictureProduct

bin_image = PictureProduct.from_file("Spectrum.png")

spectrum_table = spec  # http://odahub.io/ontology#ODAAstropyTable
spectrum_png = bin_image  # http://odahub.io/ontology#ODAPictureProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_Spectrum_spectrum_table",
        "spectrum_table_galaxy.output",
        spectrum_table,
    )
)
_oda_outs.append(
    ("out_Spectrum_spectrum_png", "spectrum_png_galaxy.output", spectrum_png)
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
