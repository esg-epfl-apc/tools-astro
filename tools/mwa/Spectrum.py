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
from oda_api.data_products import ODAAstropyTable, PictureProduct
from oda_api.json import CustomJSONEncoder
from pyvo import registry  # version >=1.4.1

src_name = "1ES 0229+200"  # http://odahub.io/ontology#AstrophysicalObject
RA = 38.202562  # http://odahub.io/ontology#PointOfInterestRA
DEC = 20.288191  # http://odahub.io/ontology#PointOfInterestDEC
T1 = "2000-10-09T13:16:00.0"  # http://odahub.io/ontology#StartTime
T2 = "2022-10-10T13:16:00.0"  # http://odahub.io/ontology#EndTime
Radius = 0.05  # http://odahub.io/ontology#AngleDegrees

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in ["src_name", "RA", "DEC", "T1", "T2", "Radius"]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

conesearch_radius = Radius  # in degrees
conesearch_center = (RA, DEC)

# the catalogue name in VizieR
CATALOGUE = "VIII/100"
# each resource in the VO has an identifier, called ivoid. For vizier catalogs,
# the VO ids can be constructed like this:
catalogue_ivoid = f"ivo://CDS.VizieR/{CATALOGUE}"
# the actual query to the registry
voresource = registry.search(ivoid=catalogue_ivoid)[0]

conesearch_records = voresource.get_service("conesearch").search(
    pos=conesearch_center,
    sr=conesearch_radius,
)
conesearch_records

h_p = (h / u.s).to(u.eV).value  # Planck constant in eV*s

conesearch_records.fieldnames
nu = []
F = []
F_err = []
for f in conesearch_records.fieldnames:
    if (f[:4] == "Fint") and (f[4] != "w") and (f[4] != "f"):
        nu.append(int(f[-3:]) * 1e6)  # in Hz
        F.append(conesearch_records[f][0] * nu[-1] * 1e-23)
        F_err.append(conesearch_records["e_" + f][0] * nu[-1] * 1e-23)

E = h_p * np.array(nu)
plt.errorbar(E, F, F_err)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("$E$, eV")
plt.ylabel("$E F_E$, erg/cm$^2$s")
plt.savefig("Spectrum.png", format="png", bbox_inches="tight")

bin_image = PictureProduct.from_file("Spectrum.png")
from astropy.table import Table

data = [E, F, F_err]
names = ("E[eV]", "Flux[erg/cm2s]", "Flux_error[erg/cm2s]")
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
