#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astroquery.ipac.irsa import Irsa
from oda_api.json import CustomJSONEncoder

src_name = "Mrk 421"  # http://odahub.io/ontology#AstrophysicalObject
RA = 166.1138083333333  # http://odahub.io/ontology#PointOfInterestRA
DEC = 38.20883277777778  # http://odahub.io/ontology#PointOfInterestDEC
T1 = "2017-03-06T13:26:48.000"  # http://odahub.io/ontology#StartTime
T2 = "2017-06-06T15:32:27.000"  # http://odahub.io/ontology#EndTime
Source_region_radius = 1.0  # http://odahub.io/ontology#AngleArcsec ; oda:label "Radius of the source region"

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in ["src_name", "RA", "DEC", "T1", "T2", "Source_region_radius"]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

c = SkyCoord(RA, DEC, unit="degree")
RA, DEC

table = Irsa.query_region(
    coordinates=c, catalog="allwise_p3as_psd", radius="0d0m10s"
)

m_wise = np.array(
    [
        table["w1mpro"][0],
        table["w2mpro"][0],
        table["w2mpro"][0],
        table["w3mpro"][0],
    ]
)
m_wise

nu_wise = np.array([8.8560e13, 6.4451e13, 2.6753e13, 1.3456e13])

m_wise_err = np.array(
    [
        table["w1sigmpro"][0],
        table["w2sigmpro"][0],
        table["w2sigmpro"][0],
        table["w3sigmpro"][0],
    ]
)
m_wise_err

F0_wise = np.array([309.54, 171.787, 31.674, 8.363]) * 1e-23
F_wise = F0_wise * 10 ** (-m_wise / 2.5) * nu_wise
F_wise_err = F0_wise * 10 ** (-(m_wise - m_wise_err) / 2.5) * nu_wise - F_wise
F_wise, F_wise_err

import astropy.units as u
from astropy.constants import h

h.to(u.eV * u.s).value
E_wise = h.to(u.eV * u.s).value * nu_wise

plt.errorbar(E_wise, F_wise, F_wise_err)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("$E$, eV")
plt.ylabel("$E^2 dN/dE$, erg/cm$^2$s")
plt.savefig("Spectrum.png", format="png", bbox_inches="tight")

from astropy.table import Table
from oda_api.data_products import ODAAstropyTable, PictureProduct

bin_image = PictureProduct.from_file("Spectrum.png")

data = [E_wise, F_wise, F_wise_err]
names = ("E[eV]", "E2dN_dE[erg/cm2s]", "E2dN_dE_err[erg/cm2s]")
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
