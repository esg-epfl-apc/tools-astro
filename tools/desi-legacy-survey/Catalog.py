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
from oda_api.json import CustomJSONEncoder

nano_maggies_to_microJy = 3.631  # constant of conversion

src_name = "Mrk 421"  # http://odahub.io/ontology#AstrophysicalObject
RA = 166.113808  # http://odahub.io/ontology#PointOfInterestRA
DEC = 38.208833  # http://odahub.io/ontology#PointOfInterestDEC
Radius = 1  # http://odahub.io/ontology#AngleMinutes
data_release = 9  # http://odahub.io/ontology#Integer

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in ["src_name", "RA", "DEC", "Radius", "data_release"]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

# https://arxiv.org/pdf/2208.08513  Table 2
# def compute_magnitude(flux, mw_transmission):
#     return 22.5 - 2.5 * np.log10(flux / mw_transmission)

# def compute_error_mag(flux, flux_ivar):
#     dflux = 1 / np.sqrt(flux_ivar)
#     # 2.5 / ln(10) * (dflux / flux)
#     return 1.0857 * dflux / flux

def clean_flux(flux, flux_ivar):
    new_flux = np.exp(np.log(flux * nano_maggies_to_microJy))
    new_flux_ivar = np.exp(
        np.log(nano_maggies_to_microJy / np.sqrt(flux_ivar))
    )
    return new_flux, new_flux_ivar

ra_s = RA
dec_s = DEC
dr = data_release
source = SkyCoord(ra_s, dec_s, unit="degree")
Radius = Angle(Radius * u.arcmin)

if source.galactic.b > 0:
    if dec_s >= 32:
        print("MzLS")
    else:
        print("DECALS")
else:
    print("DECALS")

try:
    query = DESILegacySurvey.query_region(
        coordinates=source, radius=Radius, data_release=dr
    )
except:
    raise RuntimeError("No data found. Maybe wrong RA, Dec or Data Release.")

tap_result = query

from oda_api.data_products import ODAAstropyTable

ra = tap_result["ra"]
dec = tap_result["dec"]
t = tap_result["type"]

flux_g, flux_g_err = clean_flux(
    tap_result["flux_g"], tap_result["flux_ivar_g"]
)
flux_r, flux_r_err = clean_flux(
    tap_result["flux_r"], tap_result["flux_ivar_r"]
)
flux_z, flux_z_err = clean_flux(
    tap_result["flux_z"], tap_result["flux_ivar_z"]
)
flux_w1, flux_w1_err = clean_flux(
    tap_result["flux_w1"], tap_result["flux_ivar_w1"]
)
flux_w2, flux_w2_err = clean_flux(
    tap_result["flux_w2"], tap_result["flux_ivar_w2"]
)
flux_w3, flux_w3_err = clean_flux(
    tap_result["flux_w3"], tap_result["flux_ivar_w3"]
)
flux_w4, flux_w4_err = clean_flux(
    tap_result["flux_w4"], tap_result["flux_ivar_w4"]
)

ebv = tap_result["ebv"]

try:
    flux_i, flux_i_err = clean_flux(
        tap_result["flux_i"], tap_result["flux_ivar_i"]
    )
except:
    flux_i = -99 * np.ones_like(flux_g)
    flux_i_err = -99 * np.ones_like(flux_g)
    flux_i, flux_i_err = clean_flux(flux_i, flux_i_err)

data = [
    ra,
    dec,
    t,
    flux_g,
    flux_g_err,
    flux_r,
    flux_r_err,
    flux_z,
    flux_z_err,
    flux_i,
    flux_i_err,
    flux_w1,
    flux_w1_err,
    flux_w2,
    flux_w2_err,
    flux_w3,
    flux_w3_err,
    flux_w4,
    flux_w4_err,
    ebv,
]
names = (
    "RA[deg]",
    "DEC[deg]",
    "Type",
    "flux_g[Jy]",
    "flux_g_err[Jy]",
    "flux_r[Jy]",
    "flux_r_err[Jy]",
    "flux_z[Jy]",
    "flux_z_err[Jy]",
    "flux_i[Jy]",
    "flux_i_err[Jy]",
    "flux_w1[Jy]",
    "flux_w1_err[Jy]",
    "flux_w2[Jy]",
    "flux_w2_err[Jy]",
    "flux_w3[Jy]",
    "flux_w3_err[Jy]",
    "flux_w4[Jy]",
    "flux_w4_err[Jy]",
    "ebv",
)
cat = ODAAstropyTable(Table(data, names=names))

dict_filters = json.dumps(
    {
        "filter": [
            "DECam|DECam.g",
            "DECam|DECam.r",
            "DECam|DECam.z",
            "DECam|DECam.i",
            "WISE|WISE.W1",
            "WISE|WISE.W2",
        ],
        "flux_error": [
            "flux_g_err[Jy]",
            "flux_r_err[Jy]",
            "flux_z_err[Jy]",
            "flux_i_err[Jy]",
            "flux_w1_err[Jy]",
            "flux_w2_err[Jy]",
        ],
        "flux": [
            "flux_g[Jy]",
            "flux_r[Jy]",
            "flux_z[Jy]",
            "flux_i[Jy]",
            "flux_w1[Jy]",
            "flux_w2[Jy]",
        ],
    }
)

catalog_table = cat  # http://odahub.io/ontology#ODAAstropyTable
dictionary_filters = dict_filters  # http://odahub.io/ontology#ODATextProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    ("out_Catalog_catalog_table", "catalog_table_galaxy.output", catalog_table)
)
_oda_outs.append(
    (
        "out_Catalog_dictionary_filters",
        "dictionary_filters_galaxy.output",
        dictionary_filters,
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
