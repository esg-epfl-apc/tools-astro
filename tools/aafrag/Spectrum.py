#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import copy
import json
import os
import re
import shutil

import matplotlib.pyplot as plt
import numpy as np

# import AAFragpy routine for differential cross-section
from aafragpy import (
    get_cross_section,
    get_cross_section_Kamae2006,
    get_spectrum,
)
from numpy import sqrt
from oda_api.json import CustomJSONEncoder

dN_dp = "(p/1e4)**(-1.99)*exp(-p/1e4)"  # http://odahub.io/ontology#String ; oda:label "Source spectrum dN/dp"
sec = "nu_all"  # http://odahub.io/ontology#String ; oda:label "secondary particle channel" ; oda:allowed_value "gam", "el", "pos", "nu_e", "anu_e", "nu_mu", "anu_mu", "nu_all"

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in ["dN_dp", "sec"]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

def parse_spectrum(input_str):
    dnde_str = copy.copy(input_str)

    numbers = re.findall(r"-?\d+\.?\d*(?:[ep][+-]?\d+)?", input_str)
    funcs = re.findall(r"\b\w+(?=\()", input_str)

    for num in numbers:
        input_str = input_str.replace(num, "")

    for func in funcs:
        input_str = input_str.replace(func, "")

    service_symbols = "p()[]*/.,+- "

    for symbol in service_symbols:
        input_str = input_str.replace(symbol, "")

    if input_str:
        raise ValueError("forbidden statements")

    return eval(f"lambda p: {dnde_str}")

Assumed = parse_spectrum(dN_dp)

# $$
# E^2=m^2+p^2;   EdE=pdp
# $$
# $$
# \frac{dN}{dE}=\frac{dN}{dp}\frac{dp}{dE}=\frac{dN}{dp}\frac{E}{p}
# $$
#

m_p = 0.938  # proton rest energy in GeV
E_thr = 4  # "threshold" energy in GeV

E_p = np.logspace(
    0, 10, 100
)  # energies of primary particles in GeV, 10 bins per decade
p_p = sqrt(E_p**2 - m_p**2)  # momenta of primary particles
E_k = E_p - m_p  # kinetic energy of primary particles
dNp_dEp = Assumed(p_p) * E_p / p_p  # differential energy spectrum

plt.plot(E_k, dNp_dEp * E_p**2)
plt.xscale("log")
plt.yscale("log")
plt.ylim(1e0, 1e10)

E_s = np.logspace(-3, 10, 130)  # energy binning for secondary particles

cs_matrix = get_cross_section(
    sec, "p-p", E_primaries=E_p[E_p > E_thr], E_secondaries=E_s
)
cs_matrix_Kamae2006 = get_cross_section_Kamae2006(
    sec, E_primaries=E_p[E_p <= E_thr], E_secondaries=E_s
)
# cs_matrix_Kafexhiu2014=get_cross_section_Kafexhiu2014(E_primaries=E_p[E_p<=E_thr], E_secondaries=E_s)

# spectrum from AAfrag
spec = get_spectrum(
    E_p[E_p > E_thr],
    E_s,
    cs_matrix=cs_matrix[0],
    prim_spectrum=dNp_dEp[E_p > E_thr],
)
# spectrum from Kamae
spec_Kamae2006 = get_spectrum(
    E_p[E_p <= E_thr],
    E_s,
    cs_matrix_Kamae2006[0],
    prim_spectrum=dNp_dEp[E_p <= E_thr],
)
# spec_Kafexhiu2014=get_spectrum(E_p[E_p<=E_thr],E_s,cs_matrix_Kafexhiu2014[0],prim_spectrum=dNp_dEp[E_p<=E_thr])

# plt.figure(figsize=(10,7))

plt.loglog(E_s, (spec + spec_Kamae2006) * E_s**2, "b-", label=sec)
plt.xlabel("$E$, GeV")
plt.ylabel("$E^2 \cdot dN/dE$")
plt.tick_params(axis="both", which="major")
plt.grid()

ymax = np.amax(spec * E_s**2)
m = spec * E_s**2 > 1e-5 * ymax
xmin = E_s[m][0]
xmax = E_s[m][-1]
plt.ylim(ymax / 1e5, ymax * 2)
plt.xlim(xmin, xmax)
plt.savefig("Spectrum.png")

from astropy.table import Table
from oda_api.data_products import ODAAstropyTable, PictureProduct

data = [E_s, (spec + spec_Kamae2006) * E_s**2]
names = ("E[GeV]", "E2dN_dE")
spec = ODAAstropyTable(Table(data, names=names))

bin_image = PictureProduct.from_file("Spectrum.png")

spectrum_png = bin_image  # http://odahub.io/ontology#ODAPictureProduct
spectrum_table = spec  # http://odahub.io/ontology#ODAAstropyTable

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    ("out_Spectrum_spectrum_png", "spectrum_png_galaxy.output", spectrum_png)
)
_oda_outs.append(
    (
        "out_Spectrum_spectrum_table",
        "spectrum_table_galaxy.output",
        spectrum_table,
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
