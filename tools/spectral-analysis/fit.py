#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python

# This script is generated with nb2galaxy

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

spectrum = "spectrum.fits.gz"  # http://odahub.io/ontology#POSIXPath
rmf = "rmf.fits.gz"  # http://odahub.io/ontology#POSIXPath
arf = "arf.fits.gz"  # http://odahub.io/ontology#POSIXPath
background = ""
model = "powerlaw"  # oda:allowed_value 'powerlaw' .
e_min = 25
e_max = 100

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "C_data_product_" in inp_dic.keys():
    inp_pdic = inp_dic["C_data_product_"]
else:
    inp_pdic = inp_dic
spectrum = str(inp_pdic["spectrum"])
rmf = str(inp_pdic["rmf"])
arf = str(inp_pdic["arf"])
background = str(inp_pdic["background"])
model = str(inp_pdic["model"])
e_min = int(inp_pdic["e_min"])
e_max = int(inp_pdic["e_max"])

import numpy as np
from astropy.io import fits

d = fits.open(spectrum)[1].data
rate = d["RATE"]
rate_err = d["STAT_ERR"]

ebounds_ext = fits.open(rmf)["EBOUNDS"].data
c_e1 = ebounds_ext["E_MIN"]
c_e2 = ebounds_ext["E_MAX"]

matrix_ext = fits.open(rmf)["SPECRESP MATRIX"].data
matrix = np.vstack(matrix_ext["MATRIX"])
e1 = matrix_ext["ENERG_LO"]
e2 = matrix_ext["ENERG_HI"]

from matplotlib import pylab as plt

c_de = c_e2 - c_e1

def model(e1, e2, N, slope):
    return N * (e1 / 25.0) ** slope

N = 2e-4
slope = -3.5

def convolve(p, return_rate=False):
    N, slope = p

    model_rate = np.dot(model(e1, e2, N, slope), matrix)

    m = c_e1 > e_min
    m &= c_e2 < e_max

    chi2 = (((model_rate - rate) / rate_err)[m] ** 2).sum()

    # print(p, chi2, chi2 / np.sum(m))

    if return_rate:
        return model_rate
    else:
        return chi2

from scipy.optimize import minimize

f = minimize(convolve, [N, slope])

model_rate = convolve(f.x, True)

plt.errorbar(c_e1, rate / c_de, rate_err)
plt.plot(c_e1, model_rate / c_de)

plt.axvspan(0, e_min, alpha=0.2, color="r")
plt.axvspan(e_max, 1000, alpha=0.2, color="r")

plt.loglog()

plt.savefig("fit.png")

convolve(f.x, True)

from oda_api.data_products import PictureProduct

outlcpng = PictureProduct.from_file("fit.png")

fitpng = outlcpng  # oda:ODAPictureProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(("out_fit_fitpng", "fitpng_galaxy.output", fitpng))

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
