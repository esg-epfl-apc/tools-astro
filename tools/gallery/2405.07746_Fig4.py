#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
from oda_api.api import DispatcherAPI
from oda_api.json import CustomJSONEncoder
from oda_api.token import discover_token

disp = DispatcherAPI(
    url="https://www.astro.unige.ch/mmoda//dispatch-data", instrument="mock"
)

token = discover_token()

# Parameters of the phase transition
# Temperature in GeV
T_star = 0.178  # http://odahub.io/ontology#Energy_GeV ; oda:label "Phase transition temperature"

# Numbers of relativistic degrees of freedom
g_star = 20  # http://odahub.io/ontology#Integer ; oda:label "Number of effective degrees of freedom"

# ratio of the energy density deposited in the bubbles to the radiation energy density
alpha = 1.0  # http://odahub.io/ontology#Float ; oda:label "alpha: latent heat in units of radiation density"

# beta/H : rate of the phase transition compared to Hubble rate
beta_H = 3.3  # http://odahub.io/ontology#Float ; oda:label "beta: phase transition rate in units of Hubble rate"

# fraction of turbulent energy that goes to gw N.B. arXiv:1004.4187 claims that epsilon_turb=0.05, but checks below show that it is rather 0.01
epsilon_turb = 1  # http://odahub.io/ontology#Float ; oda:label "fraction of energy release in turbulence"

# terminal velocity of bubbles
v_w = 0.999  # http://odahub.io/ontology#Float ; oda:label "ternminal velocity of bubles"
h = 0.7  # http://odahub.io/ontology#Float ; oda:label "Hubble parameter"

Parameterisation = "RoperPol2023"  # http://odahub.io/ontology#String ; ; oda:label "model parameterisation" ; oda:allowed_value "Lewicki2022","RoperPol2023"

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in [
    "T_star",
    "g_star",
    "alpha",
    "beta_H",
    "epsilon_turb",
    "v_w",
    "h",
    "Parameterisation",
]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

workdir = os.getcwd()
repo_basedir = os.environ.get("BASEDIR", os.getcwd())
data_dir = repo_basedir + "/data_2405.07746"
get_ipython().system("ls {data_dir}")   # noqa: F821

d = np.genfromtxt(data_dir + "/NANOGrav23.csv")
gammas_nano = d[:, 0]
As_nano = 10 ** d[:, 1]
ps_nano = 5 - gammas_nano

d = np.genfromtxt(data_dir + "/EPTA.csv")
gammas_epta = d[:, 0]
As_epta = 10 ** d[:, 1]
ps_epta = 5 - gammas_epta

d = np.genfromtxt(data_dir + "/PPTA.csv")
gammas_ppta = d[:, 0]
As_ppta = 10 ** d[:, 1]
ps_ppta = 5 - gammas_ppta

H0 = 70 * (u.km / u.s) / u.Mpc

par_dict = {
    "T_star": T_star,
    "alpha": alpha,
    "beta_H": beta_H,
    "epsilon_turb": epsilon_turb,
    "g_star": g_star,
    "h": 0.7,
    "instrument": "sgwb",
    "product": "Model_spectrum",
    "product_type": "Real",
    "token": token,
    "v_w": v_w,
    "Parameterisation": Parameterisation,
}

data_collection = disp.get_product(**par_dict)

data_collection.astropy_table_0.table

ff = data_collection.astropy_table_0.table["f[Hz]"]
if Parameterisation == "Lewicki2022":
    GW = data_collection.astropy_table_0.table["Omega_gw"]
    plt.plot(ff, GW, linewidth=4, color="black", alpha=0.5, label="total")
else:
    GW_s = data_collection.astropy_table_0.table["Omega_sound_waves"]
    GW_t = data_collection.astropy_table_0.table["Omega_turbulence"]
    GW = GW_s + GW_t
    plt.plot(ff, GW_t, color="blue", linestyle="dashed", label="turbulence")
    plt.plot(ff, GW_s, color="red", linestyle="dotted", label="sound waves")
    plt.plot(ff, GW, linewidth=4, color="black", alpha=0.5, label="total")

fref = (1 / u.yr).cgs.value
lgfmin = np.log10(fref / 10.0)
lgfmax = np.log10(fref / 2.0)
fff = np.logspace(lgfmin, lgfmax, 10) * u.Hz
min_nano = np.ones(len(fff))
max_nano = np.zeros(len(fff))
for i in range(len(As_nano)):
    spec = (
        2
        * pi**2
        / 3
        / H0**2
        * fff**2
        * As_nano[i] ** 2
        * (fff / fref) ** (3 - gammas_nano[i])
    ).cgs.value
    min_nano = np.minimum(spec, min_nano)
    max_nano = np.maximum(spec, max_nano)
    # plt.plot(ff,spec)
plt.fill_between(
    fff.value, min_nano, max_nano, color="red", alpha=0.5, label="NANOGrav"
)
min_epta = np.ones(len(fff))
max_epta = np.zeros(len(fff))
for i in range(len(As_epta)):
    spec = (
        2
        * pi**2
        / 3
        / H0**2
        * fff**2
        * As_epta[i] ** 2
        * (fff / fref) ** (3 - gammas_epta[i])
    ).cgs.value
    min_epta = np.minimum(spec, min_epta)
    max_epta = np.maximum(spec, max_epta)
    # plt.plot(ff,spec)
plt.fill_between(
    fff.value, min_epta, max_epta, color="blue", alpha=0.5, label="EPTA"
)
min_ppta = np.ones(len(fff))
max_ppta = np.zeros(len(fff))
for i in range(len(As_ppta)):
    spec = (
        2
        * pi**2
        / 3
        / H0**2
        * fff**2
        * As_ppta[i] ** 2
        * (fff / fref) ** (3 - gammas_ppta[i])
    ).cgs.value
    min_ppta = np.minimum(spec, min_ppta)
    max_ppta = np.maximum(spec, max_ppta)
    # plt.plot(ff,spec)
plt.fill_between(
    fff.value, min_ppta, max_ppta, color="green", alpha=0.5, label="PPTA"
)

maxGW = max(GW)
ind = np.argmax(GW)
fmax = ff[ind]
plt.xscale("log")
plt.yscale("log")
plt.ylim(maxGW / 1e5, maxGW * 10)
plt.xlim(fmax / 1e3, fmax * 1e3)
plt.grid()
plt.xscale("log")
plt.yscale("log")
plt.legend(loc="upper left")
plt.xlabel("$f$, Hz")
plt.ylabel("$\Omega_{gw}(f)$")
plt.savefig("Spectrum.png", format="png", bbox_inches="tight")

from oda_api.data_products import PictureProduct

bin_image = PictureProduct.from_file("Spectrum.png")

png = bin_image  # http://odahub.io/ontology#ODAPictureProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(("out_2405.07746_Fig4_png", "png_galaxy.output", png))

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
