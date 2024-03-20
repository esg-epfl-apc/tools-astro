#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

get_ipython().run_line_magic("load_ext", "autoreload")   # noqa: F821
get_ipython().run_line_magic("autoreload", "2")   # noqa: F821

import os

import numpy as np
import scipy.stats
from matplotlib import pyplot as plt
from oda_api.data_products import ODAAstropyTable, PictureProduct

RA = 40.669622  # http://odahub.io/ontology#PointOfInterestRA
DEC = -0.013294  # http://odahub.io/ontology#PointOfInterestDEC
Slope = 3.0  # http://odahub.io/ontology#Float

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for vn, vv in inp_pdic.items():
    if vn != "_selector":
        globals()[vn] = type(globals()[vn])(vv)

from skyllh.analyses.i3.publicdata_ps.time_integrated_ps import create_analysis
from skyllh.core.config import Config
from skyllh.core.random import RandomStateService
from skyllh.core.source_model import PointLikeSource
from skyllh.datasets.i3.PublicData_10y_ps import create_dataset_collection

cfg = Config()

if os.path.exists("20210126_PS-IC40-IC86_VII.zip") == False:
    get_ipython().system(   # noqa: F821
        "wget http://icecube.wisc.edu/data-releases/20210126_PS-IC40-IC86_VII.zip"
    )
    get_ipython().system("unzip 20210126_PS-IC40-IC86_VII.zip")   # noqa: F821

data_dir = os.getcwd() + "/icecube_10year_ps/"

dsc = create_dataset_collection(cfg=cfg, base_path=data_dir)
dsc.dataset_names

datasets = dsc["IC86_I", "IC86_II-VII"]

rss = RandomStateService(seed=1)

source = PointLikeSource(ra=np.deg2rad(RA), dec=np.deg2rad(DEC))
ana = create_analysis(cfg=cfg, datasets=datasets, source=source)
events_list = [data.exp for data in ana.data_list]
ana.initialize_trial(events_list)
(ts, x, status) = ana.unblind(minimizer_rss=rss)
print(f"TS = {ts:.3f}")
print(f'ns = {x["ns"]:.2f}')
print(f'gamma = {x["gamma"]:.2f}')

rss = RandomStateService(seed=1)
(ts, x, status) = ana.unblind(minimizer_rss=rss)

print(f"TS = {ts:.3f}")
print(f'ns = {x["ns"]:.2f}')
print(f'gamma = {x["gamma"]:.2f}')

TS_profile = []
counts = np.linspace(0, 200, 200)
for n in counts:
    TS_profile.append(2 * ana.llhratio.evaluate([n, Slope])[0])

plt.plot(counts, TS_profile)
plt.axvline(np.argmax(TS_profile), color="black", label="best fit")
tsmax = max(TS_profile)
print(max(TS_profile))
# plt.axhline(tsmax-4.5) #4.5 is for two-parameter adjustment (N_s, sigma), or (N_s, gamma), if sigma=0
# Determine the delta lambda value for the 95% quantile assuming a chi-sqaure
# distribution with 2 degrees of freedom (i.e. assuming Wilks theorem).
chi2_68_quantile = scipy.stats.chi2.ppf(0.68, df=1)
chi2_90_quantile = scipy.stats.chi2.ppf(0.90, df=1)
chi2_95_quantile = scipy.stats.chi2.ppf(0.95, df=1)

# chi2_90_quantile,chi2_95_quantile
mask_68 = TS_profile > tsmax - chi2_68_quantile
cmin = min(counts[mask])
cmax = max(counts[mask])
print("68% confidence interval:", cmin, cmax)
F_min = ana.calculate_fluxmodel_scaling_factor(cmin, [cmin, Slope])
F_max = ana.calculate_fluxmodel_scaling_factor(cmax, [cmax, Slope])
plt.axhline(
    tsmax - chi2_68_quantile, color="black", linestyle="dashed", label="68%"
)
# plt.text(0,tsmax-chi2_68_quantile+1,r'$F_{max}=$'+str(F_max)+' (GeV s cm^2)$^{-1}$ at 1 TeV')
# plt.text(0,tsmax-chi2_68_quantile-1.5,r'$F_{min}=$'+str(F_min)+' (GeV s cm^2)$^{-1}$ at 1 TeV')

mask_90 = TS_profile > tsmax - chi2_90_quantile
c_ul = max(counts[mask])
print("90% upper limit:", c_ul)
F_ul = ana.calculate_fluxmodel_scaling_factor(c_ul, [c_ul, Slope])
plt.axhline(
    tsmax - chi2_90_quantile, color="black", linestyle="dotted", label=r"90%"
)
# plt.text(0,tsmax-chi2_90_quantile-1.5,r'$F_{ul}=$'+str(F_ul)+' (GeV s cm^2)$^{-1}$ at 1 TeV')
plt.xlabel("Counts")
plt.ylabel("TS")
plt.savefig("Likelihood_profile.png", format="png", bbox_inches="tight")
plt.legend(loc="lower left")

bin_image = PictureProduct.from_file("Likelihood_profile.png")
from astropy.table import Table

data = [counts, TS_profile]
names = ("Counts", "TS")
spec_params = ODAAstropyTable(Table(data, names=names))

likelihood_profile_png = (
    bin_image  # http://odahub.io/ontology#ODAPictureProduct
)
likelihood_profile_table = (
    spec_params  # http://odahub.io/ontology#ODAAstropyTable
)

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_Likelihood_profile_likelihood_profile_png",
        "likelihood_profile_png_galaxy.output",
        likelihood_profile_png,
    )
)
_oda_outs.append(
    (
        "out_Likelihood_profile_likelihood_profile_table",
        "likelihood_profile_table_galaxy.output",
        likelihood_profile_table,
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
