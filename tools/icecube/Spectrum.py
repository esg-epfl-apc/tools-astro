#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import numpy as np
import scipy.stats
from matplotlib import pyplot as plt
from oda_api.data_products import ODAAstropyTable, PictureProduct
from oda_api.json import CustomJSONEncoder
from skyllh.analyses.i3.publicdata_ps.time_integrated_ps import create_analysis
from skyllh.core.config import Config
from skyllh.core.random import RandomStateService
from skyllh.core.source_model import PointLikeSource
from skyllh.datasets.i3.PublicData_10y_ps import create_dataset_collection

# src_name='NGC 1068' #http://odahub.io/ontology#AstrophysicalObject
# RA = 40.669622  # http://odahub.io/ontology#PointOfInterestRA
# DEC = -0.013294 # http://odahub.io/ontology#PointOfInterestDEC
src_name = "TXS 0506+056"  # http://odahub.io/ontology#AstrophysicalObject
RA = 77.35
DEC = 5.7

T1 = "2000-10-09T13:16:00.0"  # http://odahub.io/ontology#StartTime
T2 = "2022-10-10T13:16:00.0"  # http://odahub.io/ontology#EndTime

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

source = PointLikeSource(ra=np.deg2rad(RA), dec=np.deg2rad(DEC))

ana = create_analysis(cfg=cfg, datasets=datasets, source=source)

events_list = [data.exp for data in ana.data_list]
ana.initialize_trial(events_list)

rss = RandomStateService(seed=1)
(ts, x, status) = ana.unblind(minimizer_rss=rss)

print(f"TS = {ts:.3f}")
print(f'ns = {x["ns"]:.2f}')
print(f'gamma = {x["gamma"]:.2f}')

scaling_factor = ana.calculate_fluxmodel_scaling_factor(
    x["ns"], [x["ns"], x["gamma"]]
)
print(f"Flux scaling factor = {scaling_factor:.3e}")
print(
    f"that is {scaling_factor:.3e}"
    " (E/1000 GeV)^{-"
    f'{x["gamma"]:.2f}' + "} 1/(GeV s cm^2 sr)"
)

scaling_factor = ana.calculate_fluxmodel_scaling_factor(
    x["ns"], [x["ns"], x["gamma"]]
)
print(f"Flux scaling factor = {scaling_factor:.3e}")
print(
    f"that is {scaling_factor:.3e}"
    " (E/1000 GeV)^{-"
    f'{x["gamma"]:.2f}' + "} 1/(GeV s cm^2 sr)"
)

(llhratio_value, (grad_ns, grad_gamma)) = ana.llhratio.evaluate([14.58, 2.17])
print(f"llhratio_value = {llhratio_value:.3f}")
print(f"grad_ns = {grad_ns:.3f}")
print(f"grad_gamma = {grad_gamma:.3f}")

(ns_min, ns_max, ns_step) = (0, 80, 0.5)
(gamma_min, gamma_max, gamma_step) = (1.5, 4.0, 0.1)

ns_edges = np.linspace(ns_min, ns_max, int((ns_max - ns_min) / ns_step) + 1)
ns_vals = 0.5 * (ns_edges[1:] + ns_edges[:-1])

gamma_edges = np.linspace(
    gamma_min, gamma_max, int((gamma_max - gamma_min) / gamma_step + 1)
)
gamma_vals = 0.5 * (gamma_edges[1:] + gamma_edges[:-1])

delta_ts = np.empty((len(ns_vals), len(gamma_vals)), dtype=np.double)
for ns_i, ns in enumerate(ns_vals):
    for gamma_i, gamma in enumerate(gamma_vals):

        delta_ts[ns_i, gamma_i] = ana.calculate_test_statistic(
            llhratio_value, [14.58, 2.17]
        ) - ana.calculate_test_statistic(
            ana.llhratio.evaluate([ns, gamma])[0], [ns, gamma]
        )

# Determine the best fit ns and gamma values from the scan.
index_max = np.argmin(delta_ts)
ns_i_max = int(index_max / len(gamma_vals))
gamma_i_max = index_max % len(gamma_vals)
ns_best = ns_vals[ns_i_max]
gamma_best = gamma_vals[gamma_i_max]

# Determine the delta lambda value for the 95% quantile assuming a chi-sqaure
# distribution with 2 degrees of freedom (i.e. assuming Wilks theorem).
chi2_68_quantile = scipy.stats.chi2.ppf(0.68, df=2)
chi2_90_quantile = scipy.stats.chi2.ppf(0.90, df=2)
chi2_95_quantile = scipy.stats.chi2.ppf(0.95, df=2)

plt.figure(figsize=(8, 6))
# plt.pcolormesh(gamma_edges, ns_edges, delta_ts, cmap='nipy_spectral')
# cbar = plt.colorbar()
# cbar.set_label(r'$\Delta$TS')
cnt = plt.contour(
    gamma_vals, ns_vals, delta_ts, [chi2_68_quantile], colors="black"
)
cnt1 = plt.contour(
    gamma_vals, ns_vals, delta_ts, [chi2_90_quantile], colors="black"
)
cnt2 = plt.contour(
    gamma_vals, ns_vals, delta_ts, [chi2_95_quantile], colors="black"
)
plt.plot(gamma_best, ns_best, marker="x", color="black", ms=10)
plt.xlabel(r"$\gamma$")
plt.ylabel(r"$n_{\mathrm{s}}$")
plt.ylim(ns_min, ns_max)
plt.xlim(gamma_min, gamma_max)

cont = cnt.get_paths()[0].vertices
gammas = cont[:, 0]
norms = cont[:, 1]

F_norms = []
for i in range(len(norms)):
    F_norms.append(
        ana.calculate_fluxmodel_scaling_factor(norms[i], [norms[i], gammas[i]])
    )
Fbest = ana.calculate_fluxmodel_scaling_factor(ns_best, [ns_best, gamma_best])

plt.plot(gammas, F_norms)
plt.yscale("log")
plt.ylabel("Flux normalisation at 1 TeV, 1/(GeV s cm$^2$)")
plt.scatter([gamma_best], [Fbest], marker="x")
plt.xlabel("Slope")
plt.savefig("Confidence_range_68.png", format="png", bbox_inches="tight")

bin_image = PictureProduct.from_file("Confidence_range_68.png")
from astropy.table import Table

data = [gammas, F_norms]
names = ("Slopes", "Fnorms_1TeV[1/(Gev s cm2)")
spec_params = ODAAstropyTable(Table(data, names=names))

confidence_contour_png = (
    bin_image  # http://odahub.io/ontology#ODAPictureProduct
)
spec_params_astropy_table = (
    spec_params  # http://odahub.io/ontology#ODAAstropyTable
)

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_Spectrum_confidence_contour_png",
        "confidence_contour_png_galaxy.output",
        confidence_contour_png,
    )
)
_oda_outs.append(
    (
        "out_Spectrum_spec_params_astropy_table",
        "spec_params_astropy_table_galaxy.output",
        spec_params_astropy_table,
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
