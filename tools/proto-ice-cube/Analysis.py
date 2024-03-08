#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
os.system("pip install --upgrade skyllh")
import numpy as np
import scipy.stats
from matplotlib import pyplot as plt
from skyllh.analyses.i3.publicdata_ps.time_integrated_ps import create_analysis
from skyllh.core.config import Config
from skyllh.core.random import RandomStateService
from skyllh.core.source_model import PointLikeSource
from skyllh.datasets.i3.PublicData_10y_ps import create_dataset_collection

src_name = "NGC 1068"  # http://odahub.io/ontology#AstrophysicalObject
RA = 40.669622  # http://odahub.io/ontology#PointOfInterestRA
DEC = -0.013294  # http://odahub.io/ontology#PointOfInterestDEC

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
plt.pcolormesh(gamma_edges, ns_edges, delta_ts, cmap="nipy_spectral")
cbar = plt.colorbar()
cbar.set_label(r"$\Delta$TS")
plt.contour(
    gamma_vals, ns_vals, delta_ts, [chi2_68_quantile], colors="#FFFFFF"
)
plt.contour(
    gamma_vals, ns_vals, delta_ts, [chi2_90_quantile], colors="#AAAAAA"
)
plt.contour(
    gamma_vals, ns_vals, delta_ts, [chi2_95_quantile], colors="#444444"
)
plt.plot(gamma_best, ns_best, marker="x", color="white", ms=10)
plt.xlabel(r"$\gamma$")
plt.ylabel(r"$n_{\mathrm{s}}$")
plt.ylim(ns_min, ns_max)
plt.xlim(gamma_min, gamma_max)

from skyllh.core.timing import TimeLord
from skyllh.core.utils.analysis import create_trial_data_file

tl = TimeLord()

ncpu = int(os.cpu_count() / 2.0)
ncpu

rss = RandomStateService(seed=1)
(_, _, _, trials) = create_trial_data_file(
    ana=ana,
    rss=rss,
    n_trials=1e3,
    mean_n_sig=0,
    pathfilename=os.getcwd() + "/bkg_trails.npy",
    ncpu=3,
    tl=tl,
)
print(tl)

(h, be) = np.histogram(
    trials["ts"], bins=np.arange(0, np.max(trials["ts"]) + 0.1, 0.1)
)
plt.plot(
    0.5 * (be[:-1] + be[1:]), h, drawstyle="steps-mid", label="background"
)
plt.vlines(ts, 1, np.max(h), label=f"TS(TXS 0506+056)={ts:.3f}", color="red")
plt.yscale("log")
plt.xlabel("TS")
plt.ylabel("#trials per bin")
plt.legend()

from skyllh.core.utils.analysis import extend_trial_data_file

tl = TimeLord()
rss = RandomStateService(seed=2)
trials = extend_trial_data_file(
    ana=ana,
    rss=rss,
    n_trials=1e2,
    trial_data=trials,
    pathfilename=os.getcwd() + "/bkg_trails.npy",
    ncpu=3,
    tl=tl,
)

minus_log10_pval = -np.log10(len(trials[trials["ts"] > ts]) / len(trials))
print(f"-log10(p_local) = {minus_log10_pval:.2f}")

len(trials[trials["ts"] > ts])

# output gathering
_galaxy_meta_data = {}

with open(os.path.join(_galaxy_wd, "galaxy.json"), "w") as fd:
    json.dump(_galaxy_meta_data, fd)
print("*** Job finished successfully ***")
