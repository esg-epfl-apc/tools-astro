#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import numpy as np
import scipy.stats
from matplotlib import pyplot as plt
from oda_api.api import ProgressReporter
from oda_api.data_products import ODAAstropyTable, PictureProduct
from oda_api.json import CustomJSONEncoder

# src_name='NGC 1068' #http://odahub.io/ontology#AstrophysicalObject
RA = 40.669622  # http://odahub.io/ontology#PointOfInterestRA
DEC = 10.013294  # http://odahub.io/ontology#PointOfInterestDEC
Spectrum_type = "Free_slope"  # http://odahub.io/ontology#String ; oda:allowed_value "Fixed_slope","Free_slope"
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

source = PointLikeSource(ra=np.deg2rad(RA), dec=np.deg2rad(DEC))

ana = create_analysis(cfg=cfg, datasets=datasets, source=source)

events_list = [data.exp for data in ana.data_list]
ana.initialize_trial(events_list)

pr = ProgressReporter()
pr.report_progress(stage="Progress", progress=5.0)

rss = RandomStateService(seed=1)
(ts, x, status) = ana.unblind(minimizer_rss=rss)

print(f"TS = {ts:.3f}")
print(f'ns = {x["ns"]:.2f}')
print(f'gamma = {x["gamma"]:.2f}')

chi2_68_quantile = scipy.stats.chi2.ppf(0.68, df=2)
chi2_90_quantile = scipy.stats.chi2.ppf(0.90, df=2)
chi2_95_quantile = scipy.stats.chi2.ppf(0.95, df=2)
chi2_68_quantile, chi2_90_quantile, chi2_95_quantile

if Spectrum_type == "Fixed_slope":
    chi2_68_quantile = scipy.stats.chi2.ppf(0.68, df=1)
    chi2_90_quantile = scipy.stats.chi2.ppf(0.90, df=1)
    chi2_95_quantile = scipy.stats.chi2.ppf(0.95, df=1)
    TS_profile = []
    counts = np.linspace(0, 200, 200)
    for n in counts:
        TS_profile.append(2 * ana.llhratio.evaluate([n, Slope])[0])
    # plt.plot(counts,TS_profile)
    tsmax = max(TS_profile)
    print(max(TS_profile))
    cbest = counts[np.argmax(TS_profile)]
    mask = TS_profile > tsmax - chi2_68_quantile
    cmin = min(counts[mask])
    cmax = max(counts[mask])
    print(cmin, cmax)
    Fmin_68 = ana.calculate_fluxmodel_scaling_factor(
        cmin, [cmin, Slope]
    )  # 1/(GeV s cm^2 sr) at 1e3 GeV
    Fmax_68 = ana.calculate_fluxmodel_scaling_factor(
        cmax, [cmin, Slope]
    )  # 1/(GeV s cm^2 sr) at 1e3 GeV
    mask = TS_profile > tsmax - chi2_90_quantile
    cmin = min(counts[mask])
    cmax = max(counts[mask])
    Fmin_90 = ana.calculate_fluxmodel_scaling_factor(
        cmin, [cmin, Slope]
    )  # 1/(GeV s cm^2 sr) at 1e3 GeV
    Fmax_90 = ana.calculate_fluxmodel_scaling_factor(
        cmax, [cmin, Slope]
    )  # 1/(GeV s cm^2 sr) at 1e3 GeV
    Fbest = ana.calculate_fluxmodel_scaling_factor(
        cbest, [cmin, Slope]
    )  # 1/(GeV s cm^2 sr) at 1e3 GeV
    x = np.logspace(0, 2, 10)  # energy in TeV
    ymin_68 = (
        3 * Fmin_68 * (x / 1) ** (2 - Slope) * 1e3
    )  # in TeV/cm2s, all flavours
    ymax_68 = 3 * Fmax_68 * (x / 1) ** (2 - Slope) * 1e3
    ymax_90 = 3 * Fmax_90 * (x / 1) ** (2 - Slope) * 1e3
    ybest = 3 * Fbest * (x / 1) ** (2 - Slope) * 1e3

    if np.amax(TS_profile) > chi2_95_quantile:
        plt.fill_between(x, ymin_68, ymax_68, alpha=0.5, label="68% error")
        plt.plot(x, ybest, color="black")
    plt.plot(x, ymax_90, color="black", linewidth=4, label="90% UL")

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("$E$, TeV")
    # plt.ylim(1e-14,3e-10)
    plt.ylabel("$E^2 dN/dE$, TeV/(cm$^2$s), all flavours")
    plt.savefig("Spectrum.png", format="png", bbox_inches="tight")

if Spectrum_type == "Free_slope":
    chi2_68_quantile = scipy.stats.chi2.ppf(0.68, df=2)
    chi2_90_quantile = scipy.stats.chi2.ppf(0.90, df=1)
    chi2_95_quantile = scipy.stats.chi2.ppf(0.95, df=1)
    TS_map = np.zeros((50, 100))
    counts = np.linspace(0, 200, 100)
    Slopes = np.linspace(1, 5, 50)
    tsbest = 0
    slope_best = 0
    n_best = 0
    c_ul = np.zeros(len(Slopes))
    for i, Slope in enumerate(Slopes):
        pr.report_progress(stage="Progress", progress=(i / len(Slopes)))
        for j, n in enumerate(counts):
            TS_map[i, j] = 2 * ana.llhratio.evaluate([n, Slope])[0]
            if TS_map[i, j] > tsbest:
                tsbest = TS_map[i, j]
                slope_best = Slope
                n_best = n
        tsmax = np.max(TS_map[i, :])
        mask_90 = TS_map[i, :] > tsmax - chi2_90_quantile
        c_ul[i] = max(counts[mask_90])
    plt.plot(Slopes, c_ul)
    # plt.imshow(np.transpose(TS_map),origin='lower',extent=[Slopes[0],Slopes[-1],counts[0],counts[-1]],aspect=0.03,vmin=0)
    # plt.colorbar()
    # plt.scatter(slope_best,n_best,marker='x',color='red')
    print(Slope, n, tsbest, n_best, slope_best)
    cnt_68 = plt.contour(
        Slopes,
        counts,
        np.transpose(TS_map),
        [tsbest - chi2_68_quantile],
        colors="red",
    )
    cnt_90 = plt.contour(
        Slopes,
        counts,
        np.transpose(TS_map),
        [tsbest - chi2_90_quantile],
        colors="red",
    )
    cont_68 = cnt_68.get_paths()[0].vertices
    gammas_68 = cont_68[:, 0]
    norms_68 = cont_68[:, 1]

    cont_90 = cnt_90.get_paths()[0].vertices
    gammas_90 = cont_90[:, 0]
    norms_90 = cont_90[:, 1]

    F_norms_68 = []
    for i in range(len(norms_68)):
        F_norms_68.append(
            ana.calculate_fluxmodel_scaling_factor(
                norms_68[i], [norms_68[i], gammas_68[i]]
            )
        )
    Fbest = ana.calculate_fluxmodel_scaling_factor(
        n_best, [n_best, slope_best]
    )

    F_norms_90 = []
    for i in range(len(norms_90)):
        F_norms_90.append(
            ana.calculate_fluxmodel_scaling_factor(
                norms_90[i], [norms_90[i], gammas_90[i]]
            )
        )
    plt.plot(gammas_90, F_norms_90)
    F_norms_90 = []
    for i in range(len(Slopes)):
        F_norms_90.append(
            ana.calculate_fluxmodel_scaling_factor(
                c_ul[i], [c_ul[i], Slopes[i]]
            )
        )

    plt.figure()
    plt.plot(Slopes, F_norms_90)
    print(F_norms_90)
    plt.plot(gammas_68, F_norms_68)
    plt.yscale("log")
    plt.ylabel("Flux normalisation at 1 TeV, 1/(GeV s cm$^2$)")
    plt.scatter([slope_best], [Fbest], marker="x")
    plt.xlabel("Slope")
    plt.savefig("Confidence_range.png", format="png", bbox_inches="tight")

if Spectrum_type == "Free_slope":
    x = np.logspace(-1, 6, 50)

    # 3.690e-15 (E/1000 GeV)^{-2.33} 1/(GeV s cm^2 sr)
    ymax_68 = np.zeros(len(x))
    ymin_68 = np.ones(len(x))
    for i in range(len(gammas_68)):
        y = (
            3 * F_norms_68[i] * (x / 1.0) ** (-gammas_68[i]) * x**2 * 1e3
        )  # TeV *(TeV/GeV)/cm2 s
        ymax_68 = np.maximum(ymax_68, y)
        ymin_68 = np.minimum(ymin_68, y)

    ymax_90 = np.zeros(len(x))
    for i in range(len(Slopes)):
        y = (
            3 * F_norms_90[i] * (x / 1.0) ** (-Slopes[i]) * x**2 * 1e3
        )  # TeV *(TeV/GeV)/cm2 s
        ymax_90 = np.maximum(ymax_90, y)

        # plt.plot(x,y)

    ybest = 3 * Fbest * (x / 1.0) ** (-slope_best) * x**2 * 1e3
    if np.amax(TS_map) > chi2_95_quantile:
        plt.fill_between(x, ymin_68, ymax_68, alpha=0.5, label="68% error")
        plt.plot(x, ybest, color="black")
    plt.plot(x, ymax_90, color="black", linewidth=4, label="90% UL")

    plt.xscale("log")
    plt.yscale("log")
    plt.legend(loc="upper right")
    plt.xlabel("$E$, TeV")
    plt.ylabel("$E^2 dN/dE$, TeV/(cm$^2$s), all flavours")
    plt.ylim(1e-14, 3e-10)
    plt.savefig("Spectrum.png", format="png", bbox_inches="tight")

bin_image = PictureProduct.from_file("Spectrum.png")
from astropy.table import Table

data = [x, ybest, ymin_68, ymax_68, ymax_90]
names = (
    "Energy[TeV]",
    "F_best[TeV/cm2s]",
    "F_min_68[TeV/cm2s]",
    "F_max_68[TeV/cm2s]",
    "F_max_90[TeV/cm2s]",
)
spec_params = ODAAstropyTable(Table(data, names=names))

spectrum_png = bin_image  # http://odahub.io/ontology#ODAPictureProduct
spectrum_table = spec_params  # http://odahub.io/ontology#ODAAstropyTable

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
