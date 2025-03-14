#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import numpy as np
import scipy.stats
from astropy.time import Time
from matplotlib import pyplot as plt
from oda_api.data_products import ODAAstropyTable, PictureProduct
from oda_api.json import CustomJSONEncoder
from skyllh.analyses.i3.publicdata_ps.time_dependent_ps import create_analysis
from skyllh.core.config import Config
from skyllh.core.random import RandomStateService

# from skyllh.analyses.i3.publicdata_ps.time_integrated_ps import create_analysis
from skyllh.core.source_model import PointLikeSource
from skyllh.datasets.i3.PublicData_10y_ps import create_dataset_collection

# src_name='NGC 1068' #http://odahub.io/ontology#AstrophysicalObject
RA = 40.669622  # http://odahub.io/ontology#PointOfInterestRA
DEC = 0.013294  # http://odahub.io/ontology#PointOfInterestDEC
Slope = 3.0  # http://odahub.io/ontology#Float

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in ["RA", "DEC", "Slope"]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

if os.path.exists("20210126_PS-IC40-IC86_VII.zip") == False:
    get_ipython().system(   # noqa: F821
        "wget http://icecube.wisc.edu/data-releases/20210126_PS-IC40-IC86_VII.zip"
    )
    get_ipython().system("unzip 20210126_PS-IC40-IC86_VII.zip")   # noqa: F821

data_dir = os.getcwd() + "/icecube_10year_ps/"

cfg = Config()

dsc = create_dataset_collection(cfg=cfg, base_path=data_dir)
periods = [
    "IC40",
    "IC59",
    "IC79",
    "IC86_I",
    "IC86_II",
    "IC86_III",
    "IC86_IV",
    "IC86_V",
    "IC86_VI",
    "IC86_VII",
]

Tstarts = [
    "2008-04-06T00:00:00.0",
    "2009-05-20T00:00:00.0",
    "2010-06-01T00:00:00.0",
    "2011-05-13T00:00:00.0",
    "2012-04-26T00:00:00.0",
    "2013-05-02T00:00:00.0",
    "2014-04-10T00:00:00.0",
    "2015-04-24T00:00:00.0",
    "2016-05-20T00:00:00.0",
    "2017-05-18T00:00:00.0",
]
Tstops = [
    "2009-05-20T00:00:00.0",
    "2010-05-31T00:00:00.0",
    "2011-05-13T00:00:00.0",
    "2012-05-15T00:00:00.0",
    "2013-05-02T00:00:00.0",
    "2014-05-06T00:00:00.0",
    "2015-05-18T00:00:00.0",
    "2016-05-20T00:00:00.0",
    "2017-05-18T00:00:00.0",
    "2018-07-08T00:00:00.0",
]
Tstarts = Time(Tstarts, format="isot", scale="utc").mjd
Tstops = Time(Tstops, format="isot", scale="utc").mjd

source = PointLikeSource(ra=np.deg2rad(RA), dec=np.deg2rad(DEC))

chi2_68_quantile = scipy.stats.chi2.ppf(0.68, df=1)
chi2_90_quantile = scipy.stats.chi2.ppf(0.90, df=1)
chi2_95_quantile = scipy.stats.chi2.ppf(0.95, df=1)

Fmin_68 = np.zeros(len(Tstarts))
Fmax_68 = np.zeros(len(Tstarts))
Fmax_90 = np.zeros(len(Tstarts))
Fbest = np.zeros(len(Tstarts))
for ind in range(len(Tstarts)):
    print(periods[ind])
    datasets = dsc[periods[ind],]
    ana = create_analysis(
        cfg=cfg,
        datasets=datasets,
        source=source,
        refplflux_gamma=2.0,
        box={"start": Tstarts[ind], "stop": Tstops[ind]},
    )
    events_list = [data.exp for data in ana.data_list]
    ana.initialize_trial(events_list)
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
    tsmax = max(TS_profile)
    cbest = counts[np.argmax(TS_profile)]
    print(max(TS_profile), cbest)
    mask = TS_profile > tsmax - chi2_68_quantile
    cmin = min(counts[mask])
    cmax = max(counts[mask])
    print(cmin, cmax)
    Fmin_68[ind] = (
        ana.calculate_fluxmodel_scaling_factor(cmin, [cmin, Slope]) * 1e3 * 3
    )  # 1/(GeV s cm^2 sr) at 1e3 GeV
    Fmax_68[ind] = (
        ana.calculate_fluxmodel_scaling_factor(cmax, [cmin, Slope]) * 1e3 * 3
    )  # 1/(GeV s cm^2 sr) at 1e3 GeV
    Fbest[ind] = (
        ana.calculate_fluxmodel_scaling_factor(cbest, [cmin, Slope]) * 1e3 * 3
    )  # 1/(GeV s cm^2 sr) at 1e3 GeV
    mask = TS_profile > tsmax - chi2_90_quantile
    cmin = min(counts[mask])
    cmax = max(counts[mask])
    Fmax_90[ind] = (
        ana.calculate_fluxmodel_scaling_factor(cmax, [cmin, Slope]) * 1e3 * 3
    )  # 1/(GeV s cm^2 sr) at 1e3 GeV
    print("All-flavour:", Fmin_68[ind], Fmax_68[ind], Fbest[ind], Fmax_90[ind])

Tmeans = (Tstarts + Tstops) / 2.0
Tdels = Tstops - Tstarts
plt.errorbar(
    Tmeans,
    Fbest,
    yerr=[Fbest - Fmin_68, Fmax_68 - Fbest],
    xerr=Tdels / 2.0,
    color="black",
    linewidth=2,
)
plt.errorbar(
    Tmeans,
    Fmax_90,
    xerr=Tdels / 2.0,
    yerr=Fmax_90 / 5.0,
    uplims=np.ones(len(Tmeans)),
    linestyle="none",
    color="black",
    alpha=0.5,
    linewidth=2,
)
plt.xlabel("Time, MJD")
plt.ylabel(r"$E^2dN/dE$ at 1 TeV [TeV/(cm$^2$ s)], all flavours")
plt.savefig("Lightcurve.png", format="png", bbox_inches="tight")

bin_image = PictureProduct.from_file("Lightcurve.png")
from astropy.table import Table

data = [Tstarts, Tstops, Fbest, Fmin_68, Fmax_68, Fmax_90]
names = (
    "Tstart[MJD]",
    "Tstop[MJD]",
    "F[TeV/cm2s]",
    "F_min_68[TeV/cm2s]",
    "F_max_68[TeV/cm2s]",
    "F_max_90[TeV/cm2s]",
)
lightcurve = ODAAstropyTable(Table(data, names=names))

get_ipython().system(" rm -rfv {data_dir}")   # noqa: F821

lightcurve_png = bin_image  # http://odahub.io/ontology#ODAPictureProduct
lightcurve_table = lightcurve  # http://odahub.io/ontology#ODAAstropyTable

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_Lightcurve_lightcurve_png",
        "lightcurve_png_galaxy.output",
        lightcurve_png,
    )
)
_oda_outs.append(
    (
        "out_Lightcurve_lightcurve_table",
        "lightcurve_table_galaxy.output",
        lightcurve_table,
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
