#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import Angle, SkyCoord
from astropy.time import Time
from gammapy.data import DataStore
from gammapy.datasets import (
    Datasets,
    FluxPointsDataset,
    SpectrumDataset,
    SpectrumDatasetOnOff,
)
from gammapy.makers import (
    ReflectedRegionsBackgroundMaker,
    SafeMaskMaker,
    SpectrumDatasetMaker,
)

# from gammapy.makers.utils import make_theta_squared_table
from gammapy.maps import MapAxis, RegionGeom, WcsGeom
from oda_api.data_products import ODAAstropyTable, PictureProduct
from oda_api.json import CustomJSONEncoder
from regions import CircleSkyRegion

hess_data = "gammapy-datasets/1.1/hess-dl3-dr1/"
if not (os.path.exists(hess_data)):
    get_ipython().system("gammapy download datasets")   # noqa: F821

data_store = DataStore.from_dir(hess_data)

# src_name='Crab' #http://odahub.io/ontology#AstrophysicalObject
# RA = 83.628700  # http://odahub.io/ontology#PointOfInterestRA
# DEC = 22.014700 # http://odahub.io/ontology#PointOfInterestDEC
src_name = "PKS 2155-304"
RA = 329.716938  # http://odahub.io/ontology#PointOfInterestRA
DEC = -30.225588  # http://odahub.io/ontology#PointOfInterestDEC
T1 = "2000-10-09T13:16:00.0"  # http://odahub.io/ontology#StartTime
T2 = "2022-10-10T13:16:00.0"  # http://odahub.io/ontology#EndTime
Radius = 2.5  # http://odahub.io/ontology#AngleDegrees
R_s = 0.5  # http://odahub.io/ontology#AngleDegrees

Emin = 100.0  # http://odahub.io/ontology#Energy_GeV
Emax = 10000.0  # http://odahub.io/ontology#Energy_GeV
NEbins = 20  # http://odahub.io/ontology#Integer

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

T1 = Time(T1, format="isot", scale="utc").mjd
T2 = Time(T2, format="isot", scale="utc").mjd

dates1 = data_store.obs_table["DATE-OBS"]
dates2 = data_store.obs_table["DATE-END"]
times1 = data_store.obs_table["TIME-OBS"]
times2 = data_store.obs_table["TIME-END"]
OBSIDs = data_store.obs_table["OBS_ID"]
Tstart = []
Tstop = []
for i in range(len(dates1)):
    Tstart.append(
        Time(dates1[i] + "T" + times1[i], format="isot", scale="utc").mjd
    )
    Tstop.append(
        Time(dates2[i] + "T" + times2[i], format="isot", scale="utc").mjd
    )

RA_pnts = np.array(data_store.obs_table["RA_PNT"])
DEC_pnts = np.array(data_store.obs_table["DEC_PNT"])

Coords_s = SkyCoord(RA, DEC, unit="degree")
COORDS_pnts = SkyCoord(RA_pnts, DEC_pnts, unit="degree")
seps = COORDS_pnts.separation(Coords_s).deg

mask = np.where((seps < Radius) & (Tstart > T1) & (Tstop < T2))[0]
OBSlist = []
obs_ids = OBSIDs[mask]
if len(obs_ids) == 0:
    message = "No data found"
    raise RuntimeError("No data found")
obs_ids

observations = data_store.get_observations(obs_ids)

target_position = Coords_s
on_region_radius = Angle(str(R_s) + " deg")
on_region = CircleSkyRegion(center=target_position, radius=on_region_radius)
skydir = target_position.galactic
geom = WcsGeom.create(
    npix=(150, 150), binsz=0.05, skydir=skydir, proj="TAN", frame="icrs"
)

Emin = 100.0  # http://odahub.io/ontology#Energy_GeV
Emax = 10000.0  # http://odahub.io/ontology#Energy_GeV
NEbins = 20  # http://odahub.io/ontology#Integer

energy_axis = MapAxis.from_energy_bounds(
    Emin * 1e-3,
    Emax * 1e-3,
    nbin=NEbins,
    per_decade=True,
    unit="TeV",
    name="energy",
)
energy_axis_true = MapAxis.from_energy_bounds(
    0.05, 100, nbin=20, per_decade=True, unit="TeV", name="energy_true"
)

geom = RegionGeom.create(region=on_region, axes=[energy_axis])
dataset_empty = SpectrumDataset.create(
    geom=geom, energy_axis_true=energy_axis_true
)

dataset_maker = SpectrumDatasetMaker(
    containment_correction=True, selection=["counts", "exposure", "edisp"]
)
bkg_maker = ReflectedRegionsBackgroundMaker()
safe_mask_masker = SafeMaskMaker(methods=["aeff-max"], aeff_percent=10)

datasets = Datasets()

for obs_id, observation in zip(obs_ids, observations):
    dataset = dataset_maker.run(
        dataset_empty.copy(name=str(obs_id)), observation
    )
    dataset_on_off = bkg_maker.run(dataset, observation)
    # dataset_on_off = safe_mask_masker.run(dataset_on_off, observation)
    datasets.append(dataset_on_off)

print(datasets)

from pathlib import Path

path = Path("spectrum_analysis")
path.mkdir(exist_ok=True)

for dataset in datasets:
    dataset.write(
        filename=path / f"obs_{dataset.name}.fits.gz", overwrite=True
    )

datasets = Datasets()

for obs_id in obs_ids:
    filename = path / f"obs_{obs_id}.fits.gz"
    datasets.append(SpectrumDatasetOnOff.read(filename))

from gammapy.modeling import Fit
from gammapy.modeling.models import ExpCutoffPowerLawSpectralModel, SkyModel

spectral_model = ExpCutoffPowerLawSpectralModel(
    amplitude=1e-12 * u.Unit("cm-2 s-1 TeV-1"),
    index=2,
    lambda_=0.1 * u.Unit("TeV-1"),
    reference=1 * u.TeV,
)
model = SkyModel(spectral_model=spectral_model, name="crab")

datasets.models = [model]

fit_joint = Fit()
result_joint = fit_joint.run(datasets=datasets)

# we make a copy here to compare it later
model_best_joint = model.copy()

print(result_joint)

display(result_joint.models.to_parameters_table())   # noqa: F821

e_min, e_max = Emin * 1e-3, Emax * 1e-3
energy_edges = np.geomspace(e_min, e_max, NEbins) * u.TeV

from gammapy.estimators import FluxPointsEstimator

fpe = FluxPointsEstimator(
    energy_edges=energy_edges, source="crab", selection_optional="all"
)
flux_points = fpe.run(datasets=datasets)

flux_points_dataset = FluxPointsDataset(
    data=flux_points, models=model_best_joint
)
flux_points_dataset.plot_fit()
# plt.show()
plt.savefig("Spectrum.png", format="png", bbox_inches="tight")

res = flux_points.to_table(sed_type="dnde", formatted=True)
np.array(res["dnde"])

bin_image = PictureProduct.from_file("Spectrum.png")
from astropy.table import Table

Emean = np.array(res["e_ref"])
Emin = np.array(res["e_min"])
Emax = np.array(res["e_max"])
flux = Emean**2 * np.array(res["dnde"])
flux_err = Emean**2 * np.array(res["dnde_err"])
data = [Emean, Emin, Emax, flux, flux_err]
names = (
    "Emean[TeV]",
    "Emin[TeV]",
    "Emax[TeV]",
    "Flux[TeV/cm2s]",
    "Flux_error[TeV/cm2s]",
)
spec = ODAAstropyTable(Table(data, names=names))

picture_png = bin_image  # http://odahub.io/ontology#ODAPictureProduct
spectrum_astropy_table = spec  # http://odahub.io/ontology#ODAAstropyTable

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_Spectrum_gammapy_picture_png",
        "picture_png_galaxy.output",
        picture_png,
    )
)
_oda_outs.append(
    (
        "out_Spectrum_gammapy_spectrum_astropy_table",
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
