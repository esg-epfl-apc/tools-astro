#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import wget
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from gammapy.data import (
    FixedPointingInfo,
    Observation,
    PointingMode,
    observatory_locations,
)
from gammapy.datasets import MapDataset, MapDatasetEventSampler
from gammapy.irf import load_irf_dict_from_file
from gammapy.makers import MapDatasetMaker
from gammapy.maps import Map, MapAxis
from gammapy.modeling.models import (
    FoVBackgroundModel,
    Models,
    SkyModel,
    TemplateSpatialModel,
)
from numpy import cos, pi, sqrt
from oda_api.data_products import NumpyDataProduct, PictureProduct
from oda_api.json import CustomJSONEncoder

get_ipython().run_cell_magic("bash", "", "git lfs install\ngit lfs pull\n")   # noqa: F821

# We simulate point source in wobble observaiton,
# 0.4 degree off-axis

# Exposure time in hours
Texp = 2.4  # http://odahub.io/ontology#TimeIntervalHours

file_path = "3d.fits"  # http://odahub.io/ontology#POSIXPath

file_url = ""  # http://odahub.io/ontology#String

# Source flux normalisaiton F0 in 1/(TeV cm2 s) at reference energy E0
# TODO: implement flux normalisation for fits input
F0 = 4e-13  # http://odahub.io/ontology#Float
E0 = 1.0  # http://odahub.io/ontology#Energy_TeV

Emax = 30  # http://odahub.io/ontology#Energy_TeV
Emin = 0.1  # http://odahub.io/ontology#Energy_TeV

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

if len(file_url) > 0:
    if file_url.startswith("file://"):
        file_path = fileurl[len("file://") :]
    else:
        file_path = wget.download(file_url)
file_path

cube_map = Map.read(file_path)
cube_map.geom

cube_map.geom.center_skydir

# source = SkyCoord.from_name(src_name, frame='icrs', parse=False, cache=True)
source = cube_map.geom.center_skydir
DEC = float(source.dec / u.deg)
RA = float(source.ra / u.deg)

# telescope pointing will be shifted slightly
cdec = cos(DEC * pi / 180.0)
pnt_RA = RA - 0.4 / cdec
pnt_DEC = DEC
pnt = SkyCoord(pnt_RA, pnt_DEC, unit="degree")

# telescope is pointing at a fixed position in ICRS for the observation
pointing = FixedPointingInfo(fixed_icrs=pnt, mode=PointingMode.POINTING)

location = observatory_locations["cta_south"]

# irfs = load_irf_dict_from_file(path / irf_filename)
filename = "data/Prod5-North-20deg-AverageAz-4LSTs09MSTs.180000s-v0.1.fits.gz"
irfs = load_irf_dict_from_file(filename)

livetime = Texp * u.hr
observation = Observation.create(
    obs_id=1001,
    pointing=pointing,
    livetime=livetime,
    irfs=irfs,
    location=location,
)
print(observation)

mid_energy = cube_map.data.shape[0] // 2
map0 = cube_map.slice_by_idx({"energy": mid_energy})
map0.plot()
plt.show()

max_rel_energy_error = 3
selected_n_bins_per_decade = 20  # n bins per decade

energy_axis = cube_map.geom.axes["energy"]
EminMap = energy_axis.edges[0]
EmaxMap = energy_axis.edges[-1]
stepE = energy_axis.edges[1] / energy_axis.edges[0]
nbins_per_decade = int(np.round(np.log(10) / np.log(stepE)))
Emin = EminMap / max_rel_energy_error
Emax = EmaxMap * max_rel_energy_error
nbins_per_decade, Emin, Emax

# energy_axis = MapAxis.from_energy_bounds(max_rel_energy_error*Emin*u.TeV, Emax*u.TeV, nbin=selected_n_bins_per_decade, per_decade=True)
energy_axis_true = MapAxis.from_energy_bounds(
    Emin, Emax, nbin=nbins_per_decade, per_decade=True, name="energy_true"
)  # TODO: get from geom
migra_axis = MapAxis.from_bounds(
    1.0 / max_rel_energy_error,
    max_rel_energy_error,
    nbin=150,
    node_type="edges",
    name="migra",
)
# TODO: get from geom

geom = cube_map.geom

# WcsGeom.create(
#     skydir=pointing.fixed_icrs,
#     width=(2, 2),
#     binsz=0.02,
#     frame="icrs",
#     axes=[energy_axis],
# )

empty = MapDataset.create(
    geom,
    energy_axis_true=energy_axis_true,
    migra_axis=migra_axis,
    name="my-dataset",
)
maker = MapDatasetMaker(selection=["exposure", "background", "psf", "edisp"])
dataset = maker.run(empty, observation)

Path("event_sampling").mkdir(exist_ok=True)
dataset.write("./event_sampling/dataset.fits", overwrite=True)

def GetBinSpectralModel(E, bins_per_decade=20, norm=1):
    amplitude = 1e-12 * u.Unit("cm-2 s-1") * norm
    from gammapy.modeling.models import GaussianSpectralModel

    sigma = (10 ** (1 / bins_per_decade) - 1) * E
    return GaussianSpectralModel(mean=E, sigma=sigma, amplitude=amplitude)

spec = cube_map.get_spectrum()
spec

spec.plot()  # this plot shows dN/dE * E

spec.data.shape, spec.data[spec.data.shape[0] // 2, 0, 0]

energy_bins = cube_map.geom.axes["energy"].center
len(energy_bins), float(np.max(energy_bins) / u.TeV)

Npart = 5000  # TODO update
n_events_reduction_factor = 1  # suppress flux factor

int_bin_flux = (
    spec.data.flatten()
)  # we don't have to multiply by energy_bins /u.TeV since spectrum is in already multiplied by E (see above)
int_bin_flux /= (
    Npart
    / 200000
    * np.max(int_bin_flux)
    * n_events_reduction_factor
    * 20
    / len(energy_bins)
)  # roughly 100 events
int_bin_flux

bin_models = []
for i, (flux, E) in enumerate(zip(int_bin_flux, energy_bins)):
    if flux == 0:
        continue
    spectral_model_delta = GetBinSpectralModel(
        E, norm=flux
    )  # normalizing here
    spacial_template_model = TemplateSpatialModel(
        cube_map.slice_by_idx({"energy": i}),
        filename=f"cube_bin{i}.fit",
        normalize=True,
    )
    sky_bin_model = SkyModel(
        spectral_model=spectral_model_delta,
        spatial_model=spacial_template_model,
        name=f"bin_{i}",
    )
    bin_models.append(sky_bin_model)

bkg_model = FoVBackgroundModel(dataset_name="my-dataset")
models = Models(bin_models + [bkg_model])

dataset.models = models
# print(dataset.models)

sampler = MapDatasetEventSampler(random_state=0)
events = sampler.run(dataset, observation)

print(f"Source events: {(events.table['MC_ID'] > 0).sum()}")
print(f"Background events: {(events.table['MC_ID'] == 0).sum()}")

for i in range(1, len(bin_models) + 1):
    n = (events.table["MC_ID"] == i).sum()
    if n > 1:
        print(f"\tmodel {i}: {n} events")

E = events.energy / u.TeV
ras = events.radec.ra.deg
decs = events.radec.dec.deg
# plt.hist(E,bins=np.logspace(-2,2,41))

mask = events.table["MC_ID"] > 0
plt.hist(E[mask], bins=np.logspace(-2, 2, 41), alpha=0.5, label="source")
mask = events.table["MC_ID"] == 0
plt.hist(E[mask], bins=np.logspace(-2, 2, 41), alpha=0.5, label="background")

plt.xscale("log")
plt.yscale("log")
plt.legend(loc="upper right")
plt.savefig("event_spectrum.png", format="png")

ROI = 1.0
pixsize = 0.02
Npix = int(2 * ROI / pixsize) + 1

mask = np.where(E > 1)
print(len(E[mask]), len(ras[mask]), len(decs[mask]))

plt.hist2d(
    ras[mask],
    decs[mask],
    bins=[
        np.linspace(pnt_RA - ROI / cdec, pnt_RA + ROI / cdec, Npix),
        np.linspace(DEC - ROI, DEC + ROI, Npix),
    ],
)
plt.colorbar()

print(f"Save events ...")
primary_hdu = fits.PrimaryHDU()
hdu_evt = fits.BinTableHDU(events.table)
hdu_gti = fits.BinTableHDU(dataset.gti.table, name="GTI")
hdu_all = fits.HDUList([primary_hdu, hdu_evt, hdu_gti])
hdu_all.writeto(f"./events.fits", overwrite=True)
####################

hdul = fits.open("events.fits")
T_exp = hdul["EVENTS"].header["ONTIME"]
events = hdul["EVENTS"].data

coords_s = SkyCoord(RA, DEC, unit="degree")
RA_bkg = pnt_RA - (RA - pnt_RA)
DEC_bkg = pnt_DEC - (DEC - pnt_DEC)
coords_b = SkyCoord(RA_bkg, DEC_bkg, unit="degree")

Es = events["ENERGY"]
ras = events["RA"]
decs = events["DEC"]
coords = SkyCoord(ras, decs, unit="degree")
seps_s = coords.separation(coords_s).deg
seps_b = coords.separation(coords_b).deg
coords_s = SkyCoord(RA, DEC, unit="degree")

hdul = fits.open(
    "Prod5-North-20deg-AverageAz-4LSTs09MSTs.180000s-v0.1.fits.gz"
)
Aeff = hdul["EFFECTIVE AREA"].data
th_min = Aeff["THETA_LO"]
th_min
Emin_irf = Aeff["ENERG_LO"][0]
Emax_irf = Aeff["ENERG_HI"][0]
E_irf = sqrt(Emin_irf * Emax_irf)
Emin_irf
Ebins_irf = np.concatenate((Emin_irf, [Emax_irf[-1]]))
A = Aeff["EFFAREA"][0, 0]

th_cut = 0.3
mask = seps_s < th_cut
E_s = Es[mask]
h1 = plt.hist(E_s, bins=Ebins_irf, alpha=0.5)
mask = seps_b < th_cut
E_b = Es[mask]
h2 = plt.hist(E_b, bins=Ebins_irf, alpha=0.5)
plt.xscale("log")
plt.yscale("log")
Ns = h1[0]
Nb = h2[0]
Src = Ns - Nb
Src_err = sqrt(Ns + Nb)
print(Src, Src_err)

# EE=E_irf[4:-3]
# Flux=Src/(Emax_irf-Emin_irf)*E_irf**2/(A*1e4)/T_exp
# Flux_err=Src_err/(Emax_irf-Emin_irf)*E_irf**2/(A*1e4)/T_exp
# Flux_err=Flux_err+100*(Flux_err==0)
# FFlux=Flux[4:-3]
# FFlux_err=Flux_err[4:-3]

# from scipy.optimize import curve_fit
# def PL(x,Norm):
#     return Norm*(x/E0)**(2-Gamma)*exp(-tau(x))

# plt.errorbar(EE,FFlux,yerr=FFlux_err)

# popt, pcov = curve_fit(PL, EE, FFlux, sigma=FFlux_err)
# F0_best=popt[0]
# y=F0_best*(EE/E0)**(2-Gamma)*exp(-tau(EE))
# plt.yscale('log')
# plt.xscale('log')
# F0_err=sqrt(pcov[0,0])
# SN=F0_best/F0_err
# plt.plot(EE,y,label='S/N='+str(SN))
# plt.legend(loc='upper right')
# plt.ylim(1e-14,1e-9)
# plt.savefig('spectrum.png',format='png')

# theta2 plot
thbin = 0.1
nb = int(1 / thbin**2)
print(nb)
bins = np.linspace(0, 1, nb)
h1 = plt.hist(seps_s**2, bins=bins, alpha=0.5)
h2 = plt.hist(seps_b**2, bins=bins, alpha=0.5)
plt.axvline(th_cut**2)
plt.xlim(0, 0.3)

cts_s = sum(h1[0][:10])
cts_b = sum(h2[0][:10])
SN = (cts_s - cts_b) / sqrt(cts_b)
print((cts_s - cts_b) / sqrt(cts_b))
plt.text(0.15, max(h1[0][:10]), "S/N=" + str(SN))

plt.savefig("theta2.png", format="png")

events_fits = NumpyDataProduct.from_fits_file("events.fits")
# spectrum_png = PictureProduct.from_file('spectrum.png')
theta2_png = PictureProduct.from_file("theta2.png")

theta2 = theta2_png  # http://odahub.io/ontology#ODAPictureProduct
# spectrum = spectrum_png # http://odahub.io/ontology#ODAPictureProduct
events_fits = events_fits  # https://odahub.io/ontology/#Spectrum

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    ("out_model_CTA_events_from_file_theta2", "theta2_galaxy.output", theta2)
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
_simple_outs = []
_simple_outs.append(
    (
        "out_model_CTA_events_from_file_events_fits",
        "events_fits_galaxy.output",
        events_fits,
    )
)
_numpy_available = True

for _outn, _outfn, _outv in _simple_outs:
    _galaxy_outfile_name = os.path.join(_galaxy_wd, _outfn)
    if isinstance(_outv, str) and os.path.isfile(_outv):
        shutil.move(_outv, _galaxy_outfile_name)
        _galaxy_meta_data[_outn] = {"ext": "_sniff_"}
    elif _numpy_available and isinstance(_outv, np.ndarray):
        with open(_galaxy_outfile_name, "wb") as fd:
            np.savez(fd, _outv)
        _galaxy_meta_data[_outn] = {"ext": "npz"}
    else:
        with open(_galaxy_outfile_name, "w") as fd:
            json.dump(_outv, fd)
        _galaxy_meta_data[_outn] = {"ext": "expression.json"}

with open(os.path.join(_galaxy_wd, "galaxy.json"), "w") as fd:
    json.dump(_galaxy_meta_data, fd)
print("*** Job finished successfully ***")
