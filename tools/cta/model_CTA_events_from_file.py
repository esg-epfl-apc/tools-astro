#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil
import sys

from oda_api.json import CustomJSONEncoder

sys.path.append(".")
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy import wcs
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
from oda_api.data_products import (
    BinaryProduct,
    ImageDataProduct,
    PictureProduct,
)

get_ipython().run_cell_magic(   # noqa: F821
    "bash",
    "",
    'rm -r IRFS | echo "Ok"\nmkdir IRFS\ncd IRFS\nwget https://zenodo.org/records/5499840/files/cta-prod5-zenodo-fitsonly-v0.1.zip\nunzip cta-prod5-zenodo-fitsonly-v0.1.zip\ncd fits\nfor fn in *.gz ; do tar -zxvf $fn; done \n',
)

# not for run on Galaxy
# %%bash
# git lfs install
# git lfs pull

# We simulate point source in wobble observaiton,
# 0.4 degree off-axis

# Exposure time in hours
Texp = 2.4  # http://odahub.io/ontology#TimeIntervalHours

file_path = "3d.fits"  # http://odahub.io/ontology#POSIXPath

# file_url='' # http://odahub.io/ontology#String

# Source flux normalisaiton F0 in 1/(TeV cm2 s) at reference energy E0
# TODO: implement flux normalisation for fits input
# F0=4e-13 # http://odahub.io/ontology#Float
# E0=1. # http://odahub.io/ontology#Energy_TeV

Emax = 30  # http://odahub.io/ontology#Energy_TeV
Emin = 0.1  # http://odahub.io/ontology#Energy_TeV

norm_cm2_TeV_s = 1e-12  # http://odahub.io/ontology#Float
norm_energy = 1.0  # http://odahub.io/ontology#Energy_TeV

pointing_shift = 0.2  # http://odahub.io/ontology#AngleDegrees

pixsize = 0.1  # http://odahub.io/ontology#AngleDegrees

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

print("loading " + file_path)
cube_map = Map.read(file_path)
cube_map.geom

cube_map.geom.center_skydir

print("locating source")
# source = SkyCoord.from_name(src_name, frame='icrs', parse=False, cache=True)
source = cube_map.geom.center_skydir
DEC = float(source.dec / u.deg)
RA = float(source.ra / u.deg)

# telescope pointing will be shifted slightly
cdec = cos(DEC * pi / 180.0)
RA_pnt = RA - pointing_shift / cdec
DEC_pnt = DEC
pnt = SkyCoord(RA_pnt, DEC_pnt, unit="degree")

# telescope is pointing at a fixed position in ICRS for the observation
pointing = FixedPointingInfo(fixed_icrs=pnt, mode=PointingMode.POINTING)

location = observatory_locations["cta_south"]

print("loading IRFs")

# irfs = load_irf_dict_from_file(path / irf_filename)
# filename = "data/Prod5-North-20deg-AverageAz-4LSTs09MSTs.180000s-v0.1.fits.gz"
irfs_filename = (
    "IRFS/fits/Prod5-North-20deg-AverageAz-4LSTs09MSTs.180000s-v0.1.fits.gz"
)
irfs = load_irf_dict_from_file(irfs_filename)

print("Creating observation")
livetime = Texp * u.hr
observation = Observation.create(
    obs_id=1001,
    pointing=pointing,
    livetime=livetime,
    irfs=irfs,
    location=location,
)
print(observation)

# print('Sowing map')
# mid_energy = cube_map.data.shape[0]//2
# map0 = cube_map.slice_by_idx({"energy": mid_energy})
# map0.plot()
# plt.show()

def GetBinSpectralModel(
    E, bins_per_decade=20, amplitude=1e-12 * u.Unit("cm-2 s-1")
):
    # amplitude=1e-12 * u.Unit("cm-2 s-1") * norm
    from gammapy.modeling.models import GaussianSpectralModel

    sigma = (10 ** (1 / bins_per_decade) - 1) * E
    return GaussianSpectralModel(mean=E, sigma=sigma, amplitude=amplitude)

print("Calculate energy range")

# selected_n_bins_per_decade = 20 # n bins per decade
max_rel_energy_error = 3

energy_axis = cube_map.geom.axes["energy"]
EminMap = energy_axis.edges[0]
EmaxMap = energy_axis.edges[-1]
stepE = energy_axis.edges[1] / energy_axis.edges[0]
nbins_per_decade = int(np.round(np.log(10) / np.log(stepE)))
Emin = EminMap / max_rel_energy_error
Emax = EmaxMap * max_rel_energy_error
nbins_per_decade, Emin, Emax

print("Create empty dataset")

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

print("Plotting GaussianSpectralModel")
from gammapy.modeling.models import GaussianSpectralModel

meanE = 1 * u.TeV
bins_per_decade = 20
sigma = (10 ** (1 / bins_per_decade) - 1) * meanE
amplitude = 1 * u.Unit("cm-2 s-1")
gm = GaussianSpectralModel(mean=meanE, sigma=sigma, amplitude=amplitude)
ax = gm.plot(energy_bounds=(0.1, 100) * u.TeV)
ax.set_yscale("linear")
gm.integral(meanE - 3 * sigma, meanE + 3 * sigma)

print("cube_map.get_spectrum()")
spec = cube_map.get_spectrum()
spec

# print('spec.plot()')
# spec.plot() # this plot shows dN/dE * E

# spec.data.shape, spec.data[spec.data.shape[0]//2,0,0]

print("Find norm bin")

energy_bins = cube_map.geom.axes["energy"].center
len(energy_bins), float(np.max(energy_bins) / u.TeV)
norm_bin = 0
for i, E in enumerate(energy_bins):
    if E > norm_energy * u.TeV:
        norm_bin = i
        break
assert norm_bin > 0
norm_bin

print("obtain norm_bin_width")
norm_bin_width = cube_map.geom.axes["energy"].bin_width[norm_bin]
norm_bin_width

print("find norm_flux")
# Npart=5000 # TODO update
# n_events_reduction_factor = 1 # suppress flux factor

int_bin_flux = (
    spec.data.flatten()
)  # we don't have to multiply by energy_bins /u.TeV since spectrum is already multiplied by E (see above)
norm_flux = int_bin_flux[norm_bin] / norm_bin_width
norm_flux
# int_bin_flux /= (Npart/200000 * np.max(int_bin_flux) * n_events_reduction_factor * 20/len(energy_bins)) # roughly 100 events
# int_bin_flux

print("find mult")
mult = norm_cm2_TeV_s * u.Unit("cm-2 s-1 TeV-1") / norm_flux  # .decompose()
mult

print("find int_bin_flux")
int_bin_flux = mult * int_bin_flux

int_bin_flux

print("Creating bin_models")

bin_models = []
for i, (flux, E) in enumerate(zip(int_bin_flux, energy_bins)):
    # print(i)
    if flux == 0:
        print("skipping bin ", i)
        continue
    # print(flux)
    spectral_model_delta = GetBinSpectralModel(
        E, amplitude=flux
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

print("Creating bkg_model")
bkg_model = FoVBackgroundModel(dataset_name="my-dataset")
models = Models(bin_models + [bkg_model])

print("dataset.models = models")
dataset.models = models

print("Creating sampler")
sampler = MapDatasetEventSampler(random_state=0)
print("Running sampler")
events = sampler.run(dataset, observation)

print(f"Source events: {(events.table['MC_ID'] > 0).sum()}")
print(f"Background events: {(events.table['MC_ID'] == 0).sum()}")

for i in range(1, len(bin_models) + 1):
    n = (events.table["MC_ID"] == i).sum()
    if n > 1:
        print(f"\tmodel {i}: {n} events")

print(f"Save events ...")
primary_hdu = fits.PrimaryHDU()
hdu_evt = fits.BinTableHDU(events.table)
hdu_gti = fits.BinTableHDU(dataset.gti.table, name="GTI")
hdu_all = fits.HDUList([primary_hdu, hdu_evt, hdu_gti])
hdu_all.writeto(f"./events.fits", overwrite=True)

print(f"Reading events ...")
hdul = fits.open("events.fits")
ev = hdul["EVENTS"].data
ra = ev["RA"]
dec = ev["DEC"]
en = ev["ENERGY"]

cube_map.geom

[cube_map.geom.center_coord[i] / cube_map.geom.data_shape[i] for i in (0, 1)]

ra_bins, dec_bins = (int(2 * x) for x in cube_map.geom.center_pix[:2])
ra_bins, dec_bins

cube_map.geom.center_skydir

Radius = float(min(cube_map.geom.width / 2 / u.degree).decompose())

print(f"Building event image ...")
plt.close()
from matplotlib.colors import LogNorm

cube_map.geom.width[0]

Nbins = 2 * int(Radius / pixsize) + 1
ra0 = np.mean(ra)
dec0 = np.mean(dec)
from numpy import cos, pi

cdec = cos(DEC_pnt * pi / 180.0)
ra_bins = np.linspace(
    RA_pnt - Radius / cdec, RA_pnt + Radius / cdec, Nbins + 1
)
dec_bins = np.linspace(DEC_pnt - Radius, DEC_pnt + Radius, Nbins + 1)

h = plt.hist2d(ra, dec, bins=[ra_bins, dec_bins], norm=LogNorm())
image = h[0]
plt.colorbar()
plt.xlabel("RA")
plt.ylabel("Dec")

print(f"Building event image 2 ...")
plt.close()
# Create a new WCS object.  The number of axes must be set
# from the start
w = wcs.WCS(naxis=2)

w.wcs.ctype = ["RA---CAR", "DEC--CAR"]
# we need a Plate carrée (CAR) projection since histogram is binned by ra-dec
# the peculiarity here is that CAR projection produces rectilinear grid only if CRVAL2==0
# also, we will follow convention of RA increasing from right to left (CDELT1<0, need to flip an input image)
# otherwise, aladin-lite doesn't show it
w.wcs.crval = [RA_pnt, 0]
w.wcs.crpix = [Nbins / 2.0 + 0.5, 1 - dec_bins[0] / pixsize]
w.wcs.cdelt = np.array([-pixsize / cdec, pixsize])

header = w.to_header()

hdu = fits.PrimaryHDU(np.flip(image.T, axis=1), header=header)
hdu.writeto("Image.fits", overwrite=True)
hdu = fits.open("Image.fits")
im = hdu[0].data
wcs1 = wcs.WCS(hdu[0].header)
ax = plt.subplot(projection=wcs1)
lon = ax.coords["ra"]
lon.set_major_formatter("d.dd")
lat = ax.coords["dec"]
lat.set_major_formatter("d.dd")
plt.imshow(im, origin="lower")
plt.colorbar(label="Counts")

plt.scatter(
    [RA_pnt],
    [DEC_pnt],
    marker="x",
    color="white",
    alpha=0.5,
    transform=ax.get_transform("world"),
)
plt.scatter(
    [RA],
    [DEC],
    marker="+",
    color="red",
    alpha=0.5,
    transform=ax.get_transform("world"),
)
plt.grid(color="white", ls="solid")
plt.xlabel("RA")
plt.ylabel("Dec")
plt.savefig("Image.png", format="png", bbox_inches="tight")

print("building event spectrum")
plt.close()
E = (events.energy / u.TeV).decompose()
ras = events.radec.ra.deg
decs = events.radec.dec.deg
# plt.hist(E,bins=np.logspace(-2,2,41))

mask = events.table["MC_ID"] > 0
plt.hist(E[mask], bins=np.logspace(-2, 2, 41), alpha=0.5, label="source")
mask = events.table["MC_ID"] == 0
plt.hist(E[mask], bins=np.logspace(-2, 2, 41), alpha=0.5, label="background")
plt.xlabel("E, TeV")

plt.xscale("log")
plt.yscale("log")
plt.legend(loc="upper right")
plt.savefig("event_spectrum.png", format="png")

print("reading events.fits")
hdul = fits.open("events.fits")
T_exp = hdul["EVENTS"].header["ONTIME"]
saved_events = hdul["EVENTS"].data

print("reading events.fits step 2")
coords_s = SkyCoord(RA, DEC, unit="degree")
RA_bkg = RA_pnt - (RA - RA_pnt)
DEC_bkg = DEC_pnt - (DEC - DEC_pnt)
coords_b = SkyCoord(RA_bkg, DEC_bkg, unit="degree")

Es = saved_events["ENERGY"]
ras = saved_events["RA"]
decs = saved_events["DEC"]
coords = SkyCoord(ras, decs, unit="degree")
seps_s = coords.separation(coords_s).deg
seps_b = coords.separation(coords_b).deg
coords_s = SkyCoord(RA, DEC, unit="degree")

print("manually parse irfs fits file")
hdul = fits.open(irfs_filename)

Aeff = hdul["EFFECTIVE AREA"].data
th_min = Aeff["THETA_LO"]
th_min
Emin_irf = Aeff["ENERG_LO"][0]
Emax_irf = Aeff["ENERG_HI"][0]
E_irf = sqrt(Emin_irf * Emax_irf)
Emin_irf
Ebins_irf = np.concatenate((Emin_irf, [Emax_irf[-1]]))
A = Aeff["EFFAREA"][0, 0]

print("Plot event histogram in energy")
plt.close()
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

print("Plot event histogram in zenith angle")
plt.close()
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

fits_events = BinaryProduct.from_file("events.fits")
fits_image = ImageDataProduct.from_fits_file("Image.fits")
bin_image = PictureProduct.from_file("Image.png")
spec_image = PictureProduct.from_file("event_spectrum.png")
theta2_png = PictureProduct.from_file("theta2.png")

spectrum_plot = spec_image  # http://odahub.io/ontology#ODAPictureProduct
theta_plot = theta2_png  # http://odahub.io/ontology#ODAPictureProduct
picture = bin_image  # http://odahub.io/ontology#ODAPictureProduct
image = fits_image  # http://odahub.io/ontology#Image
event_list = fits_events  # http://odahub.io/ontology#ODABinaryProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_model_CTA_events_from_file_spectrum_plot",
        "spectrum_plot_galaxy.output",
        spectrum_plot,
    )
)
_oda_outs.append(
    (
        "out_model_CTA_events_from_file_theta_plot",
        "theta_plot_galaxy.output",
        theta_plot,
    )
)
_oda_outs.append(
    (
        "out_model_CTA_events_from_file_picture",
        "picture_galaxy.output",
        picture,
    )
)
_oda_outs.append(
    ("out_model_CTA_events_from_file_image", "image_galaxy.output", image)
)
_oda_outs.append(
    (
        "out_model_CTA_events_from_file_event_list",
        "event_list_galaxy.output",
        event_list,
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
