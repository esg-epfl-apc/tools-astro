#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

get_ipython().run_cell_magic(   # noqa: F821
    "bash",
    "",
    'rm -r IRFS | echo "Ok"\nmkdir IRFS\ncd IRFS\nwget https://zenodo.org/records/5499840/files/cta-prod5-zenodo-fitsonly-v0.1.zip\nunzip cta-prod5-zenodo-fitsonly-v0.1.zip\ncd fits\nfor fn in *.gz ; do tar -zxvf $fn; done \n',
)

import sys

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
from oda_api.api import ProgressReporter
from oda_api.data_products import BinaryProduct, PictureProduct

# not for run on Galaxy
# %%bash
# git lfs install
# git lfs pull

data_cube = "3d.fits"  # http://odahub.io/ontology#POSIXPath

# Source flux normalisaiton F0 in 1/(TeV cm2 s) at reference energy E0
F0 = 1e-11  # http://odahub.io/ontology#Float
E0 = 1.0  # http://odahub.io/ontology#Energy_TeV

OffAxis_angle = 0.78  # http://odahub.io/ontology#AngleDegrees

Radius_spectal_extraction = 0.2  # http://odahub.io/ontology#AngleDegrees
Radius_sky_image = 2.5  # http://odahub.io/ontology#AngleDegrees

Site = "North"  # http://odahub.io/ontology#String ; oda:allowed_value "North","South"
Telescopes_LST = True  # http://odahub.io/ontology#Boolean
Telescopes_MST = True  # http://odahub.io/ontology#Boolean
Telescopes_SST = False  # http://odahub.io/ontology#Boolean

Texp = 1.0  # http://odahub.io/ontology#TimeIntervalHours

z = 0.03  # http://odahub.io/ontology#Float

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

R_s = Radius_spectal_extraction
Radius = Radius_sky_image
LSTs = Telescopes_LST
MSTs = Telescopes_MST
SSTs = Telescopes_SST
file_path = data_cube

pr = ProgressReporter()
pr.report_progress(stage="Progress", progress=10.0)

print("loading " + file_path)
cube_map = Map.read(file_path)
cube_map.geom

print("locating source")
# source = SkyCoord.from_name(src_name, frame='icrs', parse=False, cache=True)
source = cube_map.geom.center_skydir
DEC = float(source.dec / u.deg)
RA = float(source.ra / u.deg)

CTA_south_lat = -25.0
CTA_north_lat = 18.0
if Site == "North":
    Zd = abs(DEC - CTA_north_lat)
    if Zd < 30.0:
        Zd = "20deg-"
    elif Zd < 50:
        Zd = "40deg-"
    elif Zd < 70.0:
        Zd = "60deg-"
    else:
        print("Source not visible from " + Site)
    if DEC > CTA_north_lat:
        N_S = "NorthAz-"
    else:
        N_S = "SouthAz-"
    if LSTs:
        tel = "4LSTs"
    if MSTs:
        tel += "09MSTs"
    filename = "IRFS/fits/Prod5-North-" + Zd + N_S + tel
else:
    Zd = abs(DEC - CTA_south_lat)
    if Zd < 30.0:
        Zd = "20deg-"
    elif Zd < 50:
        Zd = "40deg-"
    elif Zd < 70.0:
        Zd = "60deg-"
    else:
        print("Source not visible from " + Site)
    if DEC > CTA_south_lat:
        N_S = "NorthAz-"
    else:
        N_S = "SouthAz-"
    if MSTs:
        tel = "14MSTs"
    if SSTs:
        tel += "37MSTs"
    filename = "IRFS/fits/Prod5-South-" + Zd + N_S + tel

if Texp < 1800:
    filename += ".1800s-v0.1.fits.gz"
elif Texp < 18000:
    filename += ".18000s-v0.1.fits.gz"
else:
    filename += ".180000s-v0.1.fits.gz"
get_ipython().system("ls {filename}")   # noqa: F821

# telescope pointing will be shifted slightly
cdec = cos(DEC * pi / 180.0)
RA_pnt = RA - OffAxis_angle / cdec
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

pr.report_progress(stage="Progress", progress=20.0)

print("Find norm bin")

energy_bins = cube_map.geom.axes["energy"].center
len(energy_bins), float(np.max(energy_bins) / u.TeV)
norm_bin = 0
for i, E in enumerate(energy_bins):
    if E > E0 * u.TeV:
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
mult = F0 * u.Unit("cm-2 s-1 TeV-1") / norm_flux  # .decompose()
mult

print("find int_bin_flux")
int_bin_flux = mult * int_bin_flux

int_bin_flux

pr.report_progress(stage="Progress", progress=30.0)

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

hdul = fits.open(filename)
aeff = hdul["EFFECTIVE AREA"].data
ENERG_LO = aeff["ENERG_LO"][0]
ENERG_HI = aeff["ENERG_HI"][0]
THETA_LO = aeff["THETA_LO"][0]
THETA_HI = aeff["THETA_HI"][0]
EFFAREA = aeff["EFFAREA"][0]
ind_offaxis = len(THETA_LO[THETA_LO < OffAxis_angle] - 1)
EFAREA = EFFAREA[ind_offaxis]
HDU_EFFAREA = hdul["EFFECTIVE AREA"]
HDU_RMF = hdul["ENERGY DISPERSION"]

pr.report_progress(stage="Progress", progress=80.0)

print(f"Save events ...")
primary_hdu = fits.PrimaryHDU()
hdu_evt = fits.BinTableHDU(events.table)
hdu_gti = fits.BinTableHDU(dataset.gti.table, name="GTI")
hdu_all = fits.HDUList([primary_hdu, hdu_evt, hdu_gti, HDU_EFFAREA, HDU_RMF])
hdu_all.writeto(f"./events.fits", overwrite=True)

print(f"Reading events ...")
hdul = fits.open("events.fits")
ev = hdul["EVENTS"].data
ra = ev["RA"]
dec = ev["DEC"]
en = ev["ENERGY"]

[cube_map.geom.center_coord[i] / cube_map.geom.data_shape[i] for i in (0, 1)]

ra_bins, dec_bins = (int(2 * x) for x in cube_map.geom.center_pix[:2])
ra_bins, dec_bins

Radius = float(min(cube_map.geom.width / 2 / u.degree).decompose())

print(f"Building event image ...")
plt.close()
pixsize = 0.1
from matplotlib.colors import LogNorm

cube_map.geom.width[0]

Nbins = 2 * int(Radius / pixsize) + 1
ra0 = np.mean(ra)
dec0 = np.mean(dec)
from numpy import cos, pi

cdec = cos(DEC_pnt * pi / 180.0)
ra_bins = np.linspace(RA - Radius / cdec, RA + Radius / cdec, Nbins + 1)
dec_bins = np.linspace(DEC - Radius, DEC + Radius, Nbins + 1)

h = plt.hist2d(ra, dec, bins=[ra_bins, dec_bins], norm=LogNorm())
image = h[0]
plt.colorbar()
plt.xlabel("RA")
plt.ylabel("Dec")

print(f"Building event image 2 ...")
plt.figure()
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

coord_s = SkyCoord(RA, DEC, unit="degree")
RA_bkg = RA_pnt - (RA - RA_pnt)
DEC_bkg = DEC_pnt - (DEC - DEC_pnt)
coord_b = SkyCoord(RA_bkg, DEC_bkg, unit="degree")
coords = SkyCoord(ra, dec, unit="degree")

plt.figure()
ev_src = en[coords.separation(coord_s).deg < R_s]
ev_bkg = en[coords.separation(coord_b).deg < R_s]
ENERG_BINS = np.concatenate((ENERG_LO, [ENERG_HI[-1]]))
ENERG = sqrt(ENERG_LO * ENERG_HI)
h1 = np.histogram(ev_src, bins=ENERG_BINS)
h2 = np.histogram(ev_bkg, bins=ENERG_BINS)
cts_s = h1[0]
cts_b = h2[0]
src = cts_s - cts_b
src_err = sqrt(cts_s + cts_b)
plt.errorbar(ENERG, src, src_err)
plt.axhline(0, linestyle="dashed", color="black")
plt.xscale("log")
plt.xlabel(r"$E$, TeV")
plt.ylabel("Counts")
plt.yscale("log")
plt.ylim(0.1, 2 * max(src))
plt.savefig("Count_spectrum.png")

plt.figure()
sep_s = coords.separation(coord_s).deg
sep_b = coords.separation(coord_b).deg
plt.hist(sep_s**2, bins=np.linspace(0, 0.5, 50))
plt.hist(sep_b**2, bins=np.linspace(0, 0.5, 50))
plt.axvline(R_s**2, color="black", linestyle="dashed")
plt.xlabel(r"$\theta^2$, degrees")
plt.ylabel("Counts")
plt.savefig("Theta2_plot.png")

pr.report_progress(stage="Progress", progress=100.0)

fits_events = BinaryProduct.from_file("events.fits")
bin_image = PictureProduct.from_file("Image.png")
spec_image = PictureProduct.from_file("Count_spectrum.png")
theta2_png = PictureProduct.from_file("Theta2_plot.png")

spectrum_png = spec_image  # http://odahub.io/ontology#ODAPictureProduct
theta2_png = theta2_png  # http://odahub.io/ontology#ODAPictureProduct
image_png = bin_image  # http://odahub.io/ontology#ODAPictureProduct
event_list_fits = fits_events  # http://odahub.io/ontology#ODABinaryProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_model_cube_file_spectrum_png",
        "spectrum_png_galaxy.output",
        spectrum_png,
    )
)
_oda_outs.append(
    ("out_model_cube_file_theta2_png", "theta2_png_galaxy.output", theta2_png)
)
_oda_outs.append(
    ("out_model_cube_file_image_png", "image_png_galaxy.output", image_png)
)
_oda_outs.append(
    (
        "out_model_cube_file_event_list_fits",
        "event_list_fits_galaxy.output",
        event_list_fits,
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
