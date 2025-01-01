#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

workdir = os.getcwd()
repo_basedir = os.environ.get("BASEDIR", os.getcwd())
# global variables (DO NOT MODIFY)
npoints = 13
pathebl = (
    repo_basedir + "/dominguez_ebl_tau.txt"
)  # path with EBL model of Dominguez+11

import copy
import re

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.io import fits
from gammapy.data import Observation
from gammapy.datasets import MapDataset, MapDatasetEventSampler
from gammapy.irf import load_irf_dict_from_file
from gammapy.makers import MapDatasetMaker
from gammapy.maps import MapAxis, WcsGeom
from gammapy.modeling.models import FoVBackgroundModel, Models
from numpy import cos, exp, pi, sqrt
from oda_api.api import ProgressReporter
from oda_api.data_products import BinaryProduct, PictureProduct
from regions import CircleSkyRegion

# from gammapy.irf import load_cta_irfs

RA = 166.113809  # http://odahub.io/ontology#PointOfInterestRA
DEC = 38.208833  # http://odahub.io/ontology#PointOfInterestDEC
T1 = "2022-10-09T13:16:00.0"  # http://odahub.io/ontology#StartTime
T2 = "2022-10-10T13:16:00.0"  # http://odahub.io/ontology#EndTime

OffAxis_angle = 0.7  # http://odahub.io/ontology#AngleDegrees ; oda:label "Source off-axis angle"
Zd = 30.0  # http://odahub.io/ontology#Float ; oda:label "Zenith angle"

Texp = 1.0  # http://odahub.io/ontology#TimeIntervalHours
redshift = (
    0.13  # http://odahub.io/ontology#Float ; oda:label "Source redshift"
)
dN_dE = "2.0e-10*pow(E/1000., -1.99)*exp(-E/10000)"  # http://odahub.io/ontology#String ; oda:label "Source spectrum dN/dE 1/(TeV cm2 s1)]"

source_extension = 0.0  # http://odahub.io/ontology#Float ; oda:label "Source extension in degrees"

Radius_spectal_extraction = 0.2  # http://odahub.io/ontology#Float ; oda:group "Plotting" ; oda:label "Source radius aperture photometry"
Radius_sky_image = 4  # http://odahub.io/ontology#AngleDegrees ; oda:group "Plotting" ; oda:label "Image size"

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in [
    "RA",
    "DEC",
    "T1",
    "T2",
    "OffAxis_angle",
    "Zd",
    "Texp",
    "redshift",
    "dN_dE",
    "source_extension",
    "Radius_spectal_extraction",
    "Radius_sky_image",
]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

def parse_spectrum(input_str):
    dnde_str = copy.copy(input_str)

    numbers = re.findall(r"-?\d+\.?\d*(?:[eE][+-]?\d+)?", input_str)
    funcs = re.findall(r"\b\w+(?=\()", input_str)

    for num in numbers:
        input_str = input_str.replace(num, "")

    for func in funcs:
        input_str = input_str.replace(func, "")

    service_symbols = "E()[]*/.,+- "

    for symbol in service_symbols:
        input_str = input_str.replace(symbol, "")

    if input_str:
        raise ValueError("forbidden statements")

    return eval(f"lambda E: {dnde_str}")

# test_str = "2.0e-11*pow(E/1000., -1.99)*exp(-E/100); !wget https://scripts.com/myscript.sh"
Assumed = parse_spectrum(dN_dE)

z = redshift

Texp = Texp * 3600.0
DEC_pnt = DEC
cdec = cos(DEC * pi / 180.0)
RA_pnt1 = RA - OffAxis_angle / cdec
RA_pnt2 = RA + OffAxis_angle / cdec
Radius = Radius_sky_image
R_s = Radius_spectal_extraction

pointing1 = SkyCoord(RA_pnt1, DEC_pnt, unit="deg", frame="icrs")
pointing2 = SkyCoord(RA_pnt2, DEC_pnt, unit="deg", frame="icrs")
coord_s = SkyCoord(RA, DEC, unit="deg", frame="icrs")
RA_bkg1 = RA_pnt1 - (RA - RA_pnt1)
DEC_bkg1 = DEC_pnt - (DEC - DEC_pnt)
coord_b1 = SkyCoord(RA_bkg1, DEC_bkg1, unit="deg", frame="icrs")
RA_bkg2 = RA_pnt2 - (RA - RA_pnt2)
DEC_bkg2 = DEC_pnt - (DEC - DEC_pnt)
coord_b2 = SkyCoord(RA_bkg2, DEC_bkg2, unit="deg", frame="icrs")
offaxis = coord_s.separation(pointing1).deg
pr = ProgressReporter()
pr.report_progress(stage="Progress", progress=10.0)

Zeniths = np.array(
    [
        19.456,
        20.783,
        24.321,
        29.249,
        34.941,
        41.026,
        47.29,
        53.598,
        59.854,
        65.98,
    ]
)
Azimuths = np.array(
    [
        0.0,
        18.003,
        31.687,
        40.467,
        45.605,
        48.313,
        49.385,
        49.296,
        48.324,
        46.632,
    ]
)
len(Zeniths), len(Azimuths)

CTA_north_lat = 18.0
# Zd_min=abs(DEC-CTA_north_lat)
# if(Zd<Zd_min):
#    raise RuntimeError('Zenith angle is smaller than minimal possible for this source declination: '+str(Zd_min))
# if(Zd_min>70.):
#    raise RuntimeError('Source not visible from the LST site')

zen = Zeniths[Zd < Zeniths][0]
azi = Azimuths[Zd < Zeniths][0]
filename = (
    "IRFS_LST/irf_node_corsika_theta_"
    + str(zen)
    + "_az_"
    + str(azi)
    + "_.fits"
)
if os.path.exists(filename) == False:
    raise RuntimeError("No reponse function found")
    message = "No reponse function found!"
IRFS = load_irf_dict_from_file(filename)
hdul = fits.open(filename)
aeff = hdul["EFFECTIVE AREA"].data
ENERG_LO = aeff["ENERG_LO"][0]
ENERG_HI = aeff["ENERG_HI"][0]
THETA_LO = aeff["THETA_LO"][0]
THETA_HI = aeff["THETA_HI"][0]
EFFAREA = aeff["EFFAREA"][0]
print(offaxis)
ind_offaxis = len(THETA_LO[THETA_LO < offaxis] - 1)
EFAREA = EFFAREA[ind_offaxis]
HDU_EFFAREA = hdul["EFFECTIVE AREA"]
HDU_RMF = hdul["ENERGY DISPERSION"]

E = np.logspace(-2, 2, 20)

d = np.genfromtxt("Franceschini17.txt")
ee = d[:, 0]
z_grid = np.array([0.01, 0.03, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0])
ind = len(z_grid[z > z_grid]) - 1
coeff = (z - z_grid[ind]) / (z_grid[ind + 1] - z_grid[ind])
tau = d[:, ind + 1] + coeff * d[:, ind + 2]
tau_interp = np.interp(E, ee, tau)

def powerlaw_EBL():
    return Assumed(E * 1000) * exp(-tau_interp)

F = powerlaw_EBL()
plt.plot(E, E**2 * F)
plt.xscale("log")
plt.yscale("log")
plt.ylim(1e-14, 1e-10)

# filename = "IRFS/fits/Prod5-North-20deg-AverageAz-4LSTs09MSTs.180000s-v0.1.fits.gz"

if source_extension == 0.0:
    dic = {
        "components": [
            {
                "name": "Source1",
                "type": "SkyModel",
                "spectral": {
                    "type": "TemplateSpectralModel",
                    "parameters": [{"name": "norm", "value": 1.0}],
                    "energy": {"data": E.tolist(), "unit": "TeV"},
                    "values": {"data": F.tolist(), "unit": "1 / (cm2 TeV s)"},
                },
                "spatial": {
                    "type": "PointSpatialModel",
                    "frame": "icrs",
                    "parameters": [
                        {"name": "lon_0", "value": RA, "unit": "deg"},
                        {"name": "lat_0", "value": DEC, "unit": "deg"},
                    ],
                },
            },
            {
                "type": "FoVBackgroundModel",
                "datasets_names": ["my-dataset"],
                "spectral": {
                    "type": "PowerLawNormSpectralModel",
                    "parameters": [
                        {"name": "norm", "value": 1.0},
                        {"name": "tilt", "value": 0.0},
                        {"name": "reference", "value": 1.0, "unit": "TeV"},
                    ],
                },
            },
        ]
    }
else:
    dic = {
        "components": [
            {
                "name": "Source1",
                "type": "SkyModel",
                "spectral": {
                    "type": "TemplateSpectralModel",
                    "parameters": [{"name": "norm", "value": 1.0}],
                    "energy": {"data": E.tolist(), "unit": "TeV"},
                    "values": {"data": F.tolist(), "unit": "1 / (cm2 TeV s)"},
                },
                "spatial": {
                    "type": "GaussianSpatialModel",
                    "frame": "icrs",
                    "parameters": [
                        {"name": "lon_0", "value": RA, "unit": "deg"},
                        {"name": "lat_0", "value": DEC, "unit": "deg"},
                        {
                            "name": "sigma",
                            "value": source_extension,
                            "unit": "deg",
                        },
                    ],
                },
            },
            {
                "type": "FoVBackgroundModel",
                "datasets_names": ["my-dataset"],
                "spectral": {
                    "type": "PowerLawNormSpectralModel",
                    "parameters": [
                        {"name": "norm", "value": 1.0},
                        {"name": "tilt", "value": 0.0},
                        {"name": "reference", "value": 1.0, "unit": "TeV"},
                    ],
                },
            },
        ]
    }
modelsky = Models.from_dict(dic)

bkg_model = FoVBackgroundModel(dataset_name="my-dataset")

observation1 = Observation.create(
    obs_id="0", pointing=pointing1, livetime=str(Texp / 2.0) + " s", irfs=IRFS
)

print(f"Create the dataset")
energy_axis = MapAxis.from_energy_bounds(
    "0.012 TeV", "100 TeV", nbin=10, per_decade=True
)
energy_axis_true = MapAxis.from_energy_bounds(
    "0.001 TeV", "300 TeV", nbin=20, per_decade=True, name="energy_true"
)
migra_axis = MapAxis.from_bounds(
    0.5, 2, nbin=150, node_type="edges", name="migra"
)

geom1 = WcsGeom.create(
    skydir=pointing1,
    width=(2 * Radius, 2 * Radius),
    binsz=0.02,
    frame="icrs",
    axes=[energy_axis],
)

empty1 = MapDataset.create(
    geom1,
    energy_axis_true=energy_axis_true,
    migra_axis=migra_axis,
    name="my-dataset",
)
maker = MapDatasetMaker(selection=["exposure", "background", "psf", "edisp"])
dataset1 = maker.run(empty1, observation1)

region_sky1 = CircleSkyRegion(center=pointing1, radius=Radius * u.deg)
mask_map1 = dataset1.geoms["geom"].region_mask(region_sky1)
mod1 = modelsky.select_mask(mask_map1)

bkg_idx = np.where(np.array(modelsky.names) == "my-dataset-bkg")
mod1.append(modelsky[int(bkg_idx[0][0])])

dataset1.models = mod1

# for m in dataset.models[:-1]:
#    sep = m.spatial_model.position.separation(pointing).deg
#    print(f"This is the spatial separation of {m.name} from the pointing direction: {sep}")
pr.report_progress(stage="Progress", progress=50.0)
print("Simulate...")
sampler = MapDatasetEventSampler()
events1 = sampler.run(dataset1, observation1)

observation2 = Observation.create(
    obs_id="0", pointing=pointing2, livetime=str(Texp / 2.0) + " s", irfs=IRFS
)

geom2 = WcsGeom.create(
    skydir=pointing2,
    width=(2 * Radius, 2 * Radius),
    binsz=0.02,
    frame="icrs",
    axes=[energy_axis],
)

empty2 = MapDataset.create(
    geom2,
    energy_axis_true=energy_axis_true,
    migra_axis=migra_axis,
    name="my-dataset",
)
dataset2 = maker.run(empty2, observation2)

region_sky2 = CircleSkyRegion(center=pointing2, radius=Radius * u.deg)
mask_map2 = dataset2.geoms["geom"].region_mask(region_sky2)
mod2 = modelsky.select_mask(mask_map2)

bkg_idx = np.where(np.array(modelsky.names) == "my-dataset-bkg")
mod2.append(modelsky[int(bkg_idx[0][0])])

dataset2.models = mod2

# for m in dataset.models[:-1]:
#    sep = m.spatial_model.position.separation(pointing).deg
#    print(f"This is the spatial separation of {m.name} from the pointing direction: {sep}")
pr.report_progress(stage="Progress", progress=50.0)
print("Simulate...")
events2 = sampler.run(dataset2, observation2)

pr.report_progress(stage="Progress", progress=90.0)
print(f"Save events ...")
primary_hdu = fits.PrimaryHDU()
hdu_evt = fits.BinTableHDU(events1.table)
hdu_gti = fits.BinTableHDU(dataset1.gti.table, name="GTI")
hdu_all = fits.HDUList([primary_hdu, hdu_evt, hdu_gti, HDU_EFFAREA, HDU_RMF])
hdu_all.writeto(f"./events1.fits", overwrite=True)
primary_hdu = fits.PrimaryHDU()
hdu_evt = fits.BinTableHDU(events2.table)
hdu_gti = fits.BinTableHDU(dataset2.gti.table, name="GTI")
hdu_all = fits.HDUList([primary_hdu, hdu_evt, hdu_gti, HDU_EFFAREA, HDU_RMF])
hdu_all.writeto(f"./events2.fits", overwrite=True)

hdul = fits.open("events1.fits")
ev1 = hdul["EVENTS"].data
hdul = fits.open("events2.fits")
ev2 = hdul["EVENTS"].data

ra1 = ev1["RA"]
dec1 = ev1["DEC"]
en1 = ev1["ENERGY"]
coords1 = SkyCoord(ra1, dec1, unit="degree")
ra2 = ev2["RA"]
dec2 = ev2["DEC"]
en2 = ev2["ENERGY"]
coords2 = SkyCoord(ra2, dec2, unit="degree")

from matplotlib.colors import LogNorm

plt.figure()
pixsize = 0.1
Nbins = 2 * int(Radius / pixsize) + 1
ra0 = np.mean(ra)
dec0 = np.mean(dec)
from numpy import cos, pi

cdec = cos(DEC * pi / 180.0)
ra_bins = np.linspace(RA - Radius / cdec, RA + Radius / cdec, Nbins + 1)
dec_bins = np.linspace(DEC - Radius, DEC + Radius, Nbins + 1)

h1 = np.histogram2d(ra1, dec1, bins=[ra_bins, dec_bins])
h2 = np.histogram2d(ra2, dec2, bins=[ra_bins, dec_bins])
plt.imshow(np.transpose(h1[0] + h2[0]), norm=LogNorm())
image = h1[0] + h2[0]
plt.colorbar()
plt.xlabel("RA")
plt.ylabel("Dec")

ENERG_BINS = np.concatenate((ENERG_LO, [ENERG_HI[-1]]))
ENERG = sqrt(ENERG_LO * ENERG_HI)
m = coords1.separation(coord_s).deg < R_s
h1 = np.histogram(en1[m], bins=ENERG_BINS)
m = coords2.separation(coord_s).deg < R_s
h2 = np.histogram(en2[m], bins=ENERG_BINS)
cts_s = h1[0] + h2[0]
m = coords1.separation(coord_b1).deg < R_s
h1 = np.histogram(en1[m], bins=ENERG_BINS)
m = coords2.separation(coord_b2).deg < R_s
h2 = np.histogram(en2[m], bins=ENERG_BINS)
cts_b = h1[0] + h2[0]
src = cts_s - cts_b
src_err = sqrt(cts_s + cts_b)
plt.errorbar(
    ENERG,
    src,
    yerr=src_err,
    xerr=[ENERG - ENERG_LO, ENERG_HI - ENERG],
    linestyle="none",
    marker="o",
)

plt.xscale("log")
plt.xlabel(r"$E$, TeV")
plt.ylabel("Counts")
plt.yscale("log")
plt.ylim(0.1, 2 * max(src))
plt.savefig("Count_spectrum.png")

plt.figure()
th2_bins = np.linspace(0, 0.5, 50)
th2 = (th2_bins[:-1] + th2_bins[1:]) / 2.0
sep_s = coords1.separation(coord_s).deg

h1 = np.histogram(sep_s**2, bins=th2_bins)
sep_s = coords2.separation(coord_s).deg
h2 = np.histogram(sep_s**2, bins=th2_bins)
th2_src = h1[0] + h2[0]

sep_b = coords1.separation(coord_b1).deg
h1 = np.histogram(sep_b**2, bins=th2_bins)
sep_b = coords2.separation(coord_b2).deg
h2 = np.histogram(sep_b**2, bins=th2_bins)
th2_bkg = h1[0] + h2[0]
plt.step(th2, th2_src, where="mid")
plt.step(th2, th2_bkg, where="mid")
plt.axvline(R_s**2, color="black", linestyle="dashed")
plt.xlabel(r"$\theta^2$, degrees")
plt.ylabel("Counts")
plt.savefig("Theta2_plot.png")

# Create a new WCS object.  The number of axes must be set
# from the start
plt.figure()
w = wcs.WCS(naxis=2)

w.wcs.ctype = ["RA---CAR", "DEC--CAR"]
# we need a Plate carrée (CAR) projection since histogram is binned by ra-dec
# the peculiarity here is that CAR projection produces rectilinear grid only if CRVAL2==0
# also, we will follow convention of RA increasing from right to left (CDELT1<0, need to flip an input image)
# otherwise, aladin-lite doesn't show it
w.wcs.crval = [RA, 0]
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
    [RA_pnt1, RA_pnt2],
    [DEC_pnt, DEC_pnt],
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

fits_events = BinaryProduct.from_file("events.fits")
bin_image1 = PictureProduct.from_file("Image.png")
bin_image2 = PictureProduct.from_file("Theta2_plot.png")
bin_image3 = PictureProduct.from_file("Count_spectrum.png")
pr.report_progress(stage="Progress", progress=100.0)

image_png = bin_image1  # http://odahub.io/ontology#ODAPictureProduct
theta2_png = bin_image2  # http://odahub.io/ontology#ODAPictureProduct
spectrum_png = bin_image3  # http://odahub.io/ontology#ODAPictureProduct
event_list_fits = fits_events  # http://odahub.io/ontology#ODABinaryProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    ("out_LST1_events_image_png", "image_png_galaxy.output", image_png)
)
_oda_outs.append(
    ("out_LST1_events_theta2_png", "theta2_png_galaxy.output", theta2_png)
)
_oda_outs.append(
    (
        "out_LST1_events_spectrum_png",
        "spectrum_png_galaxy.output",
        spectrum_png,
    )
)
_oda_outs.append(
    (
        "out_LST1_events_event_list_fits",
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
