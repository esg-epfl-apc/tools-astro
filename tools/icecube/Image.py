#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import numpy as np
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.io import fits
from matplotlib import pyplot as plt
from numpy import cos, pi
from oda_api.api import ProgressReporter
from oda_api.data_products import ImageDataProduct, PictureProduct
from oda_api.json import CustomJSONEncoder

# src_name='NGC 1068' #http://odahub.io/ontology#AstrophysicalObject
RA = 40.669622  # http://odahub.io/ontology#PointOfInterestRA
DEC = -0.013294  # http://odahub.io/ontology#PointOfInterestDEC
# RA=308.65 # http://odahub.io/ontology#PointOfInterestRA
# DEC=40.9 # http://odahub.io/ontology#PointOfInterestDEC
# sigma=0.7  #http://odahub.io/ontology#AngleDegrees
Radius = 1.0  # http://odahub.io/ontology#AngleDegrees
pixel_size = 0.2  # http://odahub.io/ontology#AngleDegrees
TSmap_type = "Fixed_slope"  # http://odahub.io/ontology#String ; oda:allowed_value "Fixed_slope","Free_slope"
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
coords_s = SkyCoord(RA, DEC, unit="degree")
cdec = cos(DEC * pi / 180.0)
Npix = int(2 * (Radius / pixel_size)) + 1
print(Npix)
RA_grid = np.linspace(RA - Radius / cdec, RA + Radius / cdec, Npix)
DEC_grid = np.linspace(DEC - Radius, DEC + Radius, Npix)
TS_map = np.zeros((Npix, Npix))
ns_map = np.zeros((Npix, Npix))
gamma_map = np.zeros((Npix, Npix))

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

def process_pixel(index):
    i = int(index / Npix)
    j = index % Npix
    print(i, j)
    source = PointLikeSource(
        ra=np.deg2rad(RA_grid[i]), dec=np.deg2rad(DEC_grid[j])
    )
    ana = create_analysis(cfg=cfg, datasets=datasets, source=source)
    events_list = [data.exp for data in ana.data_list]
    ana.initialize_trial(events_list)
    (log_lambda_max, fitparam_values, status) = ana.llhratio.maximize(rss)
    (ts, x, status) = ana.unblind(rss)
    return RA_grid[i], DEC_grid[j], ts, x["ns"], x["gamma"]

def process_pixel1(index):
    i = int(index / Npix)
    j = index % Npix
    print(i, j)
    source = PointLikeSource(
        ra=np.deg2rad(RA_grid[i]), dec=np.deg2rad(DEC_grid[j])
    )
    ana = create_analysis(cfg=cfg, datasets=datasets, source=source)
    events_list = [data.exp for data in ana.data_list]
    ana.initialize_trial(events_list)
    TS_profile = []
    counts = np.linspace(0, 200, 200)
    for n in counts:
        TS_profile.append(2 * ana.llhratio.evaluate([n, 3.0])[0])
    return max(TS_profile), counts[np.argmax(TS_profile)]

# process_pixel(3)

pr = ProgressReporter()
pr.report_progress(stage="Progress", progress=5.0)

tsbest = 0
ibest = 0
jbest = 0
if TSmap_type == "Fixed_slope":
    for i, RRa in enumerate(RA_grid):
        for j, DDec in enumerate(DEC_grid):
            ind = i * Npix + j
            TS_map[i, j], ns_map[i, j] = process_pixel1(ind)
            if TS_map[i, j] > tsbest:
                tsbest = TS_map[i, j]
                ibest = i
                jbest = j
                print(RRa, DDec, tsbest)
else:
    ncpu = int(os.cpu_count())
    ncpu
    # with mp.Pool(3) as pool:
    #    res=pool.map(process_pixel, range(Npix**2))
    tsbest = 0
    for i, RRa in enumerate(RA_grid):
        pr.report_progress(stage="Progress", progress=(i / len(RA_grid)))
        for j, DDec in enumerate(DEC_grid):
            ind = i * Npix + j
            r, d, TS_map[i, j], ns_map[i, j], gamma_map[i, j] = process_pixel(
                ind
            )
            if TS_map[i, j] > tsbest:
                tsbest = TS_map[i, j]
                ibest = i
                jbest = j
                print(RRa, DDec, tsbest)

plt.imshow(
    np.flip(np.transpose(TS_map), axis=1),
    extent=(
        max(RA_grid) + pixel_size / cdec / 2.0,
        min(RA_grid) - pixel_size / cdec / 2.0,
        min(DEC_grid) - pixel_size / 2.0,
        max(DEC_grid) + pixel_size / 2.0,
    ),
    origin="lower",
    aspect=1 / cdec,
)
plt.colorbar(label="TS")
plt.scatter([RA_grid[ibest]], [DEC_grid[jbest]], marker="x", color="white")
plt.scatter([RA], [DEC], marker="+", color="red")
plt.xlabel("Right Ascension, degrees")
plt.ylabel("Declination, degrees")
plt.savefig("Image.png", format="png", bbox_inches="tight")
print(max(RA_grid), min(RA_grid), min(DEC_grid), max(DEC_grid))

# Create a new WCS object.  The number of axes must be set
# from the start
w = wcs.WCS(naxis=2)

w.wcs.ctype = ["RA---CAR", "DEC--CAR"]
# we need a Plate carrée (CAR) projection since histogram is binned by ra-dec
# the peculiarity here is that CAR projection produces rectilinear grid only if CRVAL2==0
# also, we will follow convention of RA increasing from right to left (CDELT1<0, need to flip an input image)
# otherwise, aladin-lite doesn't show it
w.wcs.crval = [RA, 0]
w.wcs.crpix = [Npix / 2.0 + 0.5, 1 - DEC_grid[0] / pixel_size]
w.wcs.cdelt = np.array([-pixel_size / cdec, pixel_size])

header = w.to_header()

hdu = fits.PrimaryHDU(np.flip(TS_map.T, axis=1), header=header)
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
plt.colorbar(label="TS")

plt.grid(color="white", ls="solid")
plt.xlabel("RA")
plt.ylabel("Dec")

fits_image = ImageDataProduct.from_fits_file("Image.fits")

bin_image = PictureProduct.from_file("Image.png")

picture = bin_image  # http://odahub.io/ontology#ODAPictureProduct
image = fits_image  # http://odahub.io/ontology#Image

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(("out_Image_picture", "picture_galaxy.output", picture))
_oda_outs.append(("out_Image_image", "image_galaxy.output", image))

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
