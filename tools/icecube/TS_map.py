#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import numpy as np
from astropy.coordinates import SkyCoord
from matplotlib import pyplot as plt
from numpy import cos, pi
from oda_api.data_products import PictureProduct
from oda_api.json import CustomJSONEncoder
from skyllh.analyses.i3.publicdata_ps.time_integrated_ps import create_analysis
from skyllh.core.config import Config
from skyllh.core.random import RandomStateService
from skyllh.core.source_model import PointLikeSource
from skyllh.datasets.i3.PublicData_10y_ps import create_dataset_collection

# src_name='NGC 1068' #http://odahub.io/ontology#AstrophysicalObject
# RA = 40.669622  # http://odahub.io/ontology#PointOfInterestRA
# DEC = -0.013294 # http://odahub.io/ontology#PointOfInterestDEC
RA = 308.2  # http://odahub.io/ontology#PointOfInterestRA
DEC = 41.0  # http://odahub.io/ontology#PointOfInterestDEC
sigma = 0.7  # http://odahub.io/ontology#AngleDegrees
Radius = 0.2  # http://odahub.io/ontology#AngleDegrees
pixel_size = 0.2  # http://odahub.io/ontology#AngleDegrees
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

import fileinput
import sys

def replaceAll(file, searchExp, replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp, replaceExp)
        sys.stdout.write(line)

get_ipython().system("cp signalpdf_template.py signalpdf.py")   # noqa: F821
replaceAll(
    "signalpdf.py",
    "sigma_sq = np.take(sigma**2, evt_idxs)",
    "sigma_sq = np.take(sigma**2, evt_idxs)+("
    + str(sigma)
    + "*np.pi/180.)**2",
)
get_ipython().system(   # noqa: F821
    "mv signalpdf.py /opt/conda/lib/python3.10/site-packages/skyllh/core"
)

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

# process_pixel(3)

ncpu = int(os.cpu_count())
ncpu
# with mp.Pool(3) as pool:
#    res=pool.map(process_pixel, range(Npix**2))
for i, RRa in enumerate(RA_grid):
    for j, DDec in enumerate(DEC_grid):
        ind = i * Npix + j
        print(ind, i, j)
        r, d, TS_map[i, j], ns_map[i, j], gamma_map[i, j] = process_pixel(ind)

plt.imshow(
    np.flip(TS_map, axis=1),
    extent=(
        max(RA_grid) + pixel_size / cdec / 2.0,
        min(RA_grid) - pixel_size / cdec / 2.0,
        min(DEC_grid) - pixel_size / 2.0,
        max(DEC_grid) + pixel_size / 2.0,
    ),
    origin="lower",
)
plt.colorbar(label="TS")
plt.xlabel("Right Ascension, degrees")
plt.ylabel("Declination, degrees")
plt.savefig("Image.png", format="png", bbox_inches="tight")
print(max(RA_grid), min(RA_grid), min(DEC_grid), max(DEC_grid))

bin_image = PictureProduct.from_file("Image.png")

picture = bin_image  # http://odahub.io/ontology#ODAPictureProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(("out_TS_map_picture", "picture_galaxy.output", picture))

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
