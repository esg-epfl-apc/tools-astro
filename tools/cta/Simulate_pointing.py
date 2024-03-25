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
from gammapy.modeling.models import Models
from oda_api.api import ProgressReporter
from oda_api.data_products import (
    BinaryProduct,
    ImageDataProduct,
    PictureProduct,
)
from regions import CircleSkyRegion

# from gammapy.irf import load_cta_irfs

RA = 166.113809  # http://odahub.io/ontology#PointOfInterestRA
DEC = 38.208833  # http://odahub.io/ontology#PointOfInterestDEC
Radius = 2.5  # http://odahub.io/ontology#AngleDegrees
# Exposure time in hours
Texp = 1.0  # http://odahub.io/ontology#TimeIntervalHours
# Source redshift
z = 0.1  # http://odahub.io/ontology#Float
# Source flux normalisaiton F0 in 1/(TeV cm2 s) at reference energy E0
F0 = 4e-13  # http://odahub.io/ontology#Float
E0 = 1.0  # http://odahub.io/ontology#Energy_TeV
Gamma = 1.75  # http://odahub.io/ontology#Float
# source extension in degrees
sigma = 0.0  # http://odahub.io/ontology#Float
RA_pnt = 167.113809  # http://odahub.io/ontology#RightAscensionDegrees
DEC_pnt = 38.208833  # http://odahub.io/ontology#DeclinationDegrees

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

Texp = Texp * 3600.0
pointing = SkyCoord(RA_pnt, DEC_pnt, unit="deg", frame="icrs")
pr = ProgressReporter()
pr.report_progress(stage="Progress", progress=10.0)

get_ipython().system(   # noqa: F821
    "ls IRFS/fits/Prod5-North-20deg-AverageAz-4LSTs09MSTs.180000s-v0.1.fits.gz"
)
filename = (
    "IRFS/fits/Prod5-North-20deg-AverageAz-4LSTs09MSTs.180000s-v0.1.fits.gz"
)
IRFS = load_irf_dict_from_file(filename)
dic = {
    "components": [
        {
            "name": "Source1",
            "type": "SkyModel",
            "spectral": {
                "type": "PowerLawSpectralModel",
                "parameters": [
                    {"name": "index", "value": Gamma},
                    {
                        "name": "amplitude",
                        "value": F0,
                        "unit": "TeV-1 s-1 cm-2",
                    },
                    {"name": "reference", "value": E0, "unit": "TeV"},
                ],
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
modelsky = Models.from_dict(dic)

observation = Observation.create(
    obs_id="0", pointing=pointing, livetime=str(Texp) + " s", irfs=IRFS
)

print(f"Create the dataset for {pointing}")
energy_axis = MapAxis.from_energy_bounds(
    "0.012 TeV", "100 TeV", nbin=10, per_decade=True
)
energy_axis_true = MapAxis.from_energy_bounds(
    "0.001 TeV", "300 TeV", nbin=20, per_decade=True, name="energy_true"
)
migra_axis = MapAxis.from_bounds(
    0.5, 2, nbin=150, node_type="edges", name="migra"
)

geom = WcsGeom.create(
    skydir=pointing,
    width=(2 * Radius, 2 * Radius),
    binsz=0.02,
    frame="icrs",
    axes=[energy_axis],
)

empty = MapDataset.create(
    geom,
    energy_axis_true=energy_axis_true,
    migra_axis=migra_axis,
    name="my-dataset",
)
maker = MapDatasetMaker(selection=["exposure", "background", "psf", "edisp"])
dataset = maker.run(empty, observation)

region_sky = CircleSkyRegion(center=pointing, radius=Radius * u.deg)
mask_map = dataset.geoms["geom"].region_mask(region_sky)
mod = modelsky.select_mask(mask_map)

bkg_idx = np.where(np.array(modelsky.names) == "my-dataset-bkg")
mod.append(modelsky[int(bkg_idx[0][0])])

dataset.models = mod

for m in dataset.models[:-1]:
    sep = m.spatial_model.position.separation(pointing).deg
    print(
        f"This is the spatial separation of {m.name} from the pointing direction: {sep}"
    )
pr.report_progress(stage="Progress", progress=50.0)
print("Simulate...")
sampler = MapDatasetEventSampler()
events = sampler.run(dataset, observation)

pr.report_progress(stage="Progress", progress=90.0)
print(f"Save events ...")
primary_hdu = fits.PrimaryHDU()
hdu_evt = fits.BinTableHDU(events.table)
hdu_gti = fits.BinTableHDU(dataset.gti.table, name="GTI")
hdu_all = fits.HDUList([primary_hdu, hdu_evt, hdu_gti])
hdu_all.writeto(f"./events.fits", overwrite=True)

hdul = fits.open("events.fits")
ev = hdul["EVENTS"].data
ra = ev["RA"]
dec = ev["DEC"]
en = ev["ENERGY"]

from matplotlib.colors import LogNorm

pixsize = 0.1
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

fits_events = BinaryProduct.from_file("events.fits")
fits_image = ImageDataProduct.from_fits_file("Image.fits")
bin_image = PictureProduct.from_file("Image.png")
pr.report_progress(stage="Progress", progress=100.0)

picture = bin_image  # http://odahub.io/ontology#ODAPictureProduct
image = fits_image  # http://odahub.io/ontology#Image
event_list = fits_events  # http://odahub.io/ontology#ODABinaryProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    ("out_Simulate_pointing_picture", "picture_galaxy.output", picture)
)
_oda_outs.append(("out_Simulate_pointing_image", "image_galaxy.output", image))
_oda_outs.append(
    (
        "out_Simulate_pointing_event_list",
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
