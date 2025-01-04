#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from numpy import cos, pi
from oda_api.data_products import ImageDataProduct, PictureProduct
from oda_api.json import CustomJSONEncoder

src_name = "Crab"  # http://odahub.io/ontology#AstrophysicalObject
RA = 83.628700  # http://odahub.io/ontology#PointOfInterestRA
DEC = 22.014700  # http://odahub.io/ontology#PointOfInterestDEC

T1 = "2000-10-09T13:16:00.0"  # http://odahub.io/ontology#StartTime
T2 = "2024-10-10T13:16:00.0"  # http://odahub.io/ontology#EndTime
Radius_search = 2.0  # http://odahub.io/ontology#AngleDegrees ; oda:label "Cone search radius"
Radius_image = 2.0  # http://odahub.io/ontology#AngleDegrees ; oda:label "Image radius" ; oda:group "Plotting"
pixsize = 0.025  # http://odahub.io/ontology#AngleDegrees ; oda:label "Pixel size" ; oda:group "Plotting"
Emin = 0.1  # http://odahub.io/ontology#Energy_TeV ; oda:label "Minimal energy" ; oda:group "Plotting"
Emax = 20  # http://odahub.io/ontology#Energy_TeV ; oda:label "Maximal energy" ; oda:group "Plotting"

Offset = 0.4  # http://odahub.io/ontology#AngleDegrees ; oda:label "Source off-axis angle"

NSB = 0  # http://odahub.io/ontology#Integer ; oda:label "Night sky background level (0-0.8)" ; allowed_value 0,1,2,3,4,5,6,7,8

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in [
    "src_name",
    "RA",
    "DEC",
    "T1",
    "T2",
    "Radius_search",
    "Radius_image",
    "pixsize",
    "Emin",
    "Emax",
    "Offset",
    "NSB",
]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

racrab = 83.628700
deccrab = 22.014700
crab = SkyCoord(racrab, deccrab, unit="degree")
source = SkyCoord(RA, DEC, unit="degree")
sep = source.separation(crab).deg
if sep > 2:
    raise ValueError("Public data release is limited to pointings around Crab")

T1 = Time(T1, format="isot", scale="utc").mjd
T2 = Time(T2, format="isot", scale="utc").mjd

workdir = os.getcwd()
repo_basedir = os.environ.get("BASEDIR", os.getcwd())
data_dir = repo_basedir + "/magic_dl3_pdr1-main/data/CrabNebula"
get_ipython().system("ls {data_dir}")   # noqa: F821

if NSB == 0:
    data_dir += "/dark"
    if Offset == 0.4:
        data_dir += "/single_offset"
    else:
        data_dir += "/multi_offset"
        if Offset == 0.2:
            data_dir += "/offset_0.20"
        elif Offset == 0.35:
            data_dir += "/offset_0.35"
        elif Offset == 0.7:
            data_dir += "/offset_0.70"
        elif Offset == 1.0:
            data_dir += "/offset_1.00"
        elif Offset == 1.4:
            data_dir += "/offset_1.40"
        else:
            raise ValueError("Offset angle value not found")
else:
    data_dir += "/moon"
    if NSB < 3:
        data_dir += "/NSB_1-2"
    elif NSB < 4:
        data_dir += "/NSB_2-3"
    elif NSB < 6:
        data_dir += "/NSB_3-5"
    elif NSB < 9:
        data_dir += "/NSB_5-8"
    else:
        raise ValueError("NSB level not found ")
get_ipython().system("ls {data_dir}")   # noqa: F821

from pathlib import Path

from gammapy.data import DataStore

path = Path(data_dir)
# data_store=DataStore.from_dir(path)
paths = list(path.rglob("*DL3*.fits"))
data_store = DataStore.from_events_files(paths)

selection = dict(
    type="sky_circle",
    frame="icrs",
    lon=str(RA) + " deg",
    lat=str(DEC) + " deg",
    radius=str(Radius_search) + " deg",
)
selected_obs_table = data_store.obs_table.select_observations(selection)
selected_obs_table

RA_pnts = selected_obs_table["RA_PNT"]
DEC_pnts = selected_obs_table["DEC_PNT"]
Tstart = selected_obs_table["TSTART"]
Tstop = selected_obs_table["TSTOP"]
DL3_fname = selected_obs_table["EVENTS_FILENAME"]

cdec = cos(DEC * pi / 180.0)
Npix = int(2 * Radius_image / pixsize) + 1
RA_bins = np.linspace(
    RA - Npix * pixsize / cdec / 2, RA + Npix * pixsize / cdec / 2, Npix + 1
)
DEC_bins = np.linspace(
    DEC - Npix * pixsize / 2, DEC + Npix * pixsize / 2, Npix + 1
)

image = np.zeros((Npix, Npix))
image_b = np.zeros((Npix, Npix))
for f in DL3_fname:
    print(f)
    hdul = fits.open(f)
    RA_pnt = hdul[1].header["RA_PNT"]
    DEC_pnt = hdul[1].header["DEC_PNT"]
    ev = hdul["EVENTS"].data
    ev_ra = ev["RA"]
    ev_dec = ev["DEC"]
    dRA = ev_ra - RA_pnt
    dDEC = ev_dec - DEC_pnt
    bkg_ra = RA_pnt - dRA
    bkg_dec = DEC_pnt - dDEC

    ev_en = ev["ENERGY"]
    ev_time = ev["TIME"]
    mask = (ev_en > Emin) & (ev_en < Emax)
    h = np.histogram2d(ev_ra[mask], ev_dec[mask], bins=[RA_bins, DEC_bins])
    image += h[0]
    h = np.histogram2d(bkg_ra[mask], bkg_dec[mask], bins=[RA_bins, DEC_bins])
    image_b += h[0]
    hdul.close()

image_s = image - image_b
image = np.transpose(image)
image_s = np.transpose(image_s)

plt.imshow(
    image,
    extent=(RA_bins[0], RA_bins[-1], DEC_bins[0], DEC_bins[-1]),
    origin="lower",
)
plt.colorbar()

# Create a new WCS object.  The number of axes must be set
# from the start
w = wcs.WCS(naxis=2)

w.wcs.ctype = ["RA---CAR", "DEC--CAR"]
# we need a Plate carrée (CAR) projection since histogram is binned by ra-dec
# the peculiarity here is that CAR projection produces rectilinear grid only if CRVAL2==0
# also, we will follow convention of RA increasing from right to left (CDELT1<0, need to flip an input image)
# otherwise, aladin-lite doesn't show it
w.wcs.crval = [RA, 0]
w.wcs.crpix = [Npix / 2.0 + 0.5, 0.5 - DEC_bins[0] / pixsize]
w.wcs.cdelt = np.array([-pixsize / cdec, pixsize])

header = w.to_header()

hdu = fits.PrimaryHDU(np.flip(image, axis=1), header=header)
hdu.writeto("Image.fits", overwrite=True)
hdu = fits.open("Image.fits")
im = hdu[0].data
wcs1 = wcs.WCS(hdu[0].header)
ax = plt.subplot(projection=wcs1)
plt.imshow(im, origin="lower")
plt.colorbar(label="Counts per pixel")
plt.scatter(
    [RA], [DEC], marker="x", color="white", transform=ax.get_transform("world")
)
plt.text(
    RA,
    DEC + 0.5 * pixsize,
    src_name,
    color="white",
    transform=ax.get_transform("world"),
)

plt.grid(color="white", ls="solid")
plt.xlabel("RA")
plt.ylabel("Dec")
plt.savefig("Image.png", format="png")

fits_image = ImageDataProduct.from_fits_file("Image.fits")
bin_image = PictureProduct.from_file("Image.png")

png = bin_image  # http://odahub.io/ontology#ODAPictureProduct
fits = fits_image  # http://odahub.io/ontology#Image

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(("out_Image_public_dl3_png", "png_galaxy.output", png))
_oda_outs.append(("out_Image_public_dl3_fits", "fits_galaxy.output", fits))

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
