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

if os.path.exists("hess_dl3_dr1.tar.gz") == False:
    get_ipython().system(   # noqa: F821
        "wget https://zenodo.org/record/1421099/files/hess_dl3_dr1.tar.gz"
    )
    get_ipython().system("tar -zxvf hess_dl3_dr1.tar.gz")   # noqa: F821

src_name = "Crab"  # http://odahub.io/ontology#AstrophysicalObject
RA = 83.628700  # http://odahub.io/ontology#PointOfInterestRA
DEC = 22.014700  # http://odahub.io/ontology#PointOfInterestDEC
T1 = "2000-10-09T13:16:00.0"  # http://odahub.io/ontology#StartTime
T2 = "2022-10-10T13:16:00.0"  # http://odahub.io/ontology#EndTime
Radius = 1.0  # http://odahub.io/ontology#AngleDegrees
pixsize = (
    0.05  # http://odahub.io/ontology#AngleDegrees ; oda:label "Pixel size"
)
Emin = 100.0  # http://odahub.io/ontology#Energy_GeV
Emax = 10000.0  # http://odahub.io/ontology#Energy_GeV

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
message = ""
RA_pnts = []
DEC_pnts = []
DL3_files = []
OBSIDs = []
Tstart = []
Tstop = []
flist = os.listdir("data")
for f in flist:
    if f[-7:] == "fits.gz":
        DL3_files.append(f)
        OBSIDs.append(int(f[20:26]))
        hdul = fits.open("data/" + f)
        RA_pnts.append(float(hdul[1].header["RA_PNT"]))
        DEC_pnts.append(float(hdul[1].header["DEC_PNT"]))
        Tstart.append(
            Time(
                hdul[1].header["DATE-OBS"] + "T" + hdul[1].header["TIME-OBS"],
                format="isot",
                scale="utc",
            ).mjd
        )
        Tstop.append(
            Time(
                hdul[1].header["DATE-END"] + "T" + hdul[1].header["TIME-END"],
                format="isot",
                scale="utc",
            ).mjd
        )
        hdul.close()

Coords_s = SkyCoord(RA, DEC, unit="degree")
COORDS_pnts = SkyCoord(RA_pnts, DEC_pnts, unit="degree")
seps = COORDS_pnts.separation(Coords_s).deg

mask = np.where((seps < Radius) & (Tstart > T1) & (Tstop < T2))[0]
OBSlist = []
for i in mask:
    OBSlist.append(DL3_files[i])
if len(OBSlist) == 0:
    message = "No data found"
    raise RuntimeError("No data found")
message

cdec = cos(DEC * pi / 180.0)
Npix = int(2 * Radius / pixsize) + 1
RA_bins = np.linspace(
    RA - Npix * pixsize / cdec / 2, RA + Npix * pixsize / cdec / 2, Npix + 1
)
DEC_bins = np.linspace(
    DEC - Npix * pixsize / 2, DEC + Npix * pixsize / 2, Npix + 1
)

image = np.zeros((Npix, Npix))
for f in OBSlist:
    hdul = fits.open("data/" + f)
    ev = hdul["EVENTS"].data
    ev_ra = ev["RA"]
    ev_dec = ev["DEC"]
    ev_en = ev["ENERGY"]
    ev_time = ev["TIME"]
    h = np.histogram2d(ev_ra, ev_dec, bins=[RA_bins, DEC_bins])
    image += h[0]
    hdul.close()

image = np.transpose(image)

plt.imshow(
    image,
    extent=(RA_bins[0], RA_bins[-1], DEC_bins[0], DEC_bins[-1]),
    origin="lower",
)
plt.colorbar()

plt.xlim(*plt.xlim()[::-1])

plt.xlabel("RA, degrees")
plt.ylabel("DEC,degrees")
plt.savefig("Image.png", format="png")

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

bin_image = PictureProduct.from_file("Image.png")
fits_image = ImageDataProduct.from_fits_file("Image.fits")

png = bin_image  # http://odahub.io/ontology#ODAPictureProduct
fits = fits_image  # http://odahub.io/ontology#Image

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(("out_Image_png", "png_galaxy.output", png))
_oda_outs.append(("out_Image_fits", "fits_galaxy.output", fits))

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
