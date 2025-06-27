#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import Angle  # Angles
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astroquery.skyview import SkyView
from matplotlib.colors import LogNorm
from oda_api.json import CustomJSONEncoder

src_name = "Mrk 421"  # http://odahub.io/ontology#AstrophysicalObject
RA = 166.113808  # http://odahub.io/ontology#PointOfInterestRA
DEC = 38.208833  # http://odahub.io/ontology#PointOfInterestDEC
Radius = 3  # http://odahub.io/ontology#AngleMinutes ; oda:label "Image size"
pixsize = (
    1.0  # http://odahub.io/ontology#AngleSeconds ; oda:label "Pixel size"
)
band = "3.4"  # http://odahub.io/ontology#String ; oda:allowed_value "3.4","4.6","12","22" ; oda:label "Band"

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in ["src_name", "RA", "DEC", "Radius", "pixsize", "band"]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

image_band = "WISE " + band
image_size = Angle(Radius * u.arcmin)
pixsize = Angle(pixsize * u.arcsec)
npix = int(2 * image_size / pixsize)
source = SkyCoord(RA, DEC, unit="degree")
npix

paths = SkyView.get_images(
    position=source, coordinates="icrs", survey=[image_band]
)

hdu = paths[0][0]
image = hdu.data

hdu.writeto("Image.fits", overwrite=True)

w = WCS(hdu.header)
sky = w.pixel_to_world(0, 0)
ra_max_image = sky.ra.degree
dec_min_image = sky.dec.degree
npix = image.shape
sky = w.pixel_to_world(npix[0] - 1, npix[1] - 1)
ra_min_image = sky.ra.degree
dec_max_image = sky.dec.degree

npix = image.shape
npix[0]

plt.figure(figsize=(17, 13))
im = plt.imshow(
    image,
    norm=LogNorm(vmax=np.max(image), vmin=np.max(image) / 1e3),
    origin="lower",
    extent=(ra_max_image, ra_min_image, dec_min_image, dec_max_image),
)
plt.grid(color="black", ls="solid")
# plt.scatter(ra,dec,color='red',alpha=0.9)
plt.scatter([RA], [DEC], linewidth=4, color="white", marker="x")

plt.xlabel("RA", fontsize=16)
plt.ylabel("DEC", fontsize=16)

plt.tick_params(axis="both", which="major", labelsize=16)
plt.colorbar(im)

plt.savefig("Image.png", format="png", bbox_inches="tight")

from oda_api.data_products import ImageDataProduct, PictureProduct

bin_image = PictureProduct.from_file("Image.png")
fits_image = ImageDataProduct.from_fits_file("Image.fits")

picture = bin_image  # http://odahub.io/ontology#ODAPictureProduct
fits = fits_image  # http://odahub.io/ontology#Image

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(("out_Image_picture", "picture_galaxy.output", picture))
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
