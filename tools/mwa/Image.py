#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import astropy.units as u
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astroquery.skyview import SkyView
from oda_api.data_products import ImageDataProduct, PictureProduct
from oda_api.json import CustomJSONEncoder

class AnalysisError(RuntimeError): ...

src_name = "1ES 0229+200"  # http://odahub.io/ontology#AstrophysicalObject
RA = 38.202562  # http://odahub.io/ontology#PointOfInterestRA
DEC = 20.288191  # http://odahub.io/ontology#PointOfInterestDEC
T1 = "2000-10-09T13:16:00.0"  # http://odahub.io/ontology#StartTime
T2 = "2022-10-10T13:16:00.0"  # http://odahub.io/ontology#EndTime
Radius = 0.1  # http://odahub.io/ontology#AngleDegrees
pixsize = 0.01  # http://odahub.io/ontology#AngleDegrees
Frequency = "GLEAM 170-231 MHz"  # http://odahub.io/ontology#String ; oda:allowed_value "GLEAM 72-103 MHz","GLEAM 103-134 MHz","GLEAM 139-170 MHz","GLEAM 170-231 MHz"

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
    "Radius",
    "pixsize",
    "Frequency",
]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

pixels = int(2 * Radius / pixsize) + 1
Radius *= u.deg
pos = str(RA) + ", " + str(DEC)
pixels

try:
    hdul = SkyView.get_images(
        position=pos, survey=[Frequency], pixels=pixels, radius=Radius
    )
except:
    raise AnalysisError("No data found")
    message = "No data found!"

hdu = hdul[0]
hdu[0].header
wcs = WCS(hdu[0].header)

image = hdu[0].data

ax = plt.subplot(projection=wcs)
im = ax.imshow(image, origin="lower")
ax.coords.grid(True, color="white", ls="solid")
plt.colorbar(im, label="Jy/beam")
plt.savefig("Image.png", format="png", bbox_inches="tight")

hdu.writeto("Image.fits", overwrite=True)
bin_image = PictureProduct.from_file("Image.png")
fits_image = ImageDataProduct.from_fits_file("Image.fits")

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
