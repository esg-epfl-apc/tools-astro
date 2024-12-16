#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import astropy.units as u

# conventional python routines
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import Angle  # Angles
from astropy.coordinates import SkyCoord  # High-level coordinates

# Conventional astronomical tools, also to be traced by Renku plugin, there is domain-specific ontology built in
from astropy.wcs import WCS

# from astroquery.desi import DESILegacySurvey
# from astroq.desi import DESILegacySurvey
from desi import DESILegacySurvey
from matplotlib.colors import LogNorm
from oda_api.json import CustomJSONEncoder

get_ipython().run_line_magic("matplotlib", "inline")   # noqa: F821
from matplotlib.colors import LogNorm

# if not(os.path.isdir('figs')):
#     os.makedirs('figs')
# if not(os.path.isdir('data')):
#     os.makedirs('data')

src_name = "Mrk 421"  # http://odahub.io/ontology#AstrophysicalObject
RA = 166.113808  # http://odahub.io/ontology#PointOfInterestRA
DEC = 38.208833  # http://odahub.io/ontology#PointOfInterestDEC
Radius = 3  # http://odahub.io/ontology#AngleMinutes
pixsize = (
    1.0  # http://odahub.io/ontology#AngleSeconds ; oda:label "Pixel size"
)
band = "g"  # http://odahub.io/ontology#String ; oda:allowed_value "g","r","i","z" ; oda:label "Band"
data_release = (
    9  # http://odahub.io/ontology#Integer ; oda:label "Data Release"
)

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
    "Radius",
    "pixsize",
    "band",
    "data_release",
]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

ra_s = RA
dec_s = DEC
image_size = Radius
image_band = band
dr = data_release
image_size = Angle(image_size * u.arcmin)
pixsize = Angle(pixsize * u.arcsec)
npix = int(2 * image_size / pixsize)
source = SkyCoord(ra_s, dec_s, unit="degree")

if dr < 10 and band == "i":
    raise RuntimeError(
        f"No data found. Data Release {dr} does not have '{band}' band."
    )

try:
    query = DESILegacySurvey.get_images(
        position=source,
        survey="dr%d" % dr,
        coordinates="icrs",
        data_release=dr,
        pixels=npix,
        radius=image_size,
        image_band=image_band,
    )
except:
    raise RuntimeError(
        f"No data found. Maybe (RA, Dec) = ({ra_s}, {dec_s}) is outside the covered region by Data Release {dr}."
    )

hdul = query[0]
hdul[0].header

# paramstring='ra='+str(ra_s)+'&dec='+str(dec_s)+'&size='+str(npix)+'&layer=ls-dr'+str(dr)+'&pixscale='+str(pixsize)+'&bands='+image_band
# suffix = hashlib.md5(paramstring.encode()).hexdigest()
# filename='data/image_legacysurvey_%s.fits'%( suffix )

# if os.path.exists(filename):
#         os.remove(filename)
# hdul.writeto(filename)

hdu = hdul[0]
image = hdu.data
w = WCS(hdu.header)
sky = w.pixel_to_world(0, 0)
ra_max_image = sky.ra.degree
dec_min_image = sky.dec.degree
sky = w.pixel_to_world(npix - 1, npix - 1)
ra_min_image = sky.ra.degree
dec_max_image = sky.dec.degree

plt.figure(figsize=(17, 13))
im = plt.imshow(
    image,
    norm=LogNorm(vmax=np.max(image), vmin=np.max(image) / 1e3),
    origin="lower",
    extent=(ra_max_image, ra_min_image, dec_min_image, dec_max_image),
)
plt.grid(color="black", ls="solid")
# plt.scatter(ra,dec,color='red',alpha=0.9)
plt.scatter([ra_s], [dec_s], color="blue", linewidth=4, alpha=0.3)

plt.xlabel("RA", fontsize=16)
plt.ylabel("DEC", fontsize=16)

plt.tick_params(axis="both", which="major", labelsize=16)
plt.colorbar(im)

plt.savefig("Image.png", format="png", bbox_inches="tight")

from oda_api.data_products import PictureProduct

bin_image = PictureProduct.from_file("Image.png")

picture = bin_image  # http://odahub.io/ontology#ODAPictureProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(("out_Image_picture", "picture_galaxy.output", picture))

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
