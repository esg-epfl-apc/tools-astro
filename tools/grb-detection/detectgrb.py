#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python

# This script is generated with nb2galaxy

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

T1 = "2023-01-16T04:53:33.9"  # http://odahub.io/ontology#StartTime
T2 = "2023-01-16T04:55:33.9"  # http://odahub.io/ontology#EndTime
detection_time_scales = "1,10"
lc_time_scale = 0.1  # https://odahub.io/ontology#TimeIntervalSeconds
background_age = 10  # oda:TimeIntervalSeconds
min_sn = 5  # https://odahub.io/ontology#SignalToNoiseRatio

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "C_data_product_" in inp_dic.keys():
    inp_pdic = inp_dic["C_data_product_"]
else:
    inp_pdic = inp_dic
T1 = str(inp_pdic["T1"])
T2 = str(inp_pdic["T2"])
detection_time_scales = str(inp_pdic["detection_time_scales"])
lc_time_scale = float(inp_pdic["lc_time_scale"])
background_age = float(inp_pdic["background_age"])
min_sn = int(inp_pdic["min_sn"])

import numpy as np
from astropy.time import Time
from matplotlib import pylab as plt
from oda_api.api import DispatcherAPI
from oda_api.data_products import LightCurveDataProduct, PictureProduct

disp = DispatcherAPI(
    url="https://www.astro.unige.ch/mmoda//dispatch-data", instrument="mock"
)

par_dict = {
    "T1": T1,
    "T2": T2,
    "T_format": "isot",
    "instrument": "spi_acs",
    "product": "spi_acs_lc",
    "product_type": "Real",
    "time_bin": lc_time_scale,
    "time_bin_format": "sec",
}

data_collection = disp.get_product(**par_dict)

lc = data_collection.spi_acs_lc_0_query.data_unit[1].data

lc["TIME"] = (
    lc["TIME"] / 24 / 3600 + (Time(T1).mjd + Time(T2).mjd) / 2.0
)  # TODO: more accurately

def rebin(x, n):
    N = int(len(x) / n)
    return np.array(x[: N * n]).reshape((N, n)).sum(1)

lc = LightCurveDataProduct.from_arrays(
    times=Time(lc["TIME"], format="mjd"), fluxes=lc["RATE"], errors=lc["ERROR"]
)

obj_results = [lc.encode()]

T0_ijd = lc.data_unit[1].data["TIME"][np.argmax(lc.data_unit[1].data["FLUX"])]

ijd2plot = lambda x: (x - T0_ijd) * 24 * 3600

m_bkg = ijd2plot(lc.data_unit[1].data["TIME"]) < background_age
np.sum(m_bkg)
bkg = np.mean(lc.data_unit[1].data["FLUX"][m_bkg])
bkg

from matplotlib import pylab as plt

plt.figure(figsize=(10, 6))

plt.errorbar(
    ijd2plot(lc.data_unit[1].data["TIME"]),
    lc.data_unit[1].data["FLUX"],
    lc.data_unit[1].data["ERROR"],
    label=f"{lc_time_scale} s",
    alpha=0.5,
)

best_detection = {"sn": None, "timescale": None, "t_max_sn": None}

for detection_time_scale in detection_time_scales.split(","):
    n = int(float(detection_time_scale.strip()) / lc_time_scale)

    r_t = rebin(lc.data_unit[1].data["TIME"], n) / n
    r_c = rebin(lc.data_unit[1].data["FLUX"], n) / n
    r_ce = rebin(lc.data_unit[1].data["ERROR"] ** 2, n) ** 0.5 / n

    sn = (r_c - bkg) / r_ce
    imax = sn.argmax()

    if best_detection["sn"] is None or best_detection["sn"] < sn[imax]:
        best_detection["sn"] = sn[imax]
        best_detection["timescale"] = detection_time_scale
        best_detection["t_max_sn"] = Time(r_t[imax], format="mjd").isot

    plt.errorbar(
        ijd2plot(r_t),
        r_c,
        r_ce,
        label=f"{detection_time_scale} s, S/N = {sn[imax]:.1f}",
    )

plt.axhline(bkg, lw=3, ls=":", c="k")

plt.grid()
plt.legend()

plt.ylabel("Flux")
plt.xlabel("Time, MJD")

plt.title(best_detection)

plt.savefig("annotated-lc.png")

Time(r_t, format="mjd").isot

bin_image = PictureProduct.from_file("annotated-lc.png")

detection_note = str(dict(best_detection=best_detection))

import numpy as np
from astropy import wcs
from astropy.io import fits
from matplotlib import pylab as plt

x, y = np.meshgrid(np.linspace(0, 10, 500), np.linspace(0, 10, 500))

image = np.exp(-((x - 5) ** 2 + (y - 5) ** 2))

pixsize = 0.1
Npix = 500
cdec = 1

RA = 0
DEC = 0

# Create a new WCS object.  The number of axes must be set
# from the start
w = wcs.WCS(naxis=2)

# Set up an "Airy's zenithal" projection
# Vector properties may be set with Python lists, or Numpy arrays
w.wcs.crpix = [image.shape[0] / 2.0, image.shape[1] / 2.0]
w.wcs.cdelt = np.array([pixsize / cdec, pixsize])
w.wcs.crval = [RA, DEC]
w.wcs.ctype = ["RA---CAR", "DEC--CAR"]
w.wcs.set_pv([(1, 1, 45.0)])

# Now, write out the WCS object as a FITS header
header = w.to_header()

# header is an astropy.io.fits.Header object.  We can use it to create a new
# PrimaryHDU and write it to a file.
hdu = fits.PrimaryHDU(image, header=header)
hdu.writeto("Image.fits", overwrite=True)
hdu = fits.open("Image.fits")
im = hdu[0].data
from astropy.wcs import WCS

wcs = WCS(hdu[0].header)
plt.subplot(projection=wcs)
plt.imshow(im, origin="lower")
plt.grid(color="white", ls="solid")
plt.xlabel("RA")
plt.ylabel("Dec")

from oda_api.data_products import ImageDataProduct

fits_image = ImageDataProduct.from_fits_file("Image.fits")

lc = lc  # http://odahub.io/ontology#LightCurve

detection_comment = detection_note  # http://odahub.io/ontology#ODATextProduct
image_output = bin_image  # http://odahub.io/ontology#ODAPictureProduct
image_fits = fits_image  # http://odahub.io/ontology#Image

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(("out_detectgrb_lc", "lc_galaxy.output", lc))
_oda_outs.append(
    (
        "out_detectgrb_detection_comment",
        "detection_comment_galaxy.output",
        detection_comment,
    )
)
_oda_outs.append(
    ("out_detectgrb_image_output", "image_output_galaxy.output", image_output)
)
_oda_outs.append(
    ("out_detectgrb_image_fits", "image_fits_galaxy.output", image_fits)
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
