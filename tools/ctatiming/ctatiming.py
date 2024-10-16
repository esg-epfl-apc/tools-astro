#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

event_file = "events.fits"  # oda:POSIXPath
tbin_s = 1000.0  # oda:Float

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

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from matplotlib import pyplot as plt
from oda_api.api import PictureProduct

events = fits.open(event_file)[1].data
evt_t = events["TIME"]
evt_e = events["ENERGY"]
evt_ra = events["RA"]
evt_dec = events["DEC"]
evt_c = SkyCoord(evt_ra, evt_dec, frame="icrs", unit="deg")

h2 = plt.hist2d(evt_c.ra.deg, evt_c.dec.deg, bins=100)

peak = np.unravel_index(h2[0].argmax(), h2[0].shape)
peak_ra = h2[1][peak[0]]
peak_dec = h2[2][peak[1]]
center_ra = np.mean(evt_ra)
center_dec = np.mean(evt_dec)

off_region_ra = center_ra - (peak_ra - center_ra)
off_region_dec = center_dec - (peak_dec - center_dec)

plt.scatter(peak_ra, peak_dec, color="red")
plt.scatter(center_ra, center_dec, color="white")
plt.scatter(off_region_ra, off_region_dec, color="cyan")

m_on = (
    SkyCoord(peak_ra, peak_dec, frame="icrs", unit="deg").separation(evt_c).deg
    < 0.5
)
m_off = (
    SkyCoord(off_region_ra, off_region_dec, frame="icrs", unit="deg")
    .separation(evt_c)
    .deg
    < 0.5
)

bysource = {}

for n, m in [("on-source", m_on), ("off-source", m_off)]:
    h = np.histogram(evt_t[m], bins=int((evt_t.max() - evt_t.min()) / tbin_s))

    rate = h[0] / tbin_s
    rate_err = np.sqrt(h[0]) / tbin_s
    t = 0.5 * (h[1][1:] + h[1][:-1])

    bysource[n] = (t, rate, rate_err)

bysource["on-source-nobkg"] = (
    bysource["off-source"][0],
    bysource["on-source"][1] - bysource["off-source"][1],
    np.sqrt(bysource["on-source"][2] ** 2 + bysource["off-source"][2] ** 2),
)

plt.subplots(2, 1, figsize=(10, 10), sharex=True, gridspec_kw={"hspace": 0})

plt.sca(plt.gcf().get_axes()[1])
plt.grid()

for n, (t, rate, rate_err) in bysource.items():
    plt.errorbar(
        t, rate / rate_err, yerr=rate_err / rate_err, fmt="o", label=n
    )

plt.ylabel("S/N")

plt.sca(plt.gcf().get_axes()[0])
plt.grid()
for n, (t, rate, rate_err) in bysource.items():
    plt.errorbar(t, rate, yerr=rate_err, fmt="o", label=n)

plt.ylabel("Rate (cts/s)")
plt.xlabel("Time (s)")
plt.legend()

plt.savefig("picture.png")

output = PictureProduct.from_file("picture.png")  # oda:ODAPictureProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(("out_ctatiming_output", "output_galaxy.output", output))

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
