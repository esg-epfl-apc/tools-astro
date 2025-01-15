#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astroplan import Observer
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
from astropy.time import Time
from oda_api.json import CustomJSONEncoder

lapalma = Observer.at_site("lapalma")

src_name = "Crab"  # http://odahub.io/ontology#AstrophysicalObject
RA = 83.628700  # http://odahub.io/ontology#PointOfInterestRA
DEC = 22.014700  # http://odahub.io/ontology#PointOfInterestDEC
T1 = "2023-10-09T13:16:00.0"  # http://odahub.io/ontology#StartTime
T2 = "2024-10-10T13:16:00.0"  # http://odahub.io/ontology#EndTime

Zdmin = 0  # http://odahub.io/ontology#AngleDegree ; oda:label "Minimal zenith"
Zdmax = (
    30  # http://odahub.io/ontology#AngleDegree ; oda:label "Maximal zenith"
)
moon_max = 0.1  # #http://odahub.io/ontology#AngleDegree ; oda:label "maximal moon (between 0 and 1)"

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in ["src_name", "RA", "DEC", "T1", "T2", "Zdmin", "Zdmax", "moon_max"]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

src_coords = SkyCoord(RA, DEC, unit="degree")
# source=FixedTarget(name=src_name, coord=src_coords)

# lst_loelevationc = EarthLocation(lat=28.761758*u.deg, lon=-17.890659*u.deg, height=2200*u.m)
lst_loc = EarthLocation(
    lat=lapalma.latitude, lon=lapalma.longitude, height=lapalma.elevation
)
lst_loc

T1 = Time(T1, format="isot", scale="utc")
T2 = Time(T2, format="isot", scale="utc")
dT = T2.mjd - T1.mjd

obstime = (
    T1 + np.linspace(0, dT, int(dT * 24 * 6)) * u.day
)  # ten-minute time steps
frame = AltAz(obstime=obstime, location=lst_loc)

moon_alt = lapalma.moon_altaz(obstime).alt

moon_phase = lapalma.moon_phase(obstime)  # pi is "new", 0 is "full"

moon_alt = moon_alt.deg

moon_phase = moon_phase.value

from numpy import pi

moon_mask = (moon_alt < 0) + ((moon_phase / pi - 1) ** 2 < moon_max**2)

sun_coords = get_sun(obstime)
sun_coords_altaz = sun_coords.transform_to(frame)
Alt_sun = sun_coords_altaz.alt.deg
Az_sun = sun_coords_altaz.az.deg

src_coords_altaz = src_coords.transform_to(frame)
Alt = src_coords_altaz.alt.deg
Az = src_coords_altaz.az.deg

Zd = 90 - Alt
mask = (Alt_sun < -5.0) & (Zd < Zdmax) & (Zd > Zdmin) & moon_mask

change = 1 * mask[1:] - 1 * mask[:-1]
m_start = change > 0
m_stop = change < 0
tstarts = obstime[:-1][m_start]
tstops = obstime[:-1][m_stop]
tstarts_mjd = tstarts.mjd
tstops_mjd = tstops.mjd
obs_length = (tstops_mjd - tstarts_mjd) * 24.0
plt.scatter(tstarts_mjd, obs_length, color="blue")
plt.plot(tstarts_mjd, obs_length, color="blue")

plt.xlabel("Time, MJD")
plt.ylabel("Available observing time, hours")
plt.savefig("Obs_time.png", format="png", bbox_inches="tight")
# plt.xlim(tstarts_mjd[0]+30,tstarts_mjd[0]+60)

from astropy.table import Table
from oda_api.data_products import ODAAstropyTable, PictureProduct

data = [tstarts, tstops, obs_length]
names = ("Tstarts", "Tstops", "Obs_length[hr]")
GTI = ODAAstropyTable(Table(data, names=names))
bin_image = PictureProduct.from_file("Obs_time.png")

GTI_png = bin_image  # http://odahub.io/ontology#ODAPictureProduct
GTI_table = GTI  # http://odahub.io/ontology#ODAAstropyTable

Table(data, names=names)[:100]

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    ("out_Visibility_predictor_GTI_png", "GTI_png_galaxy.output", GTI_png)
)
_oda_outs.append(
    (
        "out_Visibility_predictor_GTI_table",
        "GTI_table_galaxy.output",
        GTI_table,
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
