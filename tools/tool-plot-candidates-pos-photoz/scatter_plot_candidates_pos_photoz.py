#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil
import subprocess

import matplotlib.pyplot as pt
import numpy as np
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.coordinates.matrix_utilities import (
    matrix_product,
    rotation_matrix,
)
from astropy.coordinates.representation import UnitSphericalRepresentation
from astropy.io import fits
from oda_api.json import CustomJSONEncoder

## GCN Circular 39418
# name = "GRB 250221A"
# s_ra_dec = SkyCoord(ra="59.46251 deg", dec="-15.13379 deg")
# s_uncert = Angle("1.9 arcsec")
# s_real_z = 0.768

## GCN Circular 38877
## RA = 13:25:18.72, Dec = +25:37:12.7, with photometric redshift z = 0.2-0.3, resulting in a bright trace. We detect Halpha, [N II], and [S II] at z = 0.302, consistent with the Legacy Survey photometric redshift of z = 0.29 +/- 0.04 (Zhou et al., 2021).
# name = "GRB 250108B"
# s_ra_dec = SkyCoord(ra="201.32686 deg", dec="25.61541 deg")
# s_uncert = Angle("2.0 arcsec")
# s_real_z = 2.197

## GCN Circular 38814
# name = "GRB 250103A"
# s_ra_dec = SkyCoord(ra="22.083 deg", dec="-5.096 deg")
# s_uncert = Angle("4.6 arcmin")
# s_real_z = 4.01

## GCN Circular 38759
# name = "GRB 250101A"
# s_ra_dec = SkyCoord(ra="37.06963 deg", dec="19.19521 deg")
# s_uncert = Angle("3.9 arcsec")
# s_real_z = 2.49

## GCN Circular 38637
## photometric redshift determination of the extended object visible in the Legacy Survey (z = 1.20 +/- 0.48; Zhou et al. 2021, MNRAS, 501, 3309)
# name = "GRB 241217A"
# s_ra_dec = SkyCoord(ra="84.1499 deg", dec="-25.3003 deg")
# s_uncert = Angle("20 arcsec")
# s_real_z = 1.879

## GCN Circular 39585
# name = "EP250304a"
# s_ra_dec = SkyCoord(ra="208.3965 deg", dec="-42.8067 deg")
# s_uncert = Angle("20 arcsec")
# s_real_z = 0.2

## GCN Circular 39073
# name = "GRB 250129A"
# s_ra_dec = SkyCoord(ra="198.67654 deg", dec="5.03019 deg")
# s_uncert = Angle("1.9 arcsec")
# s_real_z = 2.15

name = "GRB 250129A"  # http://odahub.io/ontology#AstrophysicalObject
RA = 198.67654  # http://odahub.io/ontology#PointOfInterestRA
DEC = 5.03019  # http://odahub.io/ontology#PointOfInterestDEC
desi_URL = "https://www.astro.unige.ch/~tucci/Phosphoros/MultiBands_Catalog_1k.fits"  # http://odahub.io/ontology#FileReference
photoz_URL = "https://www.astro.unige.ch/~tucci/Phosphoros/MultiBands_Catalog_1k.fits"  # http://odahub.io/ontology#FileReference
uncertainty = 1.9  # http://odahub.io/ontology#float ; oda:label "Uncertainty"
s_real_z = 2.15  # http://odahub.io/ontology#float ; oda:label "Real Redshift"

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in [
    "name",
    "RA",
    "DEC",
    "desi_URL",
    "photoz_URL",
    "uncertainty",
    "s_real_z",
]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

workdir = os.getcwd()
path_tmp = workdir + "/tmp/"
os.makedirs(path_tmp, exist_ok=True)

# get the input catalog and save it into tmp/ directory as Input_Catalog.fits
desi_file = path_tmp + "DESI_Catalog.fits"

read_from_url = False
try:
    desi_output = subprocess.check_output(
        "cp " + desi_URL + " " + desi_file, shell=True
    ).decode()
except:
    print("NOT a file")
    read_from_url = True

if read_from_url:
    try:
        # output = subprocess.check_output('wget -nv -O ' + desi_file + desi_URL, shell=True).decode()
        desi_output = subprocess.check_output(
            "wget -O " + desi_file + ' "' + desi_URL + '"', shell=True
        ).decode()
    except:
        raise RuntimeError("File NOT found")

# get the input catalog and save it into tmp/ directory as Input_Catalog.fits
photoz_file = path_tmp + "PhotoZ_Catalog.fits"

read_from_url = False
try:
    photoz_output = subprocess.check_output(
        "cp " + photoz_URL + " " + photoz_file, shell=True
    ).decode()
except:
    print("NOT a file")
    read_from_url = True

if read_from_url:
    try:
        # output = subprocess.check_output('wget -nv -O ' + photoz_file + photoz_URL, shell=True).decode()
        photoz_output = subprocess.check_output(
            "wget -O " + photoz_file + ' "' + photoz_URL + '"', shell=True
        ).decode()
    except:
        raise RuntimeError("File NOT found")

s_ra_dec = SkyCoord(ra=f"{RA} deg", dec=f"{DEC} deg")
s_uncert = Angle(f"{uncertainty} arcsec")

ra_dec_h = fits.open(f"{desi_file}")
photoz_h = fits.open(f"{photoz_file}")

ra = ra_dec_h[1].data["RA[deg]"]
dec = ra_dec_h[1].data["DEC[deg]"]
coords = SkyCoord(ra=ra, dec=dec, unit="deg", frame="icrs")

photoz_h[1].columns, ra_dec_h[1].data.size

photoz_h[1].data["ID"]

photoz = photoz_h[1].data["Z"]

photoz

sep = s_ra_dec.separation(coords).arcmin
indx = np.argmin(sep)
print(photoz[indx], coords.ra[indx], coords.dec[indx])

def _rotate_polygon(lon, lat, lon0, lat0):
    """
    Given a polygon with vertices defined by (lon, lat), rotate the polygon
    such that the North pole of the spherical coordinates is now at (lon0,
    lat0). Therefore, to end up with a polygon centered on (lon0, lat0), the
    polygon should initially be drawn around the North pole.
    """

    # Create a representation object
    polygon = UnitSphericalRepresentation(lon=lon, lat=lat)

    # Determine rotation matrix to make it so that the circle is centered
    # on the correct longitude/latitude.
    m1 = rotation_matrix(-(0.5 * np.pi * u.radian - lat0), axis="y")
    m2 = rotation_matrix(-lon0, axis="z")
    transform_matrix = matrix_product(m2, m1)

    # Apply 3D rotation
    polygon = polygon.to_cartesian()
    polygon = polygon.transform(transform_matrix)
    polygon = UnitSphericalRepresentation.from_cartesian(polygon)

    return polygon.lon, polygon.lat

def create_circle(center, radius, resolution=100, vertex_unit=u.degree):
    longitude, latitude = center
    lon = np.linspace(0.0, 2 * np.pi, resolution + 1)[:-1] * u.radian
    lat = (
        np.repeat(0.5 * np.pi - radius.to_value(u.radian), resolution)
        * u.radian
    )
    lon, lat = _rotate_polygon(lon, lat, longitude, latitude)
    lon = lon.to_value(vertex_unit)
    lat = lat.to_value(vertex_unit)
    return lon, lat

fig, ax = pt.subplots()
lon, lat = create_circle([s_ra_dec.ra, s_ra_dec.dec], s_uncert)

ax.scatter(s_ra_dec.ra, s_ra_dec.dec, color="red", alpha=0.5, s=100)
ax.scatter(lon, lat)
ax.scatter(coords.ra[indx], coords.dec[indx], color="green", alpha=0.5, s=100)
ax.scatter(ra, dec, color="blue")

ra_lim_min = s_ra_dec.ra - 3 * s_uncert
ra_lim_max = s_ra_dec.ra + 3 * s_uncert
dec_lim_min = s_ra_dec.dec - 3 * s_uncert
dec_lim_max = s_ra_dec.dec + 3 * s_uncert

ax.set_xlim([ra_lim_min.value, ra_lim_max.value])
ax.set_ylim([dec_lim_min.value, dec_lim_max.value])

for i, txt in enumerate(photoz):
    ax.annotate(round(txt, 2), (ra[i], dec[i]))

ax.set_xlabel("RA")
ax.set_ylabel("Dec")
ax.set_title(f"zreal = {s_real_z}; uncertainty={uncertainty} arcsec")
fig.savefig(f"candidates.png")

from oda_api.data_products import PictureProduct

bin_image = PictureProduct.from_file("candidates.png")

spectrum_png = bin_image  # http://odahub.io/ontology#ODAPictureProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_scatter_plot_candidates_pos_photoz_spectrum_png",
        "spectrum_png_galaxy.output",
        spectrum_png,
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
