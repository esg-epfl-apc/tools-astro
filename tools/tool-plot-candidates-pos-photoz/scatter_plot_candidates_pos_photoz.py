#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python

# This script is generated with nb2galaxy

# flake8: noqa

import json
import os
import shutil
import subprocess

import matplotlib.pyplot as pt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.coordinates.matrix_utilities import (
    matrix_product,
    rotation_matrix,
)
from astropy.coordinates.representation import UnitSphericalRepresentation
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from desi import DESILegacySurvey
from matplotlib.colors import LogNorm
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
Radius = 1  # http://odahub.io/ontology#AngleMinutes
pixsize = (
    1.0  # http://odahub.io/ontology#AngleSeconds ; oda:label "Pixel size"
)
data_release = (
    10  # http://odahub.io/ontology#Integer ; oda:label "Data Release"
)

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "C_data_product_" in inp_dic.keys():
    inp_pdic = inp_dic["C_data_product_"]
else:
    inp_pdic = inp_dic
name = str(inp_pdic["name"])
RA = float(inp_pdic["RA"])
DEC = float(inp_pdic["DEC"])
desi_URL = str(inp_pdic["desi_URL"])
photoz_URL = str(inp_pdic["photoz_URL"])
uncertainty = float(inp_pdic["uncertainty"])
s_real_z = float(inp_pdic["s_real_z"])
Radius = float(inp_pdic["Radius"])
pixsize = float(inp_pdic["pixsize"])
data_release = int(inp_pdic["data_release"])

dr = data_release
image_size = Angle(Radius * u.arcmin)
pixsize = Angle(pixsize * u.arcsec)
npix = int(2 * image_size / pixsize)
source = SkyCoord(RA, DEC, unit="degree")

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

ra = ra_dec_h[1].data["RA"]
dec = ra_dec_h[1].data["DEC"]
coords = SkyCoord(ra=ra, dec=dec, unit="deg", frame="icrs")

photoz_h[1].columns, ra_dec_h[1].data.size

photoz_h[1].data["ID"]

photoz = photoz_h[1].data["Z"]

photoz

sep = s_ra_dec.separation(coords).arcsec
indx = np.argmin(sep)
print(photoz[indx], coords.ra[indx], coords.dec[indx])

in_circle = sep < Angle(uncertainty * u.arcsec).arcsec

dict_out = {"GRB": [], "GRB zreal": [], "photoz": [], "RA": [], "DEC": []}
for phz, ra_i, dec_i in zip(
    photoz[in_circle], coords.ra[in_circle], coords.dec[in_circle]
):
    dict_out["GRB"].append(name)
    dict_out["GRB zreal"].append(s_real_z)
    dict_out["photoz"].append(phz)
    dict_out["RA"].append(ra_i.degree)
    dict_out["DEC"].append(dec_i.degree)

t_out = Table.from_pandas(pd.DataFrame(dict_out))
print(t_out)

try:
    query = DESILegacySurvey.get_images(
        position=source,
        survey="dr%d" % dr,
        coordinates="icrs",
        data_release=dr,
        pixels=npix,
        radius=image_size,
        image_band="griz",
    )
    griz, griz_h = query[0][0].data, query[0][0].header
except:
    raise RuntimeError("ERROR: Could not get the image data from LS")

w = WCS(griz_h)
sky = w.pixel_to_world(0, 0, 0)
ra_max_image = sky[0].ra.degree
dec_min_image = sky[0].dec.degree
sky = w.pixel_to_world(npix - 1, npix - 1, npix - 1)
ra_min_image = sky[0].ra.degree
dec_max_image = sky[0].dec.degree

rgb = np.stack([griz[0], griz[1], griz[3]], axis=-1)
rgb = rgb.astype(np.float32)  # Convert to float for proper scaling
rgb = (rgb - np.min(rgb)) / (np.max(rgb) - np.min(rgb))  # Normalize

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

fig, ax = pt.subplots(1, 6, figsize=(35, 5))
lon, lat = create_circle([s_ra_dec.ra, s_ra_dec.dec], s_uncert)

ax[0].scatter(
    coords.ra[indx], coords.dec[indx], color="green", alpha=0.5, s=100
)

ra_lim_min = s_ra_dec.ra - 3 * s_uncert
ra_lim_max = s_ra_dec.ra + 3 * s_uncert
dec_lim_min = s_ra_dec.dec - 3 * s_uncert
dec_lim_max = s_ra_dec.dec + 3 * s_uncert

ax[0].set_xlim([ra_lim_min.value, ra_lim_max.value])
ax[0].set_ylim([dec_lim_min.value, dec_lim_max.value])

for i, txt in enumerate(photoz):
    ax[0].annotate(round(txt, 2), (ra[i], dec[i]))

for ax_ in ax:
    ax_.set_xlabel("RA")
    ax_.set_ylabel("Dec")
    ax_.scatter(s_ra_dec.ra, s_ra_dec.dec, color="magenta", alpha=0.5, s=100)

ax[0].scatter(lon, lat, s=4)
for ax_ in ax[1:]:
    ax_.scatter(lon, lat, s=0.5, alpha=0.5)

ax[0].scatter(ra, dec, color="red", alpha=0.5, s=25)
ax[1].scatter(ra, dec, color="magenta", alpha=0.5, s=2)

for ax_ in ax[2:]:
    ax_.scatter(ra, dec, color="red", alpha=0.5, s=2)

ax[0].set_title(
    f"name={name}; zreal = {s_real_z}; uncertainty={uncertainty} arcsec"
)

print(rgb.shape)
ax[1].set_title(f"GRZ image")
ax[1].imshow(
    rgb,
    origin="lower",
    extent=(ra_max_image, ra_min_image, dec_min_image, dec_max_image),
)

ax[2].set_title(f"G image")
ax[3].set_title(f"R image")
ax[4].set_title(f"I image")
ax[5].set_title(f"Z image")

for i, band_ in enumerate([griz[0], griz[1], griz[2], griz[3]]):
    ax[i + 2].imshow(
        band_,
        norm=LogNorm(vmax=np.max(band_), vmin=np.max(band_) / 1e3),
        origin="lower",
        extent=(ra_max_image, ra_min_image, dec_min_image, dec_max_image),
    )

fig.savefig(f"candidates.png")

from oda_api.data_products import ODAAstropyTable, PictureProduct

bin_image = PictureProduct.from_file("candidates.png")
cat = ODAAstropyTable(t_out)

spectrum_png = bin_image  # http://odahub.io/ontology#ODAPictureProduct
catalog_table = cat  # http://odahub.io/ontology#ODAAstropyTable

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
_oda_outs.append(
    (
        "out_scatter_plot_candidates_pos_photoz_catalog_table",
        "catalog_table_galaxy.output",
        catalog_table,
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
