#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from numpy import sqrt
from oda_api.api import ProgressReporter
from oda_api.data_products import PictureProduct
from oda_api.json import CustomJSONEncoder

pr = ProgressReporter()

src_name = "Crab"  # http://odahub.io/ontology#AstrophysicalObject
RA = 83.628700  # http://odahub.io/ontology#PointOfInterestRA
DEC = 22.014700  # http://odahub.io/ontology#PointOfInterestDEC

T1 = "2000-10-09T13:16:00.0"  # http://odahub.io/ontology#StartTime
T2 = "2024-10-10T13:16:00.0"  # http://odahub.io/ontology#EndTime
Radius_search = 2.0  # http://odahub.io/ontology#AngleDegrees ; oda:label "Cone search radius"
R_s = 0.2  # http://odahub.io/ontology#AngleDegrees ; oda:label "Source region radius for aperture photometry"

Emin = 0.1  # http://odahub.io/ontology#Energy_TeV ; oda:label "Minimal energy" ; oda:group "Plotting"
Emax = 20  # http://odahub.io/ontology#Energy_TeV ; oda:label "Maximal energy" ; oda:group "Plotting"
Slope = 2.7  # http://odahub.io/ontology#Float ; oda:label "Slope of the model powerlaw"
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
    "R_s",
    "Emin",
    "Emax",
    "Slope",
    "Offset",
    "NSB",
]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

racrab = 83.628700
deccrab = 22.014700
crab = SkyCoord(racrab, deccrab, unit="degree")
Coords_s = SkyCoord(RA, DEC, unit="degree")
sep = Coords_s.separation(crab).deg
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
        raise ValueError("NSB level not found")
get_ipython().system("ls {data_dir}")   # noqa: F821

pr.report_progress(stage="data selection", progress=10)

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

hdul = fits.open(DL3_fname[-1])
hdul.info()

hdul["EVENTS"].header

def met2mjd(met):
    # https://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Cicerone/Cicerone_Data/Time_in_ScienceTools.html
    mjd_ref = 52706
    sec_per_day = 86400

    return met / sec_per_day + mjd_ref

met2mjd(333862677.085005), met2mjd(336713240.122429)

ind = 0
pointing = DL3_fname[ind]
hdul = fits.open(pointing)
RMF = hdul["ENERGY DISPERSION"].data
ENERG_LO = RMF["ENERG_LO"][0]
ENERG_HI = RMF["ENERG_HI"][0]
MIGRA_LO = RMF["MIGRA_LO"][0]  # MIGRA_bins=np.linspace(0.2,5,161)
MIGRA_HI = RMF["MIGRA_HI"][0]
MIGRA = (MIGRA_LO + MIGRA_HI) / 2.0
ENERG = sqrt(ENERG_LO * ENERG_HI)
dENERG = ENERG_HI - ENERG_LO

cts_s = []
cts_b = []
Eff_area = []
Eff_area_interp = []
Texp = []
RMFs = []
Tstart = []
Tstop = []
for ind in range(len(DL3_fname)):
    pointing = DL3_fname[ind]
    hdul = fits.open(pointing)
    RA_pnt = hdul[1].header["RA_PNT"]
    DEC_pnt = hdul[1].header["DEC_PNT"]
    Texp.append(hdul[1].header["LIVETIME"])
    Tstart.append(met2mjd(hdul[1].header["TSTART"]))
    Tstop.append(met2mjd(hdul[1].header["TSTOP"]))
    dRA = RA - RA_pnt
    dDEC = DEC - DEC_pnt
    RA_b = RA_pnt - dRA
    DEC_b = DEC_pnt - dDEC
    Coords_b = SkyCoord(RA_b, DEC_b, unit="degree")
    Coords_pnt = SkyCoord(RA_pnt, DEC_pnt, unit="degree")
    dist = Coords_pnt.separation(Coords_s).deg

    RMF = hdul["ENERGY DISPERSION"].data
    mask = RMF["THETA_LO"] < dist
    ind_th = len(RMF["THETA_LO"][mask]) - 1

    AEFF = hdul["EFFECTIVE AREA"].data
    Eff_area.append(AEFF["EFFAREA"][0][ind_th])

    ev = hdul["EVENTS"].data
    ev_ra = ev["RA"]
    ev_dec = ev["DEC"]
    ev_en = ev["ENERGY"]
    ev_time = ev["TIME"]
    ev_coords = SkyCoord(ev_ra, ev_dec, unit="degree")
    sep_s = ev_coords.separation(Coords_s).deg
    sep_b = ev_coords.separation(Coords_b).deg

    mask = (sep_s < R_s) * (ev_en > Emin) * (ev_en < Emax)
    cts_s.append(len(ev_en[mask]))
    mask = (sep_b < R_s) * (ev_en > Emin) * (ev_en < Emax)
    cts_b.append(len(ev_en[mask]))
    hdul.close()

cts_s = np.array(cts_s)
cts_b = np.array(cts_b)
Eff_area = np.array(Eff_area)
Eff_area_interp = np.array(Eff_area_interp)
Texp = np.array(Texp)
Tstart = np.array(Tstart)
Tstop = np.array(Tstart)
Tmid = (Tstart + Tstop) / 2.0
dT = (-Tstart + Tstop) / 2.0

E0 = 1.0

def model_dNdE(E, N, Gam):
    return N * (E / E0) ** (-Gam)

def model_rate(E1, E2, N, Gam):
    dEE = E2 - E1
    EE = sqrt(E1 * E2)
    return model_dNdE(EE, N, Gam) * dEE

def model_cts(N, Gam, ind):
    model_ENERG = model_rate(ENERG_LO, ENERG_HI, N, Gam)
    model_counts_ENERG = model_ENERG * Eff_area[ind] * 1e4 * Texp[ind]
    return model_counts_ENERG

Norm_ref = 1e-11
plt.plot(ENERG_LO, model_cts(Norm_ref, Slope, 1))
plt.xscale("log")
plt.yscale("log")
print(sum(model_cts(Norm_ref, Slope, 1)))

model = []
for i in range(len(cts_s)):
    model.append(sum(model_cts(Norm_ref, Slope, i)))

factor = (
    Norm_ref
    * E0
    / (1 - Slope)
    * ((Emax / E0) ** (1 - Slope) - (Emin / E0) ** (1 - Slope))
)

src = cts_s - cts_b
src_err = sqrt(cts_s + cts_b)
fl = src / model * factor
fl_err = src_err / model * factor

plt.figure(figsize=(12, 5))
plt.errorbar(Tmid, fl, yerr=fl_err, xerr=dT, linestyle="none", marker="o")
plt.xlabel("Time, MJD")
plt.ylabel("Flux, 1/(cm$^2$s)")
plt.savefig("Lightcurve.png", format="png", bbox_inches="tight")

from astropy.table import Table
from oda_api.data_products import ODAAstropyTable

bin_image = PictureProduct.from_file("Lightcurve.png")

data = [Tmid, dT, cts_s, cts_b, fl, fl_err]
names = (
    "T[MJD]",
    "dT[MJD]",
    "Cts_s",
    "Cts_b",
    "Flux[1/cm2s]",
    "Flux_error[1/cm2s]",
)
lc = ODAAstropyTable(Table(data, names=names))

lightcurve_picture = bin_image  # http://odahub.io/ontology#ODAPictureProduct
lightcurve_table = lc  # http://odahub.io/ontology#ODAAstropyTable

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_Lightcurve_public_dl3_lightcurve_picture",
        "lightcurve_picture_galaxy.output",
        lightcurve_picture,
    )
)
_oda_outs.append(
    (
        "out_Lightcurve_public_dl3_lightcurve_table",
        "lightcurve_table_galaxy.output",
        lightcurve_table,
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
