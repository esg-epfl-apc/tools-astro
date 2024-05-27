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
from oda_api.data_products import ODAAstropyTable, PictureProduct
from oda_api.json import CustomJSONEncoder

if os.path.exists("hess_dl3_dr1.tar.gz") == False:
    get_ipython().system(   # noqa: F821
        "wget https://zenodo.org/record/1421099/files/hess_dl3_dr1.tar.gz"
    )
    get_ipython().system("tar -zxvf hess_dl3_dr1.tar.gz")   # noqa: F821
from oda_api.api import ProgressReporter

pr = ProgressReporter()
pr.report_progress(stage="Progress", progress=5.0)

src_name = "Crab"  # http://odahub.io/ontology#AstrophysicalObject
RA = 83.628700  # http://odahub.io/ontology#PointOfInterestRA
DEC = 22.014700  # http://odahub.io/ontology#PointOfInterestDEC
T1 = "2004-11-20T13:16:00.0"  # http://odahub.io/ontology#StartTime
T2 = "2004-12-20T13:16:00.0"  # http://odahub.io/ontology#EndTime
Radius = 2.5  # http://odahub.io/ontology#AngleDegrees
R_s = 0.2  # http://odahub.io/ontology#AngleDegrees
Emin = 1  # http://odahub.io/ontology#Energy_TeV
Emax = 100.0  # http://odahub.io/ontology#Energy_TeV
NTbins = 30  # http://odahub.io/ontology#Integer

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
Tbegs = []
for i in mask:
    OBSlist.append(DL3_files[i])
    Tbegs.append(Tstart[i])
if len(OBSlist) == 0:
    message = "No data found"
    raise RuntimeError("No data found")
message

Tbins = np.linspace(T1, T2, NTbins + 1)
Tmin = Tbins[:-1]
Tmax = Tbins[1:]
Tmean = (Tmin + Tmax) / 2.0
Tbins

flux = np.zeros(NTbins)
flux_err = np.zeros(NTbins)
flux_b = np.zeros(NTbins)
flux_b_err = np.zeros(NTbins)
Expos = np.zeros(NTbins)
src_cts = np.zeros(NTbins)
bkg_cts = np.zeros(NTbins)
for count, f in enumerate(OBSlist):
    pr.report_progress(
        stage="Progress", progress=5.0 + 95 * count / len(OBSlist)
    )
    hdul = fits.open("data/" + f)
    RA_pnt = hdul[1].header["RA_PNT"]
    DEC_pnt = hdul[1].header["DEC_PNT"]
    Texp = hdul[1].header["LIVETIME"]
    Trun_start = hdul[1].header["TSTART"]
    dRA = RA - RA_pnt
    dDEC = DEC - DEC_pnt
    RA_b = RA_pnt - dRA
    DEC_b = DEC_pnt - dDEC
    Coords_b = SkyCoord(RA_b, DEC_b, unit="degree")
    Coords_pnt = SkyCoord(RA_pnt, DEC_pnt, unit="degree")
    dist = Coords_pnt.separation(Coords_s).deg

    ev = hdul["EVENTS"].data
    ev_ra = ev["RA"]
    ev_dec = ev["DEC"]
    ev_en = ev["ENERGY"]
    ev_time = (ev["TIME"] - Trun_start) / 86400.0 + Tbegs[count]
    Nevents = len(ev)
    print(ev_time[0])
    ev_coords = SkyCoord(ev_ra, ev_dec, unit="degree")
    sep_s = ev_coords.separation(Coords_s).deg
    sep_b = ev_coords.separation(Coords_b).deg

    hdu = hdul["AEFF"].data
    EEmin = hdu["ENERG_LO"][0]
    EEmax = hdu["ENERG_HI"][0]
    EE = sqrt(EEmin * EEmax)
    EEbins = np.concatenate((EEmin, [EEmax[-1]]))
    AA = hdu["EFFAREA"][0] + 1e-10
    Thmin = hdu["THETA_LO"][0]
    Thmax = hdu["THETA_HI"][0]
    ind = np.argmin((Thmin - dist) ** 2)
    mask = EE < Emin
    ind_en = len(EE[mask])
    Expos += np.histogram(
        ev_time,
        weights=AA[ind, ind_en] * Texp * 1e4 * np.ones(Nevents) / Nevents,
        bins=Tbins,
    )[0]
    mask = np.where((sep_s < R_s) & (ev_en > Emin) & (ev_en < Emax))
    src_cts += np.histogram(ev_time[mask], bins=Tbins)[0]
    mask = np.where((sep_b < R_s) & (ev_en > Emin) & (ev_en < Emax))
    bkg_cts += np.histogram(ev_time[mask], bins=Tbins)[0]
    hdul.close()
print(src_cts)
print(bkg_cts)
print(Expos)
src = src_cts - bkg_cts
src_err = sqrt(src_cts + bkg_cts)
flux = src / (Expos + 1)
flux_err = src_err / (Expos + 1)
flux_b = bkg_cts / (Expos + 1)
flux_b_err = sqrt(bkg_cts) / (Expos + 1)
# flux_err+=np.sum(src_err/AA,axis=1)
# flux_b+=np.sum(cts2/AA,axis=1)
# flux_b_err+=np.sum(sqrt(cts2)/AA,axis=1)
print(flux)

if message == "":
    plt.errorbar(
        Tmean,
        flux,
        yerr=flux_err,
        xerr=[Tmean - Tmin, Tmax - Tmean],
        linestyle="none",
        label="source",
    )
    plt.errorbar(
        Tmean,
        flux_b,
        yerr=flux_b_err,
        xerr=[Tmean - Tmin, Tmax - Tmean],
        linestyle="none",
        label="background",
    )
    plt.xlabel("Time, MJD")
    plt.ylabel("Flux, cts/cm$^2$s")
    plt.yscale("log")
    ymin = min(min(flux - flux_err), min(flux_b - flux_b_err))
    ymax = max(max(flux + flux_err), max(flux_b + flux_b_err))
    plt.ylim(ymin / 2.0, 2 * ymax)
    plt.legend(loc="lower left")
    plt.savefig("Lightcurve.png", format="png")

if message == "":
    bin_image = PictureProduct.from_file("Lightcurve.png")
    from astropy.table import Table

    data = [Tmean, Tmin, Tmax, flux, flux_err, flux_b, flux_b_err]
    names = (
        "Tmean[MJD]",
        "Tmin[MJD]",
        "Tmax[MJD]",
        "Flux[counts/cm2s]",
        "Flux_error[counts/cm2s]",
        "Background[counts/cm2s]",
        "Background_error[counts/cm2s]",
    )
    lc = ODAAstropyTable(Table(data, names=names))

png = bin_image  # http://odahub.io/ontology#ODAPictureProduct
table = lc  # http://odahub.io/ontology#ODAAstropyTable

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(("out_Lightcurve_png", "png_galaxy.output", png))
_oda_outs.append(("out_Lightcurve_table", "table_galaxy.output", table))

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
