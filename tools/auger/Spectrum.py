#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
from oda_api.json import CustomJSONEncoder

get_ipython().run_line_magic("matplotlib", "inline")   # noqa: F821

import datetime
import os
from zipfile import ZipFile

from astropy.coordinates import SkyCoord
from astropy.time import Time
from numpy import arccos, cos, log10, pi, sin, sqrt

src_name = "Cen A"  # http://odahub.io/ontology#AstrophysicalObject
RA = 201.365063  # http://odahub.io/ontology#PointOfInterestRA
DEC = -43.019113  # http://odahub.io/ontology#PointOfInterestDEC

T1 = "2021-10-09T13:16:00.0"  # http://odahub.io/ontology#StartTime
T2 = "2021-10-13T13:16:00.0"  # http://odahub.io/ontology#EndTime

Source_region_radius = 27.0  # http://odahub.io/ontology#AngleDegrees ; oda:label "Source signal region radius"
Emin = 31.62e18  # http://odahub.io/ontology#Energy_eV ; oda:label "Minimal energy"
Emax = 316.2e18  # http://odahub.io/ontology#Energy_eV ; oda:label "Maximal energy"
NEbins = (
    2  # http://odahub.io/ontology#Integer ; oda:label "Number of energy bins"
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
    "T1",
    "T2",
    "Source_region_radius",
    "Emin",
    "Emax",
    "NEbins",
]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

if DEC > 30.0:
    raise ValueError("No exposure for this declination")

t = Time([T1, T2], format="isot", scale="utc")

t = t.mjd
Tstart = t[0]
Tstop = t[1]

coords_s = SkyCoord(RA, DEC, unit="degree")

Ebins = np.logspace(log10(Emin), log10(Emax), NEbins + 1)

if os.path.exists("Anisotropy_ ApJS_2022") == False:
    #!zenodo_get 10.5281/zenodo.6759610
    get_ipython().system(   # noqa: F821
        "wget 'https://zenodo.org/record/6759610/files/Anisotropy_ ApJS_2022.zip'"
    )
    ZipFile("Anisotropy_ ApJS_2022.zip", "r").extractall()
else:
    print("Data in place")

d = np.genfromtxt(
    "Anisotropy_ ApJS_2022/Data/AugerApJS2022_Yr_JD_UTC_Th_Ph_RA_Dec_E_Expo.dat"
)
yr = np.array(d[:, 0])
jd = d[:, 1]
utc = np.array(d[:, 2])
th = np.array(d[:, 3])
ph = d[:, 4]
ra_ev = d[:, 5]
dec_ev = d[:, 6]
E = np.array(d[:, 7]) * 1e18
expos = d[:, 8]

utc

t_ev = []
for tt in utc:
    t_ev.append(datetime.datetime.fromtimestamp(tt).isoformat())
t_ev = Time(t_ev, format="isot", scale="utc").mjd
mask = (t_ev > Tstart) & (t_ev < Tstop)
if sum(mask) > 0:
    th = th[mask]
    ph = ph[mask]
    ra_ev = ra_ev[mask]
    dec_ev = dec_ev[mask]
    E = E[mask]
    expos = expos[mask]
    tot_expos = expos[-1] - expos[0]
else:
    raise ValueError("No events found in this time interval")

coords = SkyCoord(ra_ev, dec_ev, unit="degree")
sep = coords.separation(coords_s).deg

cts_s = np.zeros(NEbins)
cts_b = np.zeros(NEbins)
mask = sep < Source_region_radius
cts_s[i] = len(E[mask])
plt.scatter(ra_ev[mask], dec_ev[mask], color="green", label="source events")
mask = (
    (dec_ev < DEC + Source_region_radius)
    & (dec_ev > DEC - Source_region_radius)
    & (sep > Source_region_radius)
)
cts_b[i] = len(E[mask])
plt.scatter(ra_ev[mask], dec_ev[mask], color="blue", label="background events")

plt.scatter([RA], [DEC], color="red")
plt.scatter(ra_ev, dec_ev, alpha=0.1)
plt.legend(loc="upper right")
plt.xlabel("RA, degrees")
plt.ylabel("DEC, degrees")

cts_s = np.zeros(NEbins)
cts_b = np.zeros(NEbins)
for i in range(NEbins):
    mask = (E > Ebins[i]) & (E < Ebins[i + 1]) & (sep < Source_region_radius)
    cts_s[i] = len(E[mask])
    plt.scatter(
        ra_ev[mask], dec_ev[mask], color="green", label="source events"
    )
    mask = (
        (E > Ebins[i])
        & (E < Ebins[i + 1])
        & (dec_ev < DEC + Source_region_radius)
        & (dec_ev > DEC - Source_region_radius)
        & (sep > Source_region_radius)
    )
    cts_b[i] = len(E[mask])
    plt.scatter(
        ra_ev[mask], dec_ev[mask], color="blue", label="background events"
    )

plt.scatter([RA], [DEC], color="red")
plt.scatter(ra_ev, dec_ev, alpha=0.1)
plt.legend(loc="upper right")

# exposure function * cos(declination) for a given zenith range
def fCos(declination, latitude, thetamin, thetamax):
    declination *= pi / 180.0
    latitude *= pi / 180.0
    thetamin *= pi / 180.0
    thetamax *= pi / 180.0
    ksiM = (cos(thetamax) - sin(latitude) * sin(declination)) / (
        cos(latitude) * cos(declination)
    )
    ksim = (cos(thetamin) - sin(latitude) * sin(declination)) / (
        cos(latitude) * cos(declination)
    )
    aM = 0
    am = 0
    if ksiM < -1:
        aM = pi
    elif (ksiM >= -1) and (ksiM <= 1):
        aM = arccos(ksiM)
    if ksim < -1:
        am = pi
    elif (ksim >= -1) and (ksim <= 1):
        am = arccos(ksim)

    return (
        cos(latitude) * cos(declination) * (sin(aM) - sin(am))
        + (aM - am) * sin(latitude) * sin(declination)
    ) * cos(declination)

def Sum(declination, thetamax_v, latitude, norm_v, thetamin, thetamax, norm_i):
    declination *= pi / 180.0
    latitude *= pi / 180.0
    thetamin *= pi / 180.0
    thetamax *= pi / 180.0
    thetamax_v *= pi / 180.0

    # vertical part
    ksiM = (cos(thetamax_v) - sin(latitude) * sin(declination)) / (
        cos(latitude) * cos(declination)
    )
    aM = 0
    if ksiM < -1:
        aM = pi
    elif (ksiM >= -1) and (ksiM <= 1):
        aM = arccos(ksiM)
    omegav = (
        (
            cos(latitude) * cos(declination) * sin(aM)
            + aM * sin(latitude) * sin(declination)
        )
    ) * norm_v
    # print(omegav)
    # inclined part

    aM = 0
    am = 0
    ksiM = (cos(thetamax) - sin(latitude) * sin(declination)) / (
        cos(latitude) * cos(declination)
    )
    ksim = (cos(thetamin) - sin(latitude) * sin(declination)) / (
        cos(latitude) * cos(declination)
    )
    if ksiM < -1:
        aM = pi
    elif (ksiM >= -1) and (ksiM <= 1):
        aM = arccos(ksiM)
    if ksim < -1:
        am = pi
    elif (ksim >= -1) and (ksim <= 1):
        am = arccos(ksim)
    omegah = (
        cos(latitude) * cos(declination) * (sin(aM) - sin(am))
        + (aM - am) * sin(latitude) * sin(declination)
    ) * norm_i

    return (omegav + omegah) * cos(declination)

# number of vertical & inclined events
nVertOfficial = 2040.0
nHorOfficial = 595.0
nVert = 0.0
nHor = 0.0

for i in range(len(E)):
    if th[i] < 60.0:
        nVert += 1
    else:
        nHor += 1

# zenith angle range
thetaMin = 60.0
thetaMax = 80.0

# Auger latitude site
augerLat = -35.2
print(nHor, nVert)

ddec = 1
Nbins = int(180 / ddec)
dec_bins = np.linspace(-90, 90, Nbins + 1)
decs = (dec_bins[1:] + dec_bins[:-1]) / 2.0
covh = np.zeros(len(decs))
covv = np.zeros(len(decs))
for i in range(len(decs)):
    covh[i] = fCos(decs[i], augerLat, thetaMin, thetaMax)
    covv[i] = fCos(decs[i], augerLat, 0, thetaMin)
covtot = covv + covh
plt.plot(decs, covh * 120)
plt.plot(decs, covv * 120)
intexpoh = sum(covh * ddec) * pi / 180.0
intexpov = sum(covv * ddec) * pi / 180.0

normHor = nHor / intexpoh / 2 / pi
normHor = nVert / intexpov / 2 / pi

print(nHor / nVert, intexpoh / intexpov)

# covh/=cos(decs[i]*pi/180.)
# covv/=cos(decs[i]*pi/180.)
covtot = covv + covh
dOm = 2 * pi * cos(decs * pi / 180.0)
plt.plot(decs, dOm * 10)
plt.plot(decs, covtot * 120, linewidth=4)

# cts=AT*dOm*F
# AT*dOm=covtot*Norm
# sum(AT*Om)=tot_exp=sum(covtot*Norm)
Norm = tot_expos / sum(covtot)
print(Norm)

# exposure as a function of declinaiton
AT = covtot * Norm / dOm

# check: total exposure in (km2 yr sr)
print(sum(AT * dOm), tot_expos)

plt.plot(decs, AT)
plt.xlabel("Declinaiton, degrees")
plt.ylabel("Exposure, km$^2$yr")

plt.plot(decs, AT * dOm)

plt.xlabel("Declinaiton, degrees")
plt.ylabel("Exposure, km$^2$yr sr")
print(sum(AT * dOm))

F = 2.2e-2  # events/km2 yr sr

plt.plot(decs, F * AT * dOm)

plt.xlabel("Declinaiton, degrees")
plt.ylabel("Exposure, km$^2$yr sr")
print(sum(AT * dOm))

h = plt.hist(dec_ev, bins=dec_bins, alpha=0.2)
cts = h[0]
plt.errorbar(decs, cts, yerr=sqrt(cts))

def Expos(dd):
    i = 0
    while dec_bins[i] < dd:
        i += 1
    return AT[i - 1]

Expos(-89.1)

Omega_s = pi * (1 - cos(Source_region_radius * pi / 180.0))
Omega_b = pi * (
    cos((DEC + Source_region_radius) * pi / 180.0)
    - cos((DEC - Source_region_radius) * pi / 180.0)
)
factor = Omega_s / (Omega_b - Omega_s)
Src = cts_s - factor * cts_b
Src_err = sqrt(cts_s + factor**2 * cts_b)
Src_err_lo = sqrt(cts_s + factor**2 * cts_b + 0.25) - 0.5
Src_err_hi = sqrt(cts_s + factor**2 * cts_b + 0.25) + 0.5
print("Source counts", cts_s)
print("Background estimate", factor * cts_b)
print(Src)
print(Src_err)
print(Src_err_lo)
print(Src_err_hi)
Exposure = Expos(DEC) * 3e7 * 1e10
E_min = Ebins[:-1]
E_max = Ebins[1:]
E_mean = sqrt(E_min * E_max)
F = Src / (E_max - E_min) * E_max * E_min / Exposure * 1.6e-12
F_err_lo = (
    Src_err_lo / (E_max - E_min) * E_max * E_min / Exposure * 1.6e-12
)  # erg.cm2s
F_err_hi = (
    Src_err_hi / (E_max - E_min) * E_max * E_min / Exposure * 1.6e-12
)  # erg.cm2s
uplims = F - F_err_lo < 0
print(F, F_err_lo, F_err_hi, uplims)
F = F * (1 - uplims) + 2 * F_err_hi * uplims
F_err_lo = F_err_lo * (1 - uplims) + 0.5 * F * uplims
F_err_hi = F_err_hi * (1 - uplims) + F * uplims
print(F, F_err_lo, F_err_hi, uplims)

plt.figure(figsize=(10, 7))
plt.errorbar(
    E_mean,
    F,
    yerr=[F_err_lo, F_err_hi],
    xerr=[E_mean - E_min, E_max - E_mean],
    uplims=uplims,
    linestyle="none",
    marker="o",
)
plt.xscale("log")
plt.yscale("log")
plt.axhline(0)
plt.xlabel("$E$, eV", fontsize=16)
plt.ylabel("$E^2 dN/dE$, erg/(cm$^2$s)", fontsize=16)
plt.tick_params(labelsize=16)
plt.ylim(0.5 * min(F - F_err_lo), 2 * max(F + F_err_hi))
plt.savefig("Auger_spectrum.png", format="png", bbox_inches="tight")

from oda_api.data_products import ODAAstropyTable, PictureProduct

bin_image = PictureProduct.from_file("Auger_spectrum.png")

from astropy.table import Table

data = [E_mean, E_min, E_max, F, F_err_lo, F_err_hi, uplims, cts_s, cts_b]
names = (
    "Emean[eV]",
    "Emin[eV]",
    "Emax[eV]",
    "Flux[erg/cm2s]",
    "Flux_error_lo[erg/cm2s]",
    "Flux_error_hi[erg/cm2s]",
    "UpLim",
    "Cts_src",
    "Cts_bkg",
)
spec = ODAAstropyTable(Table(data, names=names))

picture_png = bin_image  # http://odahub.io/ontology#ODAPictureProduct
spectrum_astropy_table = spec  # http://odahub.io/ontology#ODAAstropyTable

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    ("out_Spectrum_picture_png", "picture_png_galaxy.output", picture_png)
)
_oda_outs.append(
    (
        "out_Spectrum_spectrum_astropy_table",
        "spectrum_astropy_table_galaxy.output",
        spectrum_astropy_table,
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
