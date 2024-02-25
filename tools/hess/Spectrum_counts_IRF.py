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
from numpy import log10, sqrt
from oda_api.data_products import ODAAstropyTable, PictureProduct
from oda_api.json import CustomJSONEncoder

if os.path.exists("hess_dl3_dr1.tar.gz") == False:
    get_ipython().system(   # noqa: F821
        "wget https://zenodo.org/record/1421099/files/hess_dl3_dr1.tar.gz"
    )
    get_ipython().system("tar -zxvf hess_dl3_dr1.tar.gz")   # noqa: F821

src_name = "Crab"  # http://odahub.io/ontology#AstrophysicalObject
RA = 83.628700  # http://odahub.io/ontology#PointOfInterestRA
DEC = 22.014700  # http://odahub.io/ontology#PointOfInterestDEC
# src_name='PKS 2155-304'
# RA = 329.716938  # http://odahub.io/ontology#PointOfInterestRA
# DEC = -30.225588 # http://odahub.io/ontology#PointOfInterestDEC

T1 = "2000-10-09T13:16:00.0"  # http://odahub.io/ontology#StartTime
T2 = "2022-10-10T13:16:00.0"  # http://odahub.io/ontology#EndTime
Radius = 2.5  # http://odahub.io/ontology#AngleDegrees
R_s = 0.2  # http://odahub.io/ontology#AngleDegrees

Emin = 0.1  # http://odahub.io/ontology#Energy_TeV
Emax = 100.0  # http://odahub.io/ontology#Energy_TeV
NEbins = 30  # http://odahub.io/ontology#Integer

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

E0 = 1.0

def model_dNdE(E, N, Gam):
    return N * (E / E0) ** (Gam)

def model_rate(E1, E2, N, Gam):
    dEE = E2 - E1
    EE = sqrt(E1 * E2)
    return model_dNdE(EE, N, Gam) * dEE

Ebins = np.logspace(log10(Emin), log10(Emax), NEbins + 1)
Emins = Ebins[:-1]
Emaxs = Ebins[1:]
Emeans = sqrt(Emins * Emaxs)
lgEmeans = log10(Emeans)
dE = Ebins[1:] - Ebins[:-1]
Emeans

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
for i in mask:
    OBSlist.append(DL3_files[i])
if len(OBSlist) == 0:
    message = "No data found"
    raise RuntimeError("No data found")
offaxis = seps[mask]
Tstart = np.array(Tstart)[mask]
offaxis, Tstart

ind = 0
pointing = OBSlist[ind]
hdul = fits.open("data/" + pointing)
RMF = hdul["EDISP"].data
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
for ind in range(len(OBSlist)):
    pointing = OBSlist[ind]
    hdul = fits.open("data/" + pointing)

    RA_pnt = hdul[1].header["RA_PNT"]
    DEC_pnt = hdul[1].header["DEC_PNT"]
    Texp.append(hdul[1].header["LIVETIME"])
    dRA = RA - RA_pnt
    dDEC = DEC - DEC_pnt
    RA_b = RA_pnt - dRA
    DEC_b = DEC_pnt - dDEC
    Coords_b = SkyCoord(RA_b, DEC_b, unit="degree")
    Coords_pnt = SkyCoord(RA_pnt, DEC_pnt, unit="degree")
    dist = Coords_pnt.separation(Coords_s).deg

    RMF = hdul["EDISP"].data
    mask = RMF["THETA_LO"] < dist
    ind_th = len(RMF["THETA_LO"][mask]) - 1
    RMF_th = RMF["MATRIX"][0][ind_th]
    RMF_interp = np.zeros((len(Emeans), len(ENERG)))
    for k in range(len(ENERG)):
        dp_dErec = RMF_th[:, k] / ENERG[k]
        Erec = MIGRA * ENERG[k]
        dp_dErec_interp = (
            np.interp(Emeans, Erec, dp_dErec)
            * (Emeans > min(Erec))
            * (Emeans < max(Erec))
        )
        RMF_interp[:, k] = dp_dErec_interp
    RMFs.append(RMF_interp)

    AEFF = hdul["AEFF"].data
    Eff_area.append(AEFF["EFFAREA"][0][ind_th])
    Eff_area_interp.append(np.interp(Emeans, ENERG, Eff_area[-1]))

    ev = hdul["EVENTS"].data
    ev_ra = ev["RA"]
    ev_dec = ev["DEC"]
    ev_en = ev["ENERGY"]
    ev_time = ev["TIME"]
    ev_coords = SkyCoord(ev_ra, ev_dec, unit="degree")
    sep_s = ev_coords.separation(Coords_s).deg
    sep_b = ev_coords.separation(Coords_b).deg

    mask = sep_s < R_s
    cts_s.append(np.histogram(ev_en[mask], bins=Ebins)[0])
    mask = sep_b < R_s
    cts_b.append(np.histogram(ev_en[mask], bins=Ebins)[0])
    hdul.close()

cts_s = np.array(cts_s)
cts_b = np.array(cts_b)
Eff_area = np.array(Eff_area)
Eff_area_interp = np.array(Eff_area_interp)
Texp = np.array(Texp)

cts_s_tot = sum(cts_s)
cts_b_tot = sum(cts_b)
src_tot = cts_s_tot - cts_b_tot
src_tot_err = sqrt(cts_s_tot + cts_b_tot)
print(sum(cts_s_tot))
Expos_tot_interp = sum(Eff_area_interp * np.outer(Texp, np.ones(NEbins))) * 1e4
Expos_tot = sum(Eff_area * np.outer(Texp, np.ones(len(ENERG)))) * 1e4

flux_tot = src_tot / (Expos_tot_interp + 1) / (Emaxs - Emins) * Emaxs * Emins
flux_tot_err = (
    src_tot_err / (Expos_tot_interp + 1) / (Emaxs - Emins) * Emaxs * Emins
)
plt.errorbar(
    Emeans,
    flux_tot,
    flux_tot_err,
    xerr=[Emeans - Emins, Emaxs - Emeans],
    linestyle="none",
)

d = np.genfromtxt("Crab_spectrum.csv")
plt.plot(d[:, 0], d[:, 1])

plt.xscale("log")
plt.yscale("log")

N = 5.5e-12 * 7
Gam = -2.68
model_ENERG = model_rate(ENERG_LO, ENERG_HI, N, Gam)

def model_cts_Erec(N, Gam):
    model_ENERG = model_rate(ENERG_LO, ENERG_HI, N, Gam)
    res = np.zeros(NEbins)
    for ind in range(len(OBSlist)):
        model_counts_ENERG = model_ENERG * Eff_area[ind] * 1e4 * Texp[ind]
        for k in range(len(ENERG)):
            res += model_counts_ENERG[k] * RMFs[ind][:, k] * dE
    return res

model_counts_Emeans = model_cts_Erec(N, Gam)
plt.plot(Emeans, model_counts_Emeans)
plt.errorbar(
    Emeans,
    src_tot,
    src_tot_err,
    xerr=[Emeans - Emins, Emaxs - Emeans],
    linestyle="none",
)

plt.xscale("log")
plt.yscale("log")
plt.xlim(0.3, 30)
plt.ylim(0.1, 1e4)

flux = src_tot / (Expos_tot_interp + 1) / (Emaxs - Emins) * Emins * Emaxs
flux_err = (
    src_tot_err / (Expos_tot_interp + 1) / (Emaxs - Emins) * Emins * Emaxs
)
d = np.genfromtxt("Crab_spectrum.csv")
plt.plot(d[:, 0], d[:, 1])
plt.errorbar(
    Emeans,
    flux,
    flux_err,
    xerr=[Emeans - Emins, Emaxs - Emeans],
    linestyle="none",
)

plt.plot(Emeans, model_dNdE(Emeans, N, Gam) * Emeans**2)
plt.xscale("log")
plt.yscale("log")

def chi2(p):
    N, slope = p
    counts = model_cts_Erec(N, slope)
    print(sum(counts))
    m = Emeans > e_min
    m &= Emeans < e_max
    chi2 = (((counts[m] - src_tot[m]) / src_tot_err[m]) ** 2).sum()
    # print(N,slope,chi2)
    return chi2

N = 5.5e-12 * 7
Gam = -2.68
chi2([N, Gam])

def chi2(p):
    N, slope = p
    counts = model_cts_Erec(N, slope)
    m = Emeans > e_min
    m &= Emeans < e_max
    chi2 = (((counts[m] - src_tot[m]) / src_tot_err[m]) ** 2).sum()
    # print(N,slope,chi2)
    return chi2

chi2([6e-12, -2.75])
Norm_max = 1e-10
Norm_min = 1e-11
Norm_bins = 100
Gam_min = -2.8
Gam_max = -2.4
Gam_bins = 100
Ns = np.linspace(Norm_min, Norm_max, Norm_bins)
Gams = np.linspace(Gam_min, Gam_max, Gam_bins)
chi2_map = np.zeros((Norm_bins, Gam_bins))
Norm_best = Norm_min
Gam_best = Gam_min
chi2_best = 1e10
for i, N in enumerate(Ns):
    for j, Gam in enumerate(Gams):
        chi2_map[i, j] = chi2([N, Gam])
        if chi2_map[i, j] < chi2_best:
            Norm_best = N
            Gam_best = Gam
            chi2_best = chi2_map[i, j]
print(Norm_best, Gam_best)
# plt.imshow(chi2_map,vmax=np.amin(chi2_map)+4,origin='lower',extent=[Gams[0],Gams[-1],Ns[0],Ns[-1]],aspect=(Gams[-1]-Gams[0])/(Ns[-1]-Ns[0]))

# 68% contour from https://ui.adsabs.harvard.edu/abs/1976ApJ...208..177L/abstract for two-parameter fit

cnt = plt.contour(
    Gams, Ns, chi2_map, levels=[np.amin(chi2_map) + 2.3], colors="red"
)
plt.scatter([Gam_best], [Norm_best], marker="x", color="red")
# plt.colorbar()
print(np.amin(chi2_map))

cont = cnt.get_paths()[0].vertices
gammas = cont[:, 0]
norms = cont[:, 1]

x = np.logspace(-0.5, 1.5, 10)
ymax = np.zeros(10)
ymin = np.ones(10)
for i in range(len(gammas)):
    y = model_dNdE(x, norms[i], gammas[i]) * x**2
    ymax = np.maximum(y, ymax)
    ymin = np.minimum(y, ymin)
    # plt.plot(x,y)
plt.fill_between(x, ymin, ymax, alpha=0.2)

plt.errorbar(
    Emeans,
    flux,
    flux_err,
    xerr=[Emeans - Emins, Emaxs - Emeans],
    linestyle="none",
)

plt.xscale("log")
plt.yscale("log")
d = np.genfromtxt("Crab_spectrum.csv")
plt.plot(d[:, 0], d[:, 1])

col1 = fits.Column(name="E_MIN", format="E", unit="TeV", array=Emins)
col2 = fits.Column(name="E_MAX", format="E", unit="TeV", array=Emaxs)
cols = fits.ColDefs([col1, col2])
hdu = fits.BinTableHDU.from_columns(cols)
hdu.writeto("E_MIN_E_MAX.fits", overwrite=True)

col1 = fits.Column(name="ENERG_LO", format="E", unit="TeV", array=ENERG_LO)
col2 = fits.Column(name="ENERG_HI", format="E", unit="TeV", array=ENERG_HI)
cols = fits.ColDefs([col1, col2])
hdu = fits.BinTableHDU.from_columns(cols)
hdu.writeto("ENERG_LO_ENERG_HI.fits", overwrite=True)

hdu = fits.PrimaryHDU(RMF)
hdu.writeto("MATRIX.fits", overwrite=True)

hdu = fits.PrimaryHDU(AEFF_th)
hdu.writeto("AEFF.fits", overwrite=True)

plt.plot(Emeans, EXPOS_tot)

plt.imshow(RMF_rebin)

plt.imshow(R)

def model(e1, e2, N, slope):
    return N * (e1 / 3.0) ** slope

def convolve(p, return_counts=False):
    N, slope = p
    counts = np.dot(R, model(ENERG_LO, ENERG_HI, N, slope))

    m = Emean > e_min
    m &= Emean < e_max

    chi2 = (((counts - src) / src_err)[m] ** 2).sum()

    # print(p, chi2, chi2 / np.sum(m))

    if return_counts:
        return counts
    else:
        return chi2

N = 1e-8
slope = -2.0

f = minimize(convolve, [N, slope])
print(f)
model_cts = convolve(f.x, True)
plt.plot(Emean, model_cts)

m = Emean > e_min
m &= Emean < e_max
chi2 = (((model_cts - src) / src_err)[m] ** 2).sum()
print("reduced chi2", chi2 / np.sum(m))
# plt.plot(Emean,counts)
plt.errorbar(Emean, src, src_err)
plt.axvline(e_min, color="red")
plt.axvline(e_max, color="red")
plt.xscale("log")
plt.yscale("log")

f.x

ind = 0
pointing = OBSlist[ind]
hdul = fits.open("data/" + pointing)

RA_pnt = hdul[1].header["RA_PNT"]
DEC_pnt = hdul[1].header["DEC_PNT"]
Texp = hdul[1].header["LIVETIME"]
dRA = RA - RA_pnt
dDEC = DEC - DEC_pnt
RA_b = RA_pnt - dRA
DEC_b = DEC_pnt - dDEC
Coords_b = SkyCoord(RA_b, DEC_b, unit="degree")
Coords_pnt = SkyCoord(RA_pnt, DEC_pnt, unit="degree")
dist = Coords_pnt.separation(Coords_s).deg

Texp = hdul[1].header["LIVETIME"]

RMF = hdul["EDISP"].data
RMF.columns
mask = RMF["THETA_LO"] < dist
ind_th = len(RMF["THETA_LO"][mask])

AEFF = hdul["AEFF"].data
AEFF_th = AEFF["EFFAREA"][0][ind_th]

RMF_th = RMF["MATRIX"][0][ind_th]
ENERG_LO = RMF["ENERG_LO"][0]
ENERG_HI = RMF["ENERG_HI"][0]
MIGRA_LO = RMF["MIGRA_LO"][0]
MIGRA_HI = RMF["MIGRA_HI"][0]
Ebins = np.concatenate((MIGRA_LO, [MIGRA_HI[-1]]))

AEFF_matrix = np.outer(np.ones(len(MIGRA_LO)), AEFF_th)
R = AEFF_matrix * RMF_th * Texp

ev = hdul["EVENTS"].data
ev_ra = ev["RA"]
ev_dec = ev["DEC"]
ev_en = ev["ENERGY"]
ev_time = ev["TIME"]
ev_coords = SkyCoord(ev_ra, ev_dec, unit="degree")
sep_s = ev_coords.separation(Coords_s).deg
sep_b = ev_coords.separation(Coords_b).deg

mask = sep_s < R_s
cts_s = np.histogram(ev_en[mask], bins=Ebins)[0]
mask = sep_b < R_s
cts_b = np.histogram(ev_en[mask], bins=Ebins)[0]
hdul.close()
src = cts_s - cts_b

def model(e1, e2, N, slope):
    return N * (e1 / 25.0) ** slope

N = 2e-4
slope = -2

def convolve(p, return_rate=False):
    N, slope = p

    model_rate = np.dot(model(ENERG_LO, ENERG_HI, N, slope), R)

    m = c_e1 > e_min
    m &= c_e2 < e_max

    chi2 = (((model_rate - rate) / rate_err)[m] ** 2).sum()

    print(p, chi2, chi2 / np.sum(m))

    if return_rate:
        return model_rate
    else:
        return chi2

plt.plot(ENERG_LO, model(ENERG_LO, ENERG_HI, N, slope))
plt.xscale("log")
plt.yscale("log")

Aeff_matrix = np.outer(np.ones(160), Aeff(96))
R = Aeff_matrix * MATRIX

OBSlist[0]

cts_s = np.zeros(NEbins)
cts_b = np.zeros(NEbins)
ARF = []
RMF = []
for f in OBSlist:
    cts_s = np.zeros(NEbins)
    cts_b = np.zeros(NEbins)
    Expos = np.zeros(NEbins)
    hdul = fits.open("data/" + f)
    RA_pnt = hdul[1].header["RA_PNT"]
    DEC_pnt = hdul[1].header["DEC_PNT"]
    Texp = hdul[1].header["LIVETIME"]
    dRA = RA - RA_pnt
    dDEC = DEC - DEC_pnt
    RA_b = RA_pnt - dRA
    DEC_b = DEC_pnt - dDEC
    Coords_b = SkyCoord(RA_b, DEC_b, unit="degree")
    Coords_pnt = SkyCoord(RA_pnt, DEC_pnt, unit="degree")
    dist = Coords_pnt.separation(Coords_s).deg

    RMF.append(hdul["EDISP"])
    ARF.append(hdul["AEFF"])

    ev = hdul["EVENTS"].data
    ev_ra = ev["RA"]
    ev_dec = ev["DEC"]
    ev_en = ev["ENERGY"]
    ev_time = ev["TIME"]
    ev_coords = SkyCoord(ev_ra, ev_dec, unit="degree")
    sep_s = ev_coords.separation(Coords_s).deg
    sep_b = ev_coords.separation(Coords_b).deg

    hdu = hdul["AEFF"].data
    EEmin = hdu["ENERG_LO"][0]
    EEmax = hdu["ENERG_HI"][0]
    lgEE = log10(sqrt(EEmin * EEmax))
    lgAA = log10(hdu["EFFAREA"][0] + 1e-10)
    Thmin = hdu["THETA_LO"][0]
    Thmax = hdu["THETA_HI"][0]
    ind = np.argmin((Thmin - dist) ** 2)
    # Expos+=10**(np.interp(lgEmean,lgEE,lgAA[ind]))*Texp
    # ARF.append(10**(np.interp(lgEmean,lgEE,lgAA[ind]))*Texp)
    mask = sep_s < R_s
    cts_s += np.histogram(ev_en[mask], bins=Ebins)[0]
    mask = sep_b < R_s
    cts_b += np.histogram(ev_en[mask], bins=Ebins)[0]
    hdul.close()

BKG

bin_image = PictureProduct.from_file("Spectrum.png")
from astropy.table import Table

data = [Emean, Emin, Emax, flux, flux_err, cts_s, cts_b, Expos * 1e4]
names = (
    "Emean[TeV]",
    "Emin[TeV]",
    "Emax[TeV]",
    "Flux[TeV/cm2s]",
    "Flux_error[TeV/cm2s]",
    "Cts_s",
    "Cts_b",
    "Exposure[cm2s]",
)
spec = ODAAstropyTable(Table(data, names=names))

picture_png = bin_image  # http://odahub.io/ontology#ODAPictureProduct
spectrum_astropy_table = spec  # http://odahub.io/ontology#ODAAstropyTable

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_Spectrum_counts_IRF_picture_png",
        "picture_png_galaxy.output",
        picture_png,
    )
)
_oda_outs.append(
    (
        "out_Spectrum_counts_IRF_spectrum_astropy_table",
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
