#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import copy
import json
import os
import re
import shutil

import matplotlib.pyplot as plt
import numpy as np
from numpy import exp
from oda_api.json import CustomJSONEncoder
from scipy.integrate import quad
from scipy.interpolate import RegularGridInterpolator

src_name = "Crab"  # http://odahub.io/ontology#AstrophysicalObject
RA = 83.628700  # http://odahub.io/ontology#PointOfInterestRA
DEC = 22.014700  # http://odahub.io/ontology#PointOfInterestDEC

timeh = 20  # http://odahub.io/ontology#TimeIntervalHours ; oda:label "Observation time"; oda:descritpion "[h], time of observations"
extension = 0.0  # http://odahub.io/ontology#AngleDegrees

redshift = 0.13  # http://odahub.io/ontology#Double

ismidzd = False  #  http://odahub.io/ontology#Boolean
isSUMT = False  #  http://odahub.io/ontology#Boolean

numoff = 3  # http://odahub.io/ontology#Integer
minev = 10.0  # http://odahub.io/ontology#Double
minSBR = 0.05  # http://odahub.io/ontology#Double
PSF = 0.1  # http://odahub.io/ontology#AngleDegrees
offsetdegrad = 1.0  # http://odahub.io/ontology#Double
eplotmin = 31  # http://odahub.io/ontology#Energy_GeV ; oda:lower_limit 30.01 ; oda:upper_limit 29999.
eplotmax = 20.0e3  # http://odahub.io/ontology#Energy_GeV ; oda:lower_limit 30.01 ; oda:upper_limit 29999.
yplotmin = 1.0e-14  # http://odahub.io/ontology#Double
yplotmax = 1.0e-9  # http://odahub.io/ontology#Double
minerror = 2  # http://odahub.io/ontology#Double
drawsigma = True  # http://odahub.io/ontology#Boolean
dN_dE = "2.0e-11*pow(E/1000., -1.99)*exp(-E/1000)"  # http://odahub.io/ontology#String

pulsarmode = (
    False  # http://odahub.io/ontology#Boolean ; oda:group "Pulsar analysis"
)
pulsarOnRange = (
    0.092  # http://odahub.io/ontology#Double ; oda:group "Pulsar analysis"
)
pulsarOffRange = (
    0.25  # http://odahub.io/ontology#Double ; oda:group "Pulsar analysis"
)

isLSTmode = False  # http://odahub.io/ontology#Boolean

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
    "timeh",
    "extension",
    "redshift",
    "ismidzd",
    "isSUMT",
    "numoff",
    "minev",
    "minSBR",
    "PSF",
    "offsetdegrad",
    "eplotmin",
    "eplotmax",
    "yplotmin",
    "yplotmax",
    "minerror",
    "drawsigma",
    "dN_dE",
    "pulsarmode",
    "pulsarOnRange",
    "pulsarOffRange",
    "isLSTmode",
]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

# global variables (DO NOT MODIFY)
npoints = 13
pathebl = "dominguez_ebl_tau.txt"  # path with EBL model of Dominguez+11

version = "1.7"

def parse_spectrum(input_str):
    dnde_str = copy.copy(input_str)

    numbers = re.findall(r"-?\d+\.?\d*(?:[eE][+-]?\d+)?", input_str)
    funcs = re.findall(r"\b\w+(?=\()", input_str)

    for num in numbers:
        input_str = input_str.replace(num, "")

    for func in funcs:
        input_str = input_str.replace(func, "")

    service_symbols = "E()[]*/.,+- "

    for symbol in service_symbols:
        input_str = input_str.replace(symbol, "")

    if input_str:
        raise ValueError("forbidden statements")

    return eval(f"lambda E: {dnde_str}")

# test_str = "2.0e-11*pow(E/1000., -1.99)*exp(-E/100); !wget https://scripts.com/myscript.sh"
Assumed = parse_spectrum(dN_dE)

# Crab Nebula [Aleksic et al. 2016]
def Crab(x):
    return 3.39e-11 * pow(x / 1000.0, -2.51 - 0.21 * np.log10(x / 1000.0))

# const TString crabstr="3.39e-11 * pow(x/1000, -2.51-0.21*log10(x/1000))"; //Crab Nebula [Aleksic et al. 2016]

# Calculate SED value from energy and flux, for display
def CalcSED(x, flux):
    return 1.0e-6 * x * x * flux

def LoadEBL():
    f1 = open(pathebl, "r")
    line = f1.readline()
    firstlineEBL = list(map(float, line.split()))
    firstlineEBL[0]
    zz = firstlineEBL[1:]
    # print(nz, zz)
    f1.close()
    eblbody = np.loadtxt(pathebl, delimiter=" ", skiprows=1)
    energies = eblbody[:, 0] * 1.0e3  # in GeV  #en*=1.e3; // GeV
    taus = eblbody[:, 1:]
    if len(taus > 0):
        return zz, energies, taus

def FluxObs(xx, z, fluxint):
    en = xx
    EBLz, EBLen, EBLtaus = LoadEBL()
    # ftau = interpolate.interp2d(EBLz, EBLen, EBLtaus, kind='cubic')
    ftau = RegularGridInterpolator((EBLen, EBLz), EBLtaus, method="cubic")
    tau = ftau([en, z])
    atten = np.exp(-tau)
    return fluxint * atten

def Prepare():
    npoints_array = np.arange(0, npoints + 1)
    enbins = 100.0 * pow(10.0, (npoints_array - 2) * 0.2)  # [GeV]

    if isLSTmode:
        if ismidzd:  # values from https://doi.org/10.1051/0004-6361/202346927
            # 1st point is 50 GeV, below threshold; the last point is missing data
            crabrate = np.array(
                [
                    0,
                    0.59,
                    2.43,
                    2.65,
                    2.03,
                    1.171,
                    0.899,
                    0.806,
                    0.319,
                    0.185,
                    0.113,
                    0.084,
                    0,
                ]
            )
            bgdrate = np.array(
                [
                    0,
                    0.715,
                    1.42,
                    0.307,
                    0.093,
                    0.0229,
                    0.0093,
                    0.0096,
                    0.00264,
                    0.0007,
                    0.00138,
                    0.00148,
                    0,
                ]
            )
        else:
            # 1st point is 50 GeV; 2nd point 80 GeV; the last point is missing data
            # those numbers are taken from MC, but they agree
            # with a small bunch of data available
            crabrate = np.array(
                [
                    0.3795,
                    2.5808,
                    2.8257,
                    1.6654,
                    1.3289,
                    1.2871,
                    0.8834,
                    0.7045,
                    0.3908,
                    0.1963,
                    0.089,
                    0.0368,
                    0,
                ]
            )
            bgdrate = np.array(
                [
                    8.8870e-01,
                    1.6505e00,
                    5.7688e-01,
                    5.0804e-02,
                    2.0521e-02,
                    2.1718e-02,
                    5.6106e-03,
                    3.9491e-03,
                    3.2053e-03,
                    1.8074e-03,
                    7.3253e-04,
                    1.1352e-05,
                    0,
                ]
            )
    else:
        if isSUMT:
            if ismidzd:
                crabrate = np.array([])
                bgdrate = np.array([])
                # to be implemented
            else:
                # 2018/19 Crab data analyzed with generic ST0307 SUMT MCs
                # 1st point  ~50GeV, last point ~12TeV
                crabrate = np.array(
                    [
                        1.39684,
                        3.12657,
                        3.09145,
                        2.40268,
                        1.32915,
                        0.86180,
                        0.51666,
                        0.31533,
                        0.16207,
                        0.09279,
                        0.04624,
                        0.02345,
                        0.00874,
                    ]
                )  # [min^-1]
                bgdrate = np.array(
                    [
                        3.33321,
                        3.24046,
                        1.32361,
                        0.406588,
                        0.091944,
                        0.032226,
                        0.007277,
                        0.003123,
                        0.001487,
                        0.001464,
                        0.001231,
                        0.001152,
                        0.000957,
                    ]
                )  # [min^-1]
        else:
            if ismidzd:  # Aleksic et al 2016
                crabrate = np.array(
                    [
                        0,
                        0.404836,
                        3.17608,
                        2.67108,
                        2.86307,
                        1.76124,
                        1.43988,
                        0.944385,
                        0.673335,
                        0.316263,
                        0.200331,
                        0.0991222,
                        0.0289831,
                    ]
                )  # [min^-1], 1st point below threshold !!!
                bgdrate = np.array(
                    [
                        1.67777,
                        2.91732,
                        2.89228,
                        0.542563,
                        0.30467,
                        0.0876449,
                        0.0375621,
                        0.0197085,
                        0.0111295,
                        0.00927459,
                        0.00417356,
                        0.00521696,
                        0.000231865,
                    ]
                )  # [min^-1], 1st point below threshold !!!
            else:  # 0-30 deg values, Aleksic et al 2016
                crabrate = np.array(
                    [
                        0.818446,
                        3.01248,
                        4.29046,
                        3.3699,
                        1.36207,
                        1.21791,
                        0.880268,
                        0.579754,
                        0.299179,
                        0.166192,
                        0.0931911,
                        0.059986,
                        0.017854,
                    ]
                )  # [min^-1]
                bgdrate = np.array(
                    [
                        3.66424,
                        4.05919,
                        2.41479,
                        0.543629,
                        0.0660764,
                        0.0270313,
                        0.0132653,
                        0.00592351,
                        0.00266975,
                        0.00200231,
                        0.00141831,
                        0.00458864,
                        0.0016686,
                    ]
                )  # [min^-1]
        # offset degradation is only applied if NOT in LST mode
        crabrate *= offsetdegrad
        bgdrate *= offsetdegrad

    return npoints_array, enbins, crabrate, bgdrate

EBLz, EBLen, EBLtaus = LoadEBL()
ftau = RegularGridInterpolator((EBLen, EBLz), EBLtaus, method="cubic")
z = 0.13
tau = []
for i in range(len(EBLen)):
    tau.append(ftau([EBLen[i], z]))
plt.plot(EBLen, exp(-np.array(tau)))
z = 0.18
tau = []
for i in range(len(EBLen)):
    tau.append(ftau([EBLen[i], z]))
plt.plot(EBLen, exp(-np.array(tau)))
plt.xscale("log")
plt.yscale("log")
plt.ylim(1e-3, 2)
max(EBLen)

# // Li & Ma eq 17. significance, function abridged from MARS
def SignificanceLiMa(s, b, alpha):
    if b == 0 and s == 0:
        return 0

    b = (
        FLT_MIN if b == 0 else b
    )  # // Guarantee that a value is output even in the case b or s == 0
    s = (
        FLT_MIN if s == 0 else s
    )  # //   (which is not less allowed, possible, or meaningful than
    # //    doing it with b,s = 0.001)
    sumsb = s + b

    if sumsb < 0 or alpha <= 0:
        return -1
    l = s * np.log(s / sumsb * (alpha + 1) / alpha)
    m = b * np.log(b / sumsb * (alpha + 1))

    return -1 if l + m < 0 else np.sqrt((l + m) * 2)

def Checks():
    if extension > 1:
        raise RuntimeError(
            "Extension comparable to the size of the MAGIC camera cannot be simulated"
        )
        return False
    if (numoff <= 0) or (numoff > 7):
        raise RuntimeError(
            "Number of OFF estimation regions must be in the range 1-7"
        )
        return False
    if (extension > 0.5) and (numoff > 1):
        print(
            "For large source extensions 1 OFF estimation region (numoff) should be used"
        )
        raise RuntimeError(
            "Number of OFF estimation regions must be in the range 1-7"
        )
        return False
    if isSUMT and ismidzd:
        print("Sorry not implemented yet :-(")
        return False
    if isLSTmode and isSUMT:
        print("LST mode is not compatible with SUMT")
        return False
    if offsetdegrad > 1.00001:
        print(
            "No cheating! the performance degradation ({0}) should not be larger then 1".format(
                offsetdegrad
            )
        )
        return False
    if pulsarmode:
        if pulsarOnRange <= 0 or pulsarOnRange >= 1:
            print(
                "Pulsar mode ON phase range is {0}, and it should be in range (0,1)".format(
                    pulsarOnRange
                )
            )
            return False
        if pulsarOffRange <= 0 or pulsarOffRange >= 1:
            print(
                "Pulsar mode OFF phase range is {0}, and it should be in range (0,1)".format(
                    pulsarOffRange
                )
            )
            return False
        if redshift > 0:
            print(
                "Do you really want to observe a pulsar at redshift of {0} ??".format(
                    redshift
                )
            )
    return True

## MAIN code mss.cpp
if not Checks():
    print("exiting")
# return 0
else:
    npoints_array, enbins, crabrate, bgdrate = Prepare()

    en = []
    sed = []
    dsed = []
    sigmas = []
    detected = []
    for i, e1, e2 in zip(range(0, len(enbins)), enbins, enbins[1:]):
        intcrab, error = quad(Crab, e1, e2)
        if redshift > 0:
            intass, error = quad(
                lambda x: FluxObs(x, redshift, Assumed(x)), e1, e2
            )
        else:
            intass, error = quad(Assumed, e1, e2)
        noff = bgdrate[i] * timeh * 60  # number of off events
        noff *= (PSF * PSF + extension * extension) / (
            PSF * PSF
        )  # larger integration cut due to extension
        dnoff = np.sqrt(
            noff / numoff
        )  # error on the number of off events (computed from numoff regions)
        nexc = crabrate[i] * timeh * 60 * intass / intcrab
        dexc = np.sqrt(
            nexc + noff + dnoff * dnoff
        )  # error on the number of excess events

        noffon = 0
        if pulsarmode:
            noffon = noff * pulsarOnRange  # number of bgd events in ON phase
            noff *= pulsarOffRange  # number of bgd events in OFF phase
            dnoff = (
                np.sqrt(noff) * pulsarOnRange / pulsarOffRange
            )  # ignoring numoff for pulsars and scaling for the phase difference
            dexc = np.sqrt(
                nexc + noffon + dnoff * dnoff
            )  # error on the number of excess events

        # for tiny excesses (1.e-7 events) the function below can have numerical problems, and either way sigma should be 0 then
        sigma = 0
        if nexc > 0.01:
            if pulsarmode:
                sigma = SignificanceLiMa(
                    nexc + noffon, noff, pulsarOnRange / pulsarOffRange
                )
                noff = noffon  # needed later for SBR
            else:
                sigma = SignificanceLiMa(
                    nexc + noff, noff * numoff, 1.0 / numoff
                )

        detect = False
        if pulsarmode:
            if (sigma >= 5.0) and (nexc > minev):
                detect = True
        else:
            if sigma >= 5.0 and nexc / noff > minSBR and nexc > minev:
                detect = True
        print(
            "{0:.1f}-{1:.1f} GeV: exc. = {2:.1f}+-{3:.1f} ev., SBR={4:.2f}%, sigma = {5:.1f}".format(
                enbins[i],
                enbins[i + 1],
                nexc,
                dexc,
                100.0 * nexc / noff if noff > 0 else np.nan,
                sigma,
            ),
            " DETECTION" if detect else " ",
        )

        # print(nexc, minerror, dexc)
        if nexc > minerror * dexc:
            tmpen = np.sqrt(enbins[i] * enbins[i + 1])
            en.append(tmpen)
            if redshift > 0:
                tmpsed = CalcSED(
                    tmpen, FluxObs(tmpen, redshift, Assumed(tmpen))
                )
            else:
                tmpsed = CalcSED(tmpen, Assumed(tmpen))
            sed.append(float(tmpsed))
            dsed.append(float(tmpsed * dexc / nexc))
            sigmas.append(sigma)
            detected.append(detect)

    if len(sigmas) > 0:
        sig_comb = sum(sigmas) / np.sqrt(len(sigmas))

    print(
        "Combined significance (using the {0:d} data points shown in the SED) = {1:.1f}".format(
            len(sigmas), sig_comb
        )
    )
    inst = "MAGIC+LST1" if isLSTmode else "MAGIC"
    if sig_comb < 4:
        print(
            "The source probably will not be detected by",
            inst,
            "in ",
            timeh,
            " hours",
        )
    elif sig_comb < 6:
        print(
            "The source probably might be detected by",
            inst,
            "in ",
            timeh,
            " hours",
        )
    else:
        print(
            "The source probably will be detected by",
            inst,
            "in ",
            timeh,
            " hours",
        )

    # preparing reference SED graph
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set(
        xlim=(eplotmin, eplotmax),
        ylim=(yplotmin, yplotmax),
        xlabel="E [GeV]",
        ylabel="E$^2$ dN/dE [TeV cm$^{-2}$ s$^{-1}$]",
        title="$\Sigma\sigma_i / \sqrt{"
        + str(len(sigmas))
        + "} ="
        + str(round(sig_comb, 1))
        + "\sigma$",
    )
    x = np.logspace(np.log10(eplotmin), np.log10(eplotmax), 50)
    labeltext = "expected SED ($T_{obs}$ = " + str(timeh) + " h)"
    plt.plot(
        x, 1.0e-6 * x * x * Crab(x), "0.3", label="Crab (Aleksic et al 2016)"
    )
    if redshift > 0:
        plt.plot(
            x,
            1.0e-6
            * x
            * x
            * [float(FluxObs(ee, redshift, Assumed(ee))) for ee in x],
            "g",
            label="Assumed spectrum, z={0:.2f}".format(redshift),
        )
    else:
        plt.plot(
            x, 1.0e-6 * x * x * Assumed(x), "g", label=" Assumed spectrum"
        )
    if len(en) > 0:
        ax.errorbar(en, sed, yerr=dsed, label=labeltext, color="0", fmt="o")

    header = ""
    if isLSTmode and ismidzd:
        header = inst + " (data), mid zd"
    elif isLSTmode and not (ismidzd):
        header = inst + " (MC), low zd"
    else:
        header = inst + (", midZd " if ismidzd else ", lowZd ")
        header += ", SUMT " if isSUMT else ""

    ax.legend(loc="upper right", title=header)
    ax.grid(True, which="both", axis="x", ls="--", color="0.95")
    ax.grid(True, which="major", axis="y", ls="--", color="0.95")
    if drawsigma:
        for i in range(0, len(sigmas)):
            col = "0" if detected[i] else "0.75"
            ax.text(
                en[i],
                2 * sed[i],
                "{0:.1f}$\sigma$".format(sigmas[i]),
                rotation=90,
                size=10,
                ha="left",
                va="bottom",
                color=col,
            )
    ax.annotate("mss v" + version, xy=(0.02, 0.02), xycoords="axes fraction")
    fig.savefig("Spectrum.png", format="png", bbox_inches="tight")

from oda_api.data_products import PictureProduct

bin_image = PictureProduct.from_file("Spectrum.png")

from astropy.table import Table
from oda_api.data_products import ODAAstropyTable

data = [en, sed, dsed, sigmas]
names = (
    "E[GeV]",
    "E2dNdE[TeV/cm2s]",
    "E2dNdE_error[TeV/cm2s]",
    "detection_sigma",
)
spec = ODAAstropyTable(Table(data, names=names))

spectrum_png = bin_image  # http://odahub.io/ontology#ODAPictureProduct
spectrum_table = spec  # http://odahub.io/ontology#ODAAstropyTable

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_MAGIC_simulator_spectrum_png",
        "spectrum_png_galaxy.output",
        spectrum_png,
    )
)
_oda_outs.append(
    (
        "out_MAGIC_simulator_spectrum_table",
        "spectrum_table_galaxy.output",
        spectrum_table,
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
