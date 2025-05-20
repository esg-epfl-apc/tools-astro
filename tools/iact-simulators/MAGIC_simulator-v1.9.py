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
from scipy.interpolate import RectBivariateSpline, RegularGridInterpolator

src_name = "Crab"  # http://odahub.io/ontology#AstrophysicalObject
RA = 83.628700  # http://odahub.io/ontology#PointOfInterestRA
DEC = 22.014700  # http://odahub.io/ontology#PointOfInterestDEC

timeh = 20  # http://odahub.io/ontology#TimeIntervalHours ; oda:label "Observation time"; oda:descritpion "[h], time of observations"
extension = 0.0  # http://odahub.io/ontology#AngleDegrees ; oda:label "Source extension"

redshift = (
    0.13  # http://odahub.io/ontology#Float ; oda:label "Source redshift"
)

zenith = "lowzd"  #  http://odahub.io/ontology#String ; oda:label "Zenith angle"; oda:allowed_value  'lowzd', 'midzd', 'hizd'
isSUMT = False  #  http://odahub.io/ontology#Boolean ; oda:label "Sum Trigger?"

numoff = 3  # http://odahub.io/ontology#Integer
minev = 10.0  # http://odahub.io/ontology#Float ; oda:label "minimum number of events"
minSBR = 0.05  # http://odahub.io/ontology#Float ; oda:label "minimum ratio of excess to background"
PSF = 0.1  # http://odahub.io/ontology#AngleDegrees
offsetdegrad = 1.0  # http://odahub.io/ontology#Float ; oda:label "degradation factor (for offset >0.4 deg)"
eplotmin = 31  # http://odahub.io/ontology#Energy_GeV ; oda:group "Plotting" ; oda:label "minimal energy" ; oda:lower_limit 30.01 ; oda:upper_limit 29999.
eplotmax = 20.0e3  # http://odahub.io/ontology#Energy_GeV ; oda:group "Plotting" ; oda:label "maximal energy" ; oda:lower_limit 30.01 ; oda:upper_limit 29999.
yplotmin = 1.0e-14  # http://odahub.io/ontology#Float ; oda:group "Plotting" ; oda:label "minimal flux [TeV/(cm2 s)]"
yplotmax = 1.0e-9  # http://odahub.io/ontology#Float ; oda:group "Plotting" ; oda:label "maximal flux [TeV/(cm2 s)]"
minerror = 2  # http://odahub.io/ontology#Float ; oda:group "Plotting" ; oda:label "Minimal errorbar (signal-to-noise)"
drawsigma = True  # http://odahub.io/ontology#Boolean

dN_dE = "2.0e-11*pow(E/1000., -1.99)*exp(-E/1000)"  # http://odahub.io/ontology#String ; oda:label "Source spectrum dN/dE [1/(TeV cm2 s)]"

pulsar_mode = False  # http://odahub.io/ontology#Boolean ; oda:group "Pulsar analysis" ; oda:label "Pulsar analysis?"
on_phase_interval = 0.092  # http://odahub.io/ontology#Float ; oda:group "Pulsar analysis" ; oda:label "range of ON phases"
off_phase_interval = 0.25  # http://odahub.io/ontology#Float ; oda:group "Pulsar analysis" ; oda:label "range of OFF phases"

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
    "zenith",
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
    "pulsar_mode",
    "on_phase_interval",
    "off_phase_interval",
    "isLSTmode",
]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

pulsarmode = pulsar_mode
pulsarOnRange = on_phase_interval
pulsarOffRange = off_phase_interval

# global variables (DO NOT MODIFY)
npoints = 13
pathebl = "dominguez_ebl_tau.txt"  # path with EBL model of Dominguez+11
ismc = False
version = "1.9"

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

def FluxObs(z, xx, fluxint):
    en = xx
    EBLz, EBLen, EBLtaus = LoadEBL()
    # ftau = RegularGridInterpolator((EBLen,EBLz), EBLtaus, method='cubic')
    # tau = ftau([en,z])
    ftau = RectBivariateSpline(
        EBLz, EBLen, EBLtaus.T, kx=3, ky=3, s=0
    )  # cubic, no smoothing
    tau = ftau(z, en)[0]
    atten = np.exp(-tau)
    return fluxint * atten

def Prepare():
    npoints_array = np.arange(0, npoints + 1)
    enbins = 100.0 * pow(10.0, (npoints_array - 2) * 0.2)  # [GeV]
    crabrate = np.zeros(npoints)  # [min^-1]
    bgdrate = np.zeros(npoints)  # [min^-1]
    if isLSTmode:
        if (
            zenith == "hizd"
        ):  # MC values values computed at 59 zd angle, might by ~20% too optimistic >~1 TeV
            crabrate[0] = 0  # 50GeV, below threshold
            crabrate[1] = 0.0  # 79 GeV, below threshol
            crabrate[2] = 0.0  # 126 GeV below threshold
            crabrate[3] = 0.0  # 200 GeV, below threshold
            crabrate[4] = 0.6929  # 316 GeV
            crabrate[5] = 2.0531
            crabrate[6] = 1.5028
            crabrate[7] = 0.8963
            crabrate[8] = 0.5915
            crabrate[9] = 0.3513
            crabrate[10] = 0.2314
            crabrate[11] = 0.111  # 7 TeV
            crabrate[12] = 0.0488  # 12.6 TeV

            bgdrate[0] = 0  #  below threshold
            bgdrate[1] = 0.0
            bgdrate[2] = 0.0
            bgdrate[3] = 0.0  # 200 GeV below threshold
            bgdrate[4] = 0.349
            bgdrate[5] = 0.2163
            bgdrate[6] = 0.0393
            bgdrate[7] = 0.0126
            bgdrate[8] = 0.0037
            bgdrate[9] = 0.0013
            bgdrate[10] = 0.0018
            bgdrate[11] = 0.0005
            bgdrate[12] = 0.0006  # 12.6 TeV
            # the MC sensitivites are often too optimistic. While at low zenith MAGIC+LST1 sensitivity more or less matches
            # between the data and MCs, at mid zenith there is some ~20% difference >~1 TeV. Very likely the same happens
            # at high zenith, so this is why we artificially increase background by 40% to emulate this effect
            # and make more fair MAGIC-only vs MAGIC+LST1 comparisons at high zenith
            bgdrate[np.sqrt(enbins[:-1] * enbins[1:]) > 1000] *= 1.4

        elif (
            zenith == "midzd"
        ):  # values from https://doi.org/10.1051/0004-6361/202346927
            crabrate[0] = 0  # 50GeV, below threshold
            crabrate[1] = 0.59
            crabrate[2] = 2.43
            crabrate[3] = 2.65
            crabrate[4] = 2.03
            crabrate[5] = 1.171
            crabrate[6] = 0.899
            crabrate[7] = 0.806
            crabrate[8] = 0.319
            crabrate[9] = 0.185
            crabrate[10] = 0.113
            crabrate[11] = 0.084
            crabrate[12] = 0  # missing data

            bgdrate[0] = 0  # below threshold
            bgdrate[1] = 0.715
            bgdrate[2] = 1.42
            bgdrate[3] = 0.307
            bgdrate[4] = 0.093
            bgdrate[5] = 0.0229
            bgdrate[6] = 0.0093
            bgdrate[7] = 0.0096
            bgdrate[8] = 0.00264
            bgdrate[9] = 0.0007
            bgdrate[10] = 0.00138
            bgdrate[11] = 0.00148
            bgdrate[12] = 0  # missing data

        elif zenith == "lowzd":
            # those numbers are taken from MC, but they agree
            # with a small bunch of data available
            crabrate[0] = 0.3795  # 50GeV
            crabrate[1] = 2.5808  # 80GeV
            crabrate[2] = 2.8257
            crabrate[3] = 1.6654
            crabrate[4] = 1.3289
            crabrate[5] = 1.2871
            crabrate[6] = 0.8834
            crabrate[7] = 0.7045
            crabrate[8] = 0.3908
            crabrate[9] = 0.1963
            crabrate[10] = 0.089
            crabrate[11] = 0.0368
            crabrate[12] = 0  # missing data

            bgdrate[0] = 8.8870e-01
            bgdrate[1] = 1.6505e00
            bgdrate[2] = 5.7688e-01
            bgdrate[3] = 5.0804e-02
            bgdrate[4] = 2.0521e-02
            bgdrate[5] = 2.1718e-02
            bgdrate[6] = 5.6106e-03
            bgdrate[7] = 3.9491e-03
            bgdrate[8] = 3.2053e-03
            bgdrate[9] = 1.8074e-03
            bgdrate[10] = 7.3253e-04
            bgdrate[11] = 1.1352e-05
            bgdrate[12] = 0  # missing data

    else:
        if isSUMT:
            if zenith == "lowzd":
                # 2018/19 Crab data analyzed with generic ST0307 SUMT MCs
                # 1st point  ~50GeV, last point ~12TeV
                crabrate[0] = 1.39684  # ~50GeV
                crabrate[1] = 3.12657
                crabrate[2] = 3.09145
                crabrate[3] = 2.40268
                crabrate[4] = 1.32915
                crabrate[5] = 0.86180
                crabrate[6] = 0.51666
                crabrate[7] = 0.31533
                crabrate[8] = 0.16207
                crabrate[9] = 0.09279
                crabrate[10] = 0.04624
                crabrate[11] = 0.02345
                crabrate[12] = 0.00874  # ~12TeV

                bgdrate[0] = 3.33321
                bgdrate[1] = 3.24046
                bgdrate[2] = 1.32361
                bgdrate[3] = 0.406588
                bgdrate[4] = 0.091944
                bgdrate[5] = 0.032226
                bgdrate[6] = 0.007277
                bgdrate[7] = 0.003123
                bgdrate[8] = 0.001487
                bgdrate[9] = 0.001464
                bgdrate[10] = 0.001231
                bgdrate[11] = 0.001152
                bgdrate[12] = 0.000957
            else:
                print("Only low zenith sensitivity is available for SUMT :-(")
                # to be implemented
        else:
            if (
                zenith == "hizd"
            ):  # Crab results from 2.5 hrs of 55-62 zd data from 2016-2018. Dataset from Juliane van Scherpenberg.
                crabrate[0] = 0  # 50GeV, below threshold
                crabrate[1] = 0.0  # 79 GeV, below threshol
                crabrate[2] = 0.0  # 126 GeV below threshold
                crabrate[3] = 0.0  # 200 GeV, below threshold
                crabrate[4] = 0.503462  # 316 GeV
                crabrate[5] = 1.60232
                crabrate[6] = 2.26558
                crabrate[7] = 0.928094
                crabrate[8] = 0.698335
                crabrate[9] = 0.305662
                crabrate[10] = 0.173859
                crabrate[11] = 0.083892  # 7 TeV
                crabrate[12] = 0.069938  # 12.6 TeV

                bgdrate[0] = 0  # below threshold
                bgdrate[1] = 0.0
                bgdrate[2] = 0.0
                bgdrate[3] = 0.0  # 200 GeV below threshold
                bgdrate[4] = 0.750815
                bgdrate[5] = 0.597588
                bgdrate[6] = 0.564753
                bgdrate[7] = 0.089775
                bgdrate[8] = 0.0568584
                bgdrate[9] = 0.00954936
                bgdrate[10] = 0.00344762
                bgdrate[11] = 0.00147755
                bgdrate[12] = 0.00229841  # 12.6 TeV
            elif zenith == "midzd":  # Aleksic et al 2016
                crabrate[0] = 0  # below threshold !!!
                crabrate[1] = 0.404836
                crabrate[2] = 3.17608
                crabrate[3] = 2.67108
                crabrate[4] = 2.86307
                crabrate[5] = 1.76124
                crabrate[6] = 1.43988
                crabrate[7] = 0.944385
                crabrate[8] = 0.673335
                crabrate[9] = 0.316263
                crabrate[10] = 0.200331
                crabrate[11] = 0.0991222
                crabrate[12] = 0.0289831

                bgdrate[0] = 1.67777  # below threashold
                bgdrate[1] = 2.91732
                bgdrate[2] = 2.89228
                bgdrate[3] = 0.542563
                bgdrate[4] = 0.30467
                bgdrate[5] = 0.0876449
                bgdrate[6] = 0.0375621
                bgdrate[7] = 0.0197085
                bgdrate[8] = 0.0111295
                bgdrate[9] = 0.00927459
                bgdrate[10] = 0.00417356
                bgdrate[11] = 0.00521696
                bgdrate[12] = 0.000231865
            elif zenith == "lowzd":  # Aleksic et al 2016
                # 0-30 deg values,
                crabrate[0] = 0.818446
                crabrate[1] = 3.01248
                crabrate[2] = 4.29046
                crabrate[3] = 3.3699
                crabrate[4] = 1.36207
                crabrate[5] = 1.21791
                crabrate[6] = 0.880268
                crabrate[7] = 0.579754
                crabrate[8] = 0.299179
                crabrate[9] = 0.166192
                crabrate[10] = 0.0931911
                crabrate[11] = 0.059986
                crabrate[12] = 0.017854

                bgdrate[0] = 3.66424
                bgdrate[1] = 4.05919
                bgdrate[2] = 2.41479
                bgdrate[3] = 0.543629
                bgdrate[4] = 0.0660764
                bgdrate[5] = 0.0270313
                bgdrate[6] = 0.0132653
                bgdrate[7] = 0.00592351
                bgdrate[8] = 0.00266975
                bgdrate[9] = 0.00200231
                bgdrate[10] = 0.00141831
                bgdrate[11] = 0.00458864
                bgdrate[12] = 0.0016686
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

tau = []
ftau2 = RectBivariateSpline(
    EBLz, EBLen, EBLtaus.T, kx=3, ky=3, s=0
)  # cubic, no smoothing
for i in range(len(EBLen)):
    tau.append(ftau2(z, EBLen[i])[0])
plt.plot(EBLen, exp(-np.array(tau)))

plt.xscale("log")
plt.yscale("log")
plt.ylim(1e-3, 2)
max(EBLen)

get_ipython().run_line_magic("pinfo", "RectBivariateSpline")   # noqa: F821

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
    if isSUMT and (zenith != "lowzd"):
        print("Only low zenith sensitivity is available for SUMT :-(")
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
    nexc_all = 0
    noff_all = 0
    best_int_e = -1
    best_int_sigma = -1
    pulsarOnOffRatio = pulsarOnRange / pulsarOffRange

    for i, e1, e2 in reversed(
        list(zip(range(0, len(enbins)), enbins, enbins[1:]))
    ):
        intcrab, error = quad(Crab, e1, e2)
        if redshift > 0:
            intass, error = quad(
                lambda x: FluxObs(redshift, x, Assumed(x)), e1, e2
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
                np.sqrt(noff) * pulsarOnOffRatio
            )  # ignoring numoff for pulsars and scaling for the phase difference
            dexc = np.sqrt(
                nexc + noffon + dnoff * dnoff
            )  # error on the number of excess events
        nexc_all += nexc
        noff_all += noff

        # for tiny excesses (1.e-7 events) the function below can have numerical problems, and either way sigma should be 0 then
        sigma = 0
        if nexc > 0.01:
            if pulsarmode:
                sigma = SignificanceLiMa(nexc + noffon, noff, pulsarOnOffRatio)
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

        sigma_all = -1
        if nexc_all > 0.01:
            if pulsarmode:
                sigma_all = SignificanceLiMa(
                    nexc_all + noff_all * pulsarOnOffRatio,
                    noff_all,
                    pulsarOnOffRatio,
                )
            else:
                sigma_all = SignificanceLiMa(
                    nexc_all + noff_all, noff_all * numoff, 1.0 / numoff
                )
                print(
                    "Integral significance > {0:.1f} GeV = {1:.2f} with {2:.1f} excess".format(
                        enbins[i], sigma_all, nexc_all
                    ),
                    (
                        " "
                        if pulsarmode
                        else ", SBR={0:.2f}".format(nexc_all / noff_all)
                    ),
                )
            if (nexc_all > minev) and (
                pulsarmode or (nexc_all > minSBR * noff_all)
            ):
                if sigma_all > best_int_sigma:
                    best_int_sigma = sigma_all
                    best_int_e = enbins[i]

        # print(nexc, minerror, dexc)
        if nexc > minerror * dexc:
            tmpen = np.sqrt(enbins[i] * enbins[i + 1])
            en.append(tmpen)
            if redshift > 0:
                tmpsed = CalcSED(
                    tmpen, FluxObs(redshift, tmpen, Assumed(tmpen))[0]
                )
            else:
                tmpsed = CalcSED(tmpen, Assumed(tmpen))
            sed.append(float(tmpsed))
            dsed.append(float(tmpsed * dexc / nexc))
            sigmas.append(sigma)
            detected.append(detect)

    inst = "MAGIC+LST1" if isLSTmode else "MAGIC"

    # preparing reference SED graph
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set(
        xlim=(eplotmin, eplotmax),
        ylim=(yplotmin, yplotmax),
        xlabel="E [GeV]",
        ylabel="E$^2$ dN/dE [TeV cm$^{-2}$ s$^{-1}$]",
    )
    x = np.logspace(np.log10(eplotmin), np.log10(eplotmax), 50)
    labeltext = "expected SED ($T_{obs}$ = " + str(timeh) + " h)"
    plt.plot(
        x, 1.0e-6 * x * x * Crab(x), "0.3", label="Crab (Aleksic et al 2016)"
    )
    if redshift > 0:
        plt.plot(
            x,
            1.0e-6 * x * x * FluxObs(redshift, x, Assumed(x)),
            "limegreen",
            label=src_name + " (Assumed, z={0:.2f})".format(redshift),
        )
    else:
        plt.plot(
            x,
            1.0e-6 * x * x * Assumed(x),
            "limegreen",
            label=src_name + " (Assumed)",
        )
    if len(en) > 0:
        ax.errorbar(en, sed, yerr=dsed, label=labeltext, color="0", fmt="o")
        handles, labels = ax.get_legend_handles_labels()
    else:
        handles, labels = ax.get_legend_handles_labels()
    #    scatter_proxy = mlines.Line2D([], [] , color='0', marker='o', linestyle='None', markersize=3)
    #    handles.append(scatter_proxy)
    #    labels.append(labeltext)

    # zdtag = "lowZd";
    # if (zenith == Zeniths.midzd): zdtag="midZd"
    # if (zenith == Zeniths.hizd): zdtag="highZd"
    zdtag = zenith

    header = ""
    header = inst + " " + zdtag + (", SUMT " if isSUMT else "")

    ax.legend(handles=handles, labels=labels, loc="upper right", title=header)
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
    ax.annotate(
        "mss v"
        + version
        + (", (MC based)" if ismc else "")
        + (f", offset degr.={offsetdegrad}" if offsetdegrad < 0.99 else ""),
        xy=(0.02, 0.02),
        xycoords="axes fraction",
    )
    if drawsigma and (best_int_sigma > 0):
        ax.set(
            title="$\sigma$ (>{0:.1f} GeV) = {1:.2f}".format(
                best_int_e, best_int_sigma
            )
        )
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
        "out_MAGIC_simulator_v1.9_spectrum_png",
        "spectrum_png_galaxy.output",
        spectrum_png,
    )
)
_oda_outs.append(
    (
        "out_MAGIC_simulator_v1.9_spectrum_table",
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
