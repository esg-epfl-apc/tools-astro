#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import matplotlib.pyplot as plt

# Loading necessary packages
import numpy as np
from oda_api.json import CustomJSONEncoder

get_ipython().run_line_magic("matplotlib", "inline")   # noqa: F821
from random import random

import astropy.constants as const
import astropy.units as u
from astropy.constants import c
from numpy import log, pi, sqrt
from oda_api.data_products import ODAAstropyTable, PictureProduct

# Parameters of the phase transition
# Temperature in GeV
T_star = 0.178  # http://odahub.io/ontology#Energy_GeV

# Numbers of relativistic degrees of freedom
g_star = 20  # http://odahub.io/ontology#Integer

# ratio of the energy density deposited in the bubbles to the radiation energy density
alpha = 1.0  # http://odahub.io/ontology#Float

# beta/H : rate of the phase transition compared to Hubble rate
beta_H = 3.3

# fraction of turbulent energy that goes to gw N.B. arXiv:1004.4187 claims that epsilon_turb=0.05, but checks below show that it is rather 0.01
epsilon_turb = 1  # http://odahub.io/ontology#Float

# terminal velocity of bubbles
v_w = 0.999  # http://odahub.io/ontology#Float
h = 0.7  # http://odahub.io/ontology#Float

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

c_s = 3 ** (-0.5)  # speed of sound
F0_GW = 1.64e-5 / h**2 * (100 / g_star) ** (1 / 3.0)

# Eq.20 of arXiv:1512.06239, corrected according to Appendix A of arXiv:1004.4187

def ka_v(alpha_v, v):
    zeta = ((2.0 / 3.0 * alpha_v + alpha_v**2) ** 0.5 + (1 / 3.0) ** 0.5) / (
        1 + alpha_v
    )
    kappa_a = (
        6.9
        * v ** (6.0 / 5.0)
        * alpha_v
        / (1.36 - 0.037 * alpha_v**0.5 + alpha_v)
    )
    kappa_b = alpha_v ** (2.0 / 5.0) / (
        0.017 + (0.997 + alpha_v) ** (2.0 / 5.0)
    )
    kappa_c = alpha_v**0.5 / (0.135 + (0.98 + alpha_v) ** 0.5)
    kappa_d = alpha_v / (0.73 + 0.083 * alpha_v**0.6 + alpha_v)
    deltak = -0.9 * log(alpha_v**0.5 / (1 + alpha_v**0.5))
    if v < c_s:
        return (
            c_s ** (11.0 / 5.0)
            * kappa_a
            * kappa_b
            / (
                (c_s ** (11.0 / 5.0) - v ** (11.0 / 5.0)) * kappa_b
                + v * c_s ** (6.0 / 5.0) * kappa_a
            )
        )
    elif v > zeta:
        return (
            (zeta - 1) ** 3
            * (zeta / v) ** (5.0 / 2)
            * kappa_c
            * kappa_d
            / (
                ((zeta - 1) ** 3 - (v - 1) ** 3)
                * zeta ** (5.0 / 2.0)
                * kappa_c
                + (v - 1) ** 3 * kappa_d
            )
        )
    else:
        return (
            kappa_b
            + (v - c_s) * deltak
            + (v - c_s) ** 3
            / (zeta - c_s) ** 3
            * (kappa_c - kappa_b - (zeta - c_s) * deltak)
        )

kappa_v = ka_v(alpha, v_w)
kappa_therm = 1 - kappa_v
print(
    "Fraction of released energy in bulk motion:",
    kappa_v,
    ", Thermal energy fraction:",
    kappa_therm,
)
Omega_star = (kappa_v) * alpha / (1 + alpha)
print("Omega in bulk motion:", Omega_star)

# Comoving Hubble rate at the phase transition Eq. 11 of arXiv:1512.06239
def H_star(T_star):
    return (
        16.5e-6
        * (T_star / (100.0 * u.GeV))
        * (g_star / 100) ** (1.0 / 6.0)
        * u.Hz
    )

Hstar = H_star(T_star * u.GeV)
logHstar = np.log10(Hstar / u.Hz)

fmin = logHstar.value - 5
fmax = logHstar.value + 5
ff = np.logspace(fmin, fmax, 101)

# HL model formula
def GW_sound(f, T_star, alpha, beta_H, v_w):
    Hstar = H_star(T_star)
    kappa_v = ka_v(alpha, v_w)
    K = kappa_v * alpha / (1 + alpha)
    lambda_star = (8 * pi) ** (1 / 3) * max([v_w, c_s]) / (beta_H * Hstar) * c
    Delta_w = sqrt((v_w - c_s) ** 2) / v_w
    s2 = Delta_w * lambda_star * f / c
    # print((c/lambda_star).cgs,(c/(Delta_w*lambda_star)).cgs)
    s1 = lambda_star * f / c
    M = (
        16
        * (1 + Delta_w ** (-3)) ** (2 / 3.0)
        * (Delta_w * s1) ** 3
        / ((1 + s1**3) ** (2.0 / 3.0) * (3 + s2**2) ** 2)
    )
    factor = (K * lambda_star * Hstar) ** 2 / (
        sqrt(K) + lambda_star * Hstar / c
    )
    mu = 4.78 - 6.27 * Delta_w + 3.34 * Delta_w**2
    ff = f / u.Hz
    dlnf = (ff[1] - ff[0]) / ff[0]
    mu = sum(M) * dlnf
    B = 1e-2 / mu
    return (3 * B * factor / c**2 * F0_GW * M).cgs

# Eq 1 of the new Overleaf, from Alberto

def GW_turb_Theo(f, T_star, alpha, beta_H, v_w, epsilon_turb):
    (const.c / H_star(T_star)).to(u.Mpc)
    kappa_v = ka_v(alpha, v_w)
    Omega_star = (kappa_v) * alpha / (1 + alpha)
    lambda_star = (
        (8 * pi) ** (1 / 3) * max([v_w, c_s]) / (beta_H * H_star(T_star))
    ) * c  # characteristic light-crossing distance scale

    Delta_w = (
        sqrt((v_w - c_s) ** 2) / v_w
    )  # <1, the width of the bubble wall in units of R_star
    R_peak = (
        Delta_w * lambda_star
    )  # the GW spectrum peaks ad the frequency correspoding to the bubble wall width
    u_star = sqrt(0.75 * epsilon_turb * Omega_star)  # alfven velocity
    dt_fin = 2 * lambda_star / u_star / c  # eddy processing time

    T_GW = np.log10(1 + (H_star(T_star) * dt_fin / (2 * pi))) * (
        f < 1 / dt_fin
    ) + np.log10(1 + (H_star(T_star) / (2 * pi * f))) * (f >= 1 / dt_fin)
    T_GW = T_GW**2
    alpha_pi = 2.15
    P_pi = (1 + ((lambda_star * f) / (2.2 * c)) ** alpha_pi) ** (
        -11 / (3 * alpha_pi)
    )
    M = (
        (lambda_star * f / c) ** 3
        * (4 * pi**2 * T_GW * P_pi)
        / (0.84 * (lambda_star * H_star(T_star) / c) ** 2)
    )
    res = (
        3
        * (1.75 * 10 ** (-3))
        * (Omega_star * epsilon_turb) ** 2
        * (lambda_star * H_star(T_star) / c) ** 2
        * 1.64
        * 10 ** (-5)
        * (100 / g_star) ** (1 / 3.0)
        / h**2
    )
    return res * M

def GW_turb_Andrii(f, T_star, alpha, beta_H, v_w, epsilon_turb):
    Hstar = H_star(T_star)

    # Eq. 1 of 2307.10744
    lambda_star = (
        (8 * pi) ** (1 / 3) * max([v_w, c_s]) / (beta_H * Hstar)
    ) * c  # characteristic light-crossing distance scale

    kappa_v = ka_v(alpha, v_w)
    K = kappa_v * alpha / (1 + alpha)

    # Eq. 10 of 2307.10744
    Omega_star = epsilon_turb * K

    # Eq. 13 of 2307.10744 and text after it
    u_star = sqrt(0.75 * Omega_star)
    dt_fin = 2 * lambda_star / u_star / c
    s3 = dt_fin * f
    s1 = f * lambda_star / c
    # print('u_star=',u_star,'lambda_star=',lambda_star.cgs)

    # Eq. 15 of 2307.10744
    T_GW = np.log(1 + Hstar * dt_fin / (2 * pi)) * (s3 < 1) + np.log(
        1 + lambda_star * Hstar / c / (2 * pi * s1)
    ) * (s3 >= 1)
    T_GW = T_GW**2

    # Eq. 17 of 2307.10744
    alpha_pi = 2.15
    s_pi = 2.2
    P_pi = (1 + (s1 / s_pi) ** alpha_pi) ** (-11 / (3 * alpha_pi))

    # Eq 18,19,20 of 2307.10744
    A = 2e-3 * 1.4 * 0.6

    # Eq. 14 of 2307.10744
    Sturb = (
        4
        * pi**2
        * s1**3
        * T_GW
        / (lambda_star * Hstar / c) ** 2
        * P_pi
        / 1.4
        / 0.6
    )

    # Eq 9 of 2307.10744
    res = (
        3 * A * Omega_star**2 * (lambda_star * Hstar / c) ** 2 * F0_GW * Sturb
    )
    return res

H0 = 70 * (u.km / u.s) / u.Mpc

d = np.genfromtxt("NANOGrav23.csv")
gammas_nano = d[:, 0]
As_nano = 10 ** d[:, 1]
ps_nano = 5 - gammas_nano
plt.plot(gammas_nano, As_nano)

d = np.genfromtxt("EPTA.csv")
gammas_epta = d[:, 0]
As_epta = 10 ** d[:, 1]
ps = 5 - gammas_epta
plt.plot(gammas_epta, As_epta)

plt.yscale("log")

GW_t = GW_turb_Theo(
    ff * u.Hz, T_star * u.GeV, alpha, beta_H, v_w, epsilon_turb
)
# plt.plot(ff,GW_t,color='magenta')
GW_t = GW_turb_Andrii(
    ff * u.Hz, T_star * u.GeV, alpha, beta_H, v_w, epsilon_turb
)
plt.plot(ff, GW_t, color="blue", linestyle="dashed", label="turbulence")
GW_s = GW_sound(ff * u.Hz, T_star * u.GeV, alpha, beta_H, v_w)
plt.plot(ff, GW_s, color="red", linestyle="dotted", label="sound waves")
GW = GW_s + GW_t
plt.plot(ff, GW, linewidth=4, color="black", alpha=0.5, label="total")

maxGW = max(GW)
ind = np.argmax(GW)
fmax = ff[ind]
plt.xscale("log")
plt.yscale("log")
plt.ylim(maxGW / 1e5, maxGW * 10)
plt.xlim(fmax / 1e3, fmax * 1e3)
plt.grid()
butterfly2 = np.array(
    [
        [8.533691739758777e-8, 0.0000016780550400704105],
        [8.361136018023004e-8, 1.0740085826829729e-7],
        [8.499066585218291e-9, 3.882820258009008e-9],
        [1.029481229115365e-9, 3.6152599901104956e-11],
        [1.0353166841868053e-9, 2.766052794917924e-10],
        [7.3944540651210095e-9, 6.632930188159042e-9],
    ]
)
plt.fill(butterfly2[:, 0], butterfly2[:, 1], label="NANOGrav'23")
plt.xscale("log")
plt.yscale("log")
plt.axvline(fref.cgs.value)
plt.legend(loc="upper left")
plt.xlabel("$f$, Hz")
plt.ylabel("$\Omega_{gw}(f)$")
plt.savefig("Spectrum.png", format="png", bbox_inches="tight")

ind = len(ff[ff < fref.cgs.value])
p_model = (GW[ind + 1] - GW[ind - 1]) / (ff[ind + 1] - ff[ind - 1]) / GW[
    ind
] * ff[ind] + 1
GW_model = GW[ind]
A_model = sqrt(GW_model / (2 * pi**2 / 3 / H0**2 * fref**2).cgs.value)
g_model = 5 - p_model
print(A_model, p_model)

plt.plot(gammas_nano, As_nano)
plt.scatter([g_model], [A_model])
plt.yscale("log")

# spec=(2*pi**2/3/H0**2*ff**2*As[i]**2).cgs.value

lgfmin = np.log10(fref.cgs.value / 10)
lgfmax = np.log10(fref.cgs.value / 2)

fff = np.logspace(lgfmin, lgfmax, 10) * u.Hz
GW_interp = np.interp(fff.value, ff, GW)
# plt.plot(fff,GW_interp)
fref = 1 / u.yr
print(fref.cgs)
spec_min = np.ones(len(fff))
spec_max = np.zeros(len(fff))
for i in range(len(As)):
    spec = (
        2
        * pi**2
        / 3
        / H0**2
        * fff**2
        * As[i] ** 2
        * (fff / fref) ** (3 - gammas[i])
    ).cgs.value
    spec_min = np.minimum(spec, spec_min)
    spec_max = np.maximum(spec, spec_max)
    # plt.plot(ff,spec)
plt.plot(fff, spec_min, linewidth=4)
plt.plot(fff, spec_max, linewidth=4)

Tgrid = np.logspace(-3, 1, 31)
agrid = np.logspace(-1, 2, 31)
bgrid = np.logspace(0, 2, 31)
Tgood = []
agood = []
bgood = []
for T in Tgrid:
    print(T)
    for a in agrid:
        for b in bgrid:
            GW_t = GW_turb_Andrii(
                ff * u.Hz, T * u.GeV, a, b, v_w, epsilon_turb
            )
            GW_s = GW_sound(ff * u.Hz, T * u.GeV, a, b, v_w)
            GW = GW_s + GW_t
            GW_interp = np.interp(fff.value, ff, GW)
            if sum((GW_interp < spec_max) * (GW_interp > spec_min)) == len(
                fff
            ):
                agood.append(a)
                bgood.append(b)
                Tgood.append(T)
                plt.plot(fff, GW_interp)

plt.xscale("log")
plt.yscale("log")
spec

plt.scatter(agood, bgood)
plt.xscale("log")
plt.yscale("log")
plt.xlim(0.1, 100)
plt.ylim(1, 100)

from random import random

Tgood_nano = []
agood_nano = []
bgood_nano = []

nsamples = 20000
for i in range(nsamples):
    if i % 1000 == 0:
        print(i, len(agood_nano))
    T = 10.0 ** (-2 + 3 * random())
    a = 10.0 ** (-1 + 3 * random())
    b = 10.0 ** (0 + 2 * random())
    GW_t = GW_turb_Andrii(ff * u.Hz, T * u.GeV, a, b, v_w, epsilon_turb)
    GW_s = GW_sound(ff * u.Hz, T * u.GeV, a, b, v_w)
    GW = GW_s + GW_t
    GW_interp = np.interp(fff.value, ff, GW)
    if sum((GW_interp < spec_max) * (GW_interp > spec_min)) == len(fff):
        agood_nano.append(a)
        bgood_nano.append(b)
        Tgood_nano.append(T)
        plt.plot(fff, GW_interp)

print(len(agood_nano))
plt.xscale("log")
plt.yscale("log")

plt.scatter(agood_nano, bgood_nano)
plt.xscale("log")
plt.yscale("log")
plt.xlim(0.1, 100)
plt.ylim(1, 100)

plt.scatter(Tgood_nano, bgood_nano)
plt.xscale("log")
plt.yscale("log")
plt.xlim(0.01, 10)
plt.ylim(1, 100)

from random import random

Tgood_epta = []
agood_epta = []
bgood_epta = []

nsamples = 20000
for i in range(nsamples):
    if i % 1000 == 0:
        print(i, len(agood_epta))
    T = 10.0 ** (-2 + 3 * random())
    a = 10.0 ** (-1 + 3 * random())
    b = 10.0 ** (0 + 2 * random())
    GW_t = GW_turb_Andrii(ff * u.Hz, T * u.GeV, a, b, v_w, epsilon_turb)
    GW_s = GW_sound(ff * u.Hz, T * u.GeV, a, b, v_w)
    GW = GW_s + GW_t
    GW_interp = np.interp(fff.value, ff, GW)
    if sum((GW_interp < spec_max) * (GW_interp > spec_min)) == len(fff):
        agood_epta.append(a)
        bgood_epta.append(b)
        Tgood_epta.append(T)
        plt.plot(fff, GW_interp)

print(len(agood_epta))
plt.xscale("log")
plt.yscale("log")

plt.scatter(agood_epta, bgood_epta)
plt.scatter(agood_nano, bgood_nano)
plt.xscale("log")
plt.yscale("log")
plt.xlim(0.1, 100)
plt.ylim(1, 100)

plt.scatter(Tgood_epta, bgood_epta)
plt.scatter(Tgood_nano, bgood_nano)
plt.xscale("log")
plt.yscale("log")
plt.xlim(0.01, 10)
plt.ylim(1, 100)

d = np.genfromtxt("PPTA.csv")
gammas_ppta = d[:, 0]
As_ppta = 10 ** d[:, 1]
ps_ppta = 5 - gammas_ppta
H0 = 70 * (u.km / u.s) / u.Mpc
plt.plot(gammas_ppta, As_ppta)
plt.plot(gammas_epta, As_epta)

plt.yscale("log")

from random import random

Tgood_ppta = []
agood_ppta = []
bgood_ppta = []

nsamples = 20000
for i in range(nsamples):
    if i % 1000 == 0:
        print(i, len(agood_ppta))
    T = 10.0 ** (-2 + 3 * random())
    a = 10.0 ** (-1 + 3 * random())
    b = 10.0 ** (0 + 2 * random())
    GW_t = GW_turb_Andrii(ff * u.Hz, T * u.GeV, a, b, v_w, epsilon_turb)
    GW_s = GW_sound(ff * u.Hz, T * u.GeV, a, b, v_w)
    GW = GW_s + GW_t
    GW_interp = np.interp(fff.value, ff, GW)
    if sum((GW_interp < spec_max) * (GW_interp > spec_min)) == len(fff):
        agood_ppta.append(a)
        bgood_ppta.append(b)
        Tgood_ppta.append(T)
        plt.plot(fff, GW_interp)

print(len(agood_ppta))
plt.xscale("log")
plt.yscale("log")

plt.scatter(agood_ppta, bgood_ppta)
plt.scatter(agood_epta, bgood_epta)
plt.scatter(agood_nano, bgood_nano)
plt.xscale("log")
plt.yscale("log")
plt.xlim(0.1, 100)
plt.ylim(1, 100)

plt.scatter(Tgood_ppta, bgood_ppta)
plt.scatter(Tgood_epta, bgood_epta)
plt.scatter(Tgood_nano, bgood_nano)
plt.xscale("log")
plt.yscale("log")
plt.xlim(0.01, 10)
plt.ylim(1, 100)

fff = np.logspace(lgfmin, lgfmax, 10) * u.Hz
GW_interp = np.interp(fff.value, ff, GW)
# plt.plot(fff,GW_interp)
fref = 1 / u.yr
print(fref.cgs)
spec_min = np.ones(len(fff))
spec_max = np.zeros(len(fff))
for i in range(len(As)):
    spec = (
        2
        * pi**2
        / 3
        / H0**2
        * fff**2
        * As[i] ** 2
        * (fff / fref) ** (3 - gammas[i])
    ).cgs.value
    spec_min = np.minimum(spec, spec_min)
    spec_max = np.maximum(spec, spec_max)
    # plt.plot(ff,spec)
plt.plot(fff, spec_min, linewidth=4)
plt.plot(fff, spec_max, linewidth=4)

Tgrid = np.logspace(-3, 1, 31)
agrid = np.logspace(-1, 2, 31)
bgrid = np.logspace(0, 2, 31)
Tgood = []
agood = []
bgood = []
for T in Tgrid:
    print(T)
    for a in agrid:
        for b in bgrid:
            GW_t = GW_turb_Andrii(
                ff * u.Hz, T * u.GeV, a, b, v_w, epsilon_turb
            )
            GW_s = GW_sound(ff * u.Hz, T * u.GeV, a, b, v_w)
            GW = GW_s + GW_t
            GW_interp = np.interp(fff.value, ff, GW)
            if sum((GW_interp < spec_max) * (GW_interp > spec_min)) == len(
                fff
            ):
                agood.append(a)
                bgood.append(b)
                Tgood.append(T)
                plt.plot(fff, GW_interp)

plt.xscale("log")
plt.yscale("log")

plt.scatter(agood, bgood)
plt.xscale("log")
plt.yscale("log")
plt.xlim(0.1, 100)
plt.ylim(1, 100)

bin_image = PictureProduct.from_file("Spectrum.png")
from astropy.table import Table

data = [ff, GW_s, GW_t]
names = ("f[Hz]", "Omega_sound_waves", "Omega_turbulence")
spectrum = ODAAstropyTable(Table(data, names=names))

png = bin_image  # http://odahub.io/ontology#ODAPictureProduct
astropy_table = spectrum  # http://odahub.io/ontology#ODAAstropyTable

ffref = fref.cgs.value
print(ffref)
fgrid = np.array([ffref / 1.1, ffref, ffref * 1.1])

print(GW_sound(ffref * np.ones(10), T_star * u.GeV, alpha, beta_H, v_w))

GW_t = GW_turb_Theo(
    ff * u.Hz, T_star * u.GeV, alpha, beta_H, v_w, epsilon_turb
)
plt.plot(ff, GW_t, color="magenta")
GW_t = GW_turb_Andrii(
    ff * u.Hz, T_star * u.GeV, alpha, beta_H, v_w, epsilon_turb
)
plt.plot(ff, h**2 * GW_t, color="blue", linestyle="dashed", label="turbulence")
GW_s = GW_sound(ff * u.Hz, T_star * u.GeV, alpha, beta_H, v_w)
plt.plot(ff, h**2 * GW_s, color="red", linestyle="dotted", label="sound waves")
GW = GW_s + GW_t
plt.plot(ff, h**2 * GW, linewidth=4, color="black", alpha=0.5, label="total")

print(GW_sound(fref, T_star * u.GeV, alpha, beta_H, v_w))

maxGW = max(GW)
ind = np.argmax(GW)
fmax = ff[ind]
plt.xscale("log")
plt.yscale("log")
plt.ylim(3e-15, 5e-9)
plt.xlim(1e-6, 0.1)
plt.grid()
butterfly2 = np.array(
    [
        [8.533691739758777e-8, 0.0000016780550400704105],
        [8.361136018023004e-8, 1.0740085826829729e-7],
        [8.499066585218291e-9, 3.882820258009008e-9],
        [1.029481229115365e-9, 3.6152599901104956e-11],
        [1.0353166841868053e-9, 2.766052794917924e-10],
        [7.3944540651210095e-9, 6.632930188159042e-9],
    ]
)
plt.fill(butterfly2[:, 0], butterfly2[:, 1])
plt.xscale("log")
plt.yscale("log")
plt.legend(loc="upper left")
plt.xlabel("$f$, Hz")
plt.ylabel("$h_0^2\Omega_{gw}(f)$")
plt.savefig("Spectrum.png", format="png", bbox_inches="tight")

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(("out_GW_limits_png", "png_galaxy.output", png))
_oda_outs.append(
    (
        "out_GW_limits_astropy_table",
        "astropy_table_galaxy.output",
        astropy_table,
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
