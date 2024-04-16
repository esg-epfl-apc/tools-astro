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
import astropy.units as u
from astropy.constants import c
from numpy import exp, log, pi, sqrt
from oda_api.api import ProgressReporter
from oda_api.data_products import PictureProduct

# Parameters of the phase transition

# fraction of turbulent energy that goes to gw N.B. arXiv:1004.4187 claims that epsilon_turb=0.05, but checks below show that it is rather 0.01
epsilon_turb = 1  # http://odahub.io/ontology#Float

# Numbers of relativistic degrees of freedom
g_star = 20  # http://odahub.io/ontology#Integer

Tmax_yrs = 10.0  # http://odahub.io/ontology#Float
Tmin_yrs = 2.0  # http://odahub.io/ontology#Float
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
pr = ProgressReporter()
pr.report_progress(stage="Progress", progress=1.0)

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

T_star = 0.2

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

    Omega_B = Omega_star / 2.0
    Omega_gamma = 2 / g_star
    B = 3e-6 * (Omega_B / Omega_gamma) ** 0.5
    return res, lambda_star, B

H0 = 70 * (u.km / u.s) / u.Mpc

d = np.genfromtxt("NANOGrav23.csv")
gammas_nano = d[:, 0]
As_nano = 10 ** d[:, 1]
ps_nano = 5 - gammas_nano

d = np.genfromtxt("EPTA.csv")
gammas_epta = d[:, 0]
As_epta = 10 ** d[:, 1]
ps_epta = 5 - gammas_epta

d = np.genfromtxt("PPTA.csv")
gammas_ppta = d[:, 0]
As_ppta = 10 ** d[:, 1]
ps_ppta = 5 - gammas_ppta

fref = (1 / u.yr).cgs.value
lgfmin = np.log10(fref / Tmax_yrs)
lgfmax = np.log10(fref / Tmin_yrs)
fff = np.logspace(lgfmin, lgfmax, 10) * u.Hz
min_nano = np.ones(len(fff))
max_nano = np.zeros(len(fff))
for i in range(len(As_nano)):
    spec = (
        2
        * pi**2
        / 3
        / H0**2
        * fff**2
        * As_nano[i] ** 2
        * (fff / fref) ** (3 - gammas_nano[i])
    ).cgs.value
    min_nano = np.minimum(spec, min_nano)
    max_nano = np.maximum(spec, max_nano)
    # plt.plot(ff,spec)
plt.fill_between(
    fff.value, min_nano, max_nano, color="red", alpha=0.5, label="NANOGrav"
)
min_epta = np.ones(len(fff))
max_epta = np.zeros(len(fff))
for i in range(len(As_epta)):
    spec = (
        2
        * pi**2
        / 3
        / H0**2
        * fff**2
        * As_epta[i] ** 2
        * (fff / fref) ** (3 - gammas_epta[i])
    ).cgs.value
    min_epta = np.minimum(spec, min_epta)
    max_epta = np.maximum(spec, max_epta)
    # plt.plot(ff,spec)
plt.fill_between(
    fff.value, min_epta, max_epta, color="blue", alpha=0.5, label="EPTA"
)
min_ppta = np.ones(len(fff))
max_ppta = np.zeros(len(fff))
for i in range(len(As_ppta)):
    spec = (
        2
        * pi**2
        / 3
        / H0**2
        * fff**2
        * As_ppta[i] ** 2
        * (fff / fref) ** (3 - gammas_ppta[i])
    ).cgs.value
    min_ppta = np.minimum(spec, min_ppta)
    max_ppta = np.maximum(spec, max_ppta)
    # plt.plot(ff,spec)
plt.fill_between(
    fff.value, min_ppta, max_ppta, color="green", alpha=0.5, label="PPTA"
)

# PBH constraints
def PBH(alpha):
    return (5.2 * (1.0 - exp(-1.1 * (abs(alpha - 1.0)) ** (1 / 3))) + 1.0) * (
        alpha > 1
    )

Tgrid = np.logspace(-2.6, 1, 81)
agrid = np.logspace(-1, 2, 61)
bgrid = np.logspace(0, 2, 41)
amin_b_nano_eps1 = 100.0 * np.ones(len(bgrid))
amax_b_nano_eps1 = 0.0 * np.ones(len(bgrid))
amin_b_epta_eps1 = 100.0 * np.ones(len(bgrid))
amax_b_epta_eps1 = 0.0 * np.ones(len(bgrid))
amin_b_ppta_eps1 = 100.0 * np.ones(len(bgrid))
amax_b_ppta_eps1 = 0.0 * np.ones(len(bgrid))
amin_T_nano_eps1 = 100.0 * np.ones(len(Tgrid))
amax_T_nano_eps1 = 0.0 * np.ones(len(Tgrid))
amin_T_epta_eps1 = 100.0 * np.ones(len(Tgrid))
amax_T_epta_eps1 = 0.0 * np.ones(len(Tgrid))
amin_T_ppta_eps1 = 100.0 * np.ones(len(Tgrid))
amax_T_ppta_eps1 = 0.0 * np.ones(len(Tgrid))

Tmin_b_nano_eps1 = 100.0 * np.ones(len(bgrid))
Tmax_b_nano_eps1 = 0.0 * np.ones(len(bgrid))
Tmin_b_epta_eps1 = 100.0 * np.ones(len(bgrid))
Tmax_b_epta_eps1 = 0.0 * np.ones(len(bgrid))
Tmin_b_ppta_eps1 = 100.0 * np.ones(len(bgrid))
Tmax_b_ppta_eps1 = 0.0 * np.ones(len(bgrid))
B_nano_eps1 = []
lam_nano_eps1 = []
B_ppta_eps1 = []
lam_ppta_eps1 = []
B_epta_eps1 = []
lam_epta_eps1 = []
for i, T in enumerate(Tgrid):
    print(T, len(B_nano_eps1), len(B_epta_eps1), len(B_ppta_eps1))
    pr.report_progress(stage="Progress", progress=1 + 94.0 * (i / len(Tgrid)))
    for j, a in enumerate(agrid):
        for k, b in enumerate(bgrid):
            if b > PBH(a):
                GW_t = GW_turb_Andrii(
                    ff * u.Hz, T * u.GeV, a, b, v_w, epsilon_turb
                )
                lam = GW_t[1].cgs.value
                B = GW_t[2]
                GW_t = GW_t[0]
                GW_s = GW_sound(ff * u.Hz, T * u.GeV, a, b, v_w)
                GW = GW_s + GW_t
                GW_interp = np.interp(fff.value, ff, GW)
                if sum((GW_interp < max_nano) * (GW_interp > min_nano)) == len(
                    fff
                ):
                    if a < amin_b_nano_eps1[k]:
                        amin_b_nano_eps1[k] = a
                    if a > amax_b_nano_eps1[k]:
                        amax_b_nano_eps1[k] = a
                    if a < amin_T_nano_eps1[i]:
                        amin_T_nano_eps1[i] = a
                    if a > amax_T_nano_eps1[i]:
                        amax_T_nano_eps1[i] = a
                    if T < Tmin_b_nano_eps1[k]:
                        Tmin_b_nano_eps1[k] = T
                    if T > Tmax_b_nano_eps1[k]:
                        Tmax_b_nano_eps1[k] = T
                    B_nano_eps1.append(B)
                    lam_nano_eps1.append(lam)
                if sum((GW_interp < max_epta) * (GW_interp > min_epta)) == len(
                    fff
                ):
                    if a < amin_b_epta_eps1[k]:
                        amin_b_epta_eps1[k] = a
                    if a > amax_b_epta_eps1[k]:
                        amax_b_epta_eps1[k] = a
                    if a < amin_T_epta_eps1[i]:
                        amin_T_epta_eps1[i] = a
                    if a > amax_T_epta_eps1[i]:
                        amax_T_epta_eps1[i] = a
                    if T < Tmin_b_epta_eps1[k]:
                        Tmin_b_epta_eps1[k] = T
                    if T > Tmax_b_epta_eps1[k]:
                        Tmax_b_epta_eps1[k] = T
                    B_epta_eps1.append(B)
                    lam_epta_eps1.append(lam)
                if sum((GW_interp < max_ppta) * (GW_interp > min_ppta)) == len(
                    fff
                ):
                    if a < amin_b_ppta_eps1[k]:
                        amin_b_ppta_eps1[k] = a
                    if a > amax_b_ppta_eps1[k]:
                        amax_b_ppta_eps1[k] = a
                    if a < amin_T_ppta_eps1[i]:
                        amin_T_ppta_eps1[i] = a
                    if a > amax_T_ppta_eps1[i]:
                        amax_T_ppta_eps1[i] = a
                    if T < Tmin_b_ppta_eps1[k]:
                        Tmin_b_ppta_eps1[k] = T
                    if T > Tmax_b_ppta_eps1[k]:
                        Tmax_b_ppta_eps1[k] = T
                    B_ppta_eps1.append(B)
                    lam_ppta_eps1.append(lam)

lam_bins = np.logspace(-7, -5, 41)
B_bins = np.logspace(-7, -5, 41)

h = plt.hist2d(
    np.array(lam_nano_eps1) / 3e24,
    B_nano_eps1,
    bins=[lam_bins, B_bins],
    vmax=1,
)
plt.colorbar()
cnt = plt.contour(
    sqrt(lam_bins[1:] * lam_bins[:-1]),
    sqrt(B_bins[1:] * B_bins[:-1]),
    np.transpose(h[0]),
    levels=[1],
    colors="white",
)
cont = cnt.get_paths()[0].vertices
B_nanos_eps1 = cont[:, 1]
lam_nanos_eps1 = cont[:, 0]

plt.xscale("log")
plt.yscale("log")

lam_bins = np.logspace(-7, -5, 41)
B_bins = np.logspace(-7, -5, 41)

h = plt.hist2d(
    np.array(lam_epta_eps1) / 3e24,
    B_epta_eps1,
    bins=[lam_bins, B_bins],
    vmax=1,
)
plt.colorbar()
cnt = plt.contour(
    sqrt(lam_bins[1:] * lam_bins[:-1]),
    sqrt(B_bins[1:] * B_bins[:-1]),
    np.transpose(h[0]),
    levels=[1],
    colors="white",
)
cont = cnt.get_paths()[0].vertices
B_eptas_eps1 = cont[:, 1]
lam_eptas_eps1 = cont[:, 0]

plt.xscale("log")
plt.yscale("log")

lam_bins = np.logspace(-7, -5, 41)
B_bins = np.logspace(-7, -5, 41)

h = plt.hist2d(
    np.array(lam_ppta_eps1) / 3e24,
    B_ppta_eps1,
    bins=[lam_bins, B_bins],
    vmax=1,
)
plt.colorbar()
cnt = plt.contour(
    sqrt(lam_bins[1:] * lam_bins[:-1]),
    sqrt(B_bins[1:] * B_bins[:-1]),
    np.transpose(h[0]),
    levels=[1],
    colors="white",
)
cont = cnt.get_paths()[0].vertices
B_pptas_eps1 = cont[:, 1]
lam_pptas_eps1 = cont[:, 0]

plt.xscale("log")
plt.yscale("log")

import numpy as np
from matplotlib import pyplot as plt

fig, ax = plt.subplots(figsize=(7, 7))
x = [1e-8, 1e3]
y1 = [3e-6, 3e-6]
y2 = [3e-5, 3e-5]
plt.fill_between(x, y1, y2, color="grey", alpha=0.5)
ax.fill(lam_nanos_eps1, B_nanos_eps1, color="red", alpha=0.3, label="NANOGrav")
ax.fill(lam_eptas_eps1, B_eptas_eps1, color="blue", alpha=0.3, label="EPTA")
ax.fill(lam_pptas_eps1, B_pptas_eps1, color="green", alpha=0.3, label="PPTA")
ax.set_xscale("log")
ax.set_yscale("log")

x = np.logspace(-8, 3, 10)
y = 10 ** (-8.5) * x
ax.plot(x, y, color="grey", linewidth=4, linestyle="dashed")
y = 10 ** (-6.8) * x
ax.plot(x, y, color="grey", linewidth=4)
x = np.logspace(-6.5, -3.1, 10)
y = 2e-6 * (x / x[0]) ** (-5 / 4.0)
ax.annotate(
    "",
    xy=(x[-1], y[-1]),
    xytext=(x[0], y[0]),
    arrowprops={"arrowstyle": "->", "color": "black"},
)
x = np.logspace(-5.7, -2.4, 10)
y = 7e-6 * (x / x[0]) ** (-5 / 4.0)
ax.annotate(
    "",
    xy=(x[-1], y[-1]),
    xytext=(x[0], y[0]),
    arrowprops={"arrowstyle": "->", "color": "black"},
)

x = np.logspace(-6.5, -3.8, 10)
y = 2e-6 * (x / x[0]) ** (-5 / 2.0)
ax.annotate(
    "",
    xy=(x[-1], y[-1]),
    xytext=(x[0], y[0]),
    arrowprops={"arrowstyle": "->", "color": "black", "linestyle": "dashed"},
)
x = np.logspace(-5.7, -3.1, 10)
y = 7e-6 * (x / x[0]) ** (-5 / 2.0)
ax.annotate(
    "",
    xy=(x[-1], y[-1]),
    xytext=(x[0], y[0]),
    arrowprops={"arrowstyle": "->", "color": "black", "linestyle": "dashed"},
)

d = np.genfromtxt("RoperPol22.csv")
ax.plot(d[:, 0], d[:, 1], color="black", linewidth=4, label="2201.05630")

d = np.genfromtxt("MAGIC.csv")
ax.fill_between(d[:, 0], np.zeros(len(d)), d[:, 1], color="grey", alpha=0.5)
ax.set_xlim(1e-8, 1e3)
ax.set_ylim(5e-18, 2e-5)
ax.set_xlabel(r"$\lambda_B$, Mpc")
ax.set_ylabel(r"$B$, G")
plt.text(1, 3e-17, "MAGIC '22", color="black")
plt.text(
    1e-6,
    0.1e-14,
    "Alfvenic larges processed eddies",
    color="black",
    rotation=42,
)

plt.text(
    0.15e-7,
    0.01e-13,
    "Helicity fluctuaitons / reconnection controlled",
    color="black",
    rotation=41.5,
)
ax.legend(loc="lower left")

y = [1.5e-12 * 1e1, 1.5e-12, 1.5e-12]
x = [0.3e-3, 0.03, 1e3]
plt.plot(x, y, color="olive", linewidth=4, linestyle="dotted")
ax.annotate(
    "",
    xytext=(1, 1.0e-13),
    xy=(1, 1.5e-12),
    arrowprops={"arrowstyle": "->", "color": "olive", "linewidth": 4},
)
plt.text(0.5, 2e-12, "CTA", color="olive")

x = [0.7e-3, 1e3]
y = [3e-11, 3e-11]
plt.plot(x, y, color="olive", linewidth=4, linestyle="dotted")
ax.annotate(
    "",
    xytext=(1, 3e-10),
    xy=(1, 3e-11),
    arrowprops={"arrowstyle": "->", "color": "olive", "linewidth": 4},
)
plt.text(0.5, 1.1e-11, "CMB", color="olive")

# inset axes....
x1, x2, y1, y2 = 2e-7, 3e-6, 7e-7, 1e-5  # subregion of the original image
axins = ax.inset_axes(
    [0.5, 0.7, 0.47, 0.4], xlim=(x1, x2), ylim=(y1, y2)
)  # xticklabels=[], yticklabels=[])
# axins.imshow(Z2, extent=extent, origin="lower")
axins.fill(lam_nanos_eps1, B_nanos_eps1, color="red", alpha=0.3, linewidth=0)
axins.fill(lam_eptas_eps1, B_eptas_eps1, color="blue", alpha=0.3, linewidth=0)
axins.fill(lam_pptas_eps1, B_pptas_eps1, color="green", alpha=0.3, linewidth=0)
axins.set_xscale("log")
axins.set_yscale("log")

ax.indicate_inset_zoom(axins, edgecolor="black")
plt.savefig("B_lambdaB.png", bbox_inches="tight")

fig = plt.figure()
mask = Tmax_b_nano_eps1 > 0
plt.fill_between(
    bgrid[mask],
    Tmin_b_nano_eps1[mask],
    Tmax_b_nano_eps1[mask],
    alpha=0.3,
    color="red",
    label="NANOGrav",
)
mask = Tmax_b_epta_eps1 > 0
plt.fill_between(
    bgrid[mask],
    Tmin_b_epta_eps1[mask],
    Tmax_b_epta_eps1[mask],
    alpha=0.3,
    color="blue",
    label="EPTA",
)
mask = Tmax_b_ppta_eps1 > 0
plt.fill_between(
    bgrid[mask],
    Tmin_b_ppta_eps1[mask],
    Tmax_b_ppta_eps1[mask],
    alpha=0.3,
    color="green",
    label="PPTA",
)
d = np.genfromtxt("Ellis.csv")
lgT = d[:, 0]
lgb = d[:, 1]
plt.plot(10**lgb, 10**lgT, color="black", label="2308.08546")
d = np.genfromtxt("NANO_bubble.csv")
lgT = d[:, 0]
lgb = -d[:, 1]
plt.plot(
    10**lgb, 10**lgT, color="black", linestyle="dashed", label="2306.162196"
)
d = np.genfromtxt("NANO_sound.csv")
lgT = d[:, 0]
lgb = -d[:, 1]
plt.plot(10**lgb, 10**lgT, color="black", linestyle="dashed")
plt.axhline(0.160, color="black", linestyle="dotted")

plt.xscale("log")
plt.yscale("log")
plt.ylabel(r"$T$ [GeV]")
plt.xlabel(r"$\beta/H$")
plt.xlim(0.8, 80)
plt.ylim(0.001, 10)
plt.legend(loc="lower left")
plt.savefig("T_Beta.png", bbox_inches="tight")

fig = plt.figure()
mask = amax_b_nano_eps1 > 0
plt.fill_between(
    bgrid[mask],
    amin_b_nano_eps1[mask],
    amax_b_nano_eps1[mask],
    alpha=0.3,
    color="red",
    label="NANOGrav",
)
mask = amax_b_epta_eps1 > 0
plt.fill_between(
    bgrid[mask],
    amin_b_epta_eps1[mask],
    amax_b_epta_eps1[mask],
    alpha=0.3,
    color="blue",
    label="EPTA",
)
mask = amax_b_ppta_eps1 > 0
plt.fill_between(
    bgrid[mask],
    amin_b_ppta_eps1[mask],
    amax_b_ppta_eps1[mask],
    alpha=0.3,
    color="green",
    label="PPTA",
)
bpbh = PBH(agrid)
x = [2, 3]
y = [1, 1]
y1 = [2, 2]
plt.fill_between(x, y, y1, color="white", linewidth=0)
plt.fill(bpbh, agrid, color="white")
plt.xscale("log")
plt.yscale("log")
plt.xlim(0.8, 80)
plt.ylim(0.2, 80)
plt.xlabel(r"$\beta/H$")
plt.ylabel(r"$\alpha$")
plt.legend(loc="upper left")
plt.savefig("Alpha_Beta.png", bbox_inches="tight")

fig = plt.figure()
mask = amax_T_nano_eps1 > 0
plt.fill_between(
    Tgrid[mask],
    amin_T_nano_eps1[mask],
    amax_T_nano_eps1[mask],
    alpha=0.3,
    color="red",
    label="NANOGrav",
)
mask = amax_T_epta_eps1 > 0
plt.fill_between(
    Tgrid[mask],
    amin_T_epta_eps1[mask],
    amax_T_epta_eps1[mask],
    alpha=0.3,
    color="blue",
    label="EPTA",
    linewidth=0,
)
mask = amax_T_ppta_eps1 > 0
plt.fill_between(
    Tgrid[mask],
    amin_T_ppta_eps1[mask],
    amax_T_ppta_eps1[mask],
    alpha=0.3,
    color="green",
    label="PPTA",
    linewidth=0,
)
d = np.genfromtxt("NANO_alpha_T.csv")
lgT = d[:, 0]
lga = d[:, 1]
plt.plot(
    10**lgT, 10**lga, color="black", label="2306.16219", linestyle="dashed"
)
d = np.genfromtxt("NANO_alpha_T1.csv")
lgT = d[:, 0]
lga = d[:, 1]
plt.plot(10**lgT, 10**lga, color="black", linestyle="dashed")
plt.axvline(0.16, color="black", linestyle="dotted")
plt.xscale("log")
plt.yscale("log")
plt.xlim(0.001, 10)
plt.ylim(0.2, 80)
plt.xlabel(r"$T$ [GeV]")
plt.ylabel(r"$\alpha$")
plt.legend(loc="upper left")
plt.savefig("Alpha_T.png", bbox_inches="tight")

bin_image1 = PictureProduct.from_file("B_lambdaB.png")
bin_image2 = PictureProduct.from_file("T_Beta.png")
bin_image3 = PictureProduct.from_file("Alpha_Beta.png")
bin_image4 = PictureProduct.from_file("Alpha_T.png")

B_lambdaB_png = bin_image1  # http://odahub.io/ontology#ODAPictureProduct
T_Beta_png = bin_image2  # http://odahub.io/ontology#ODAPictureProduct
Alpha_Beta_png = bin_image3  # http://odahub.io/ontology#ODAPictureProduct
Alpha_T_png = bin_image4  # http://odahub.io/ontology#ODAPictureProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_Phase_transition_parameters_B_lambdaB_png",
        "B_lambdaB_png_galaxy.output",
        B_lambdaB_png,
    )
)
_oda_outs.append(
    (
        "out_Phase_transition_parameters_T_Beta_png",
        "T_Beta_png_galaxy.output",
        T_Beta_png,
    )
)
_oda_outs.append(
    (
        "out_Phase_transition_parameters_Alpha_Beta_png",
        "Alpha_Beta_png_galaxy.output",
        Alpha_Beta_png,
    )
)
_oda_outs.append(
    (
        "out_Phase_transition_parameters_Alpha_T_png",
        "Alpha_T_png_galaxy.output",
        Alpha_T_png,
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
