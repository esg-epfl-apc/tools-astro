#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import matplotlib.pyplot as plt

# we first load necessary packages
import numpy as np
from oda_api.json import CustomJSONEncoder

get_ipython().run_line_magic("matplotlib", "inline")   # noqa: F821
import astropy.constants as const
import astropy.units as u
from astropy.constants import c
from numpy import log, pi, sqrt
from oda_api.data_products import ODAAstropyTable, PictureProduct

# Parameters of the phase transition
# Temperature in GeV
T_star = 0.3  # http://odahub.io/ontology#Energy_GeV

# Numbers of relativistic degrees of freedom
g_star = 50  # http://odahub.io/ontology#Integer
# g_star_s=50 # http://odahub.io/ontology#Integer

# beta/H : rate of the phase transition compared to Hubble rate
beta_H = 1.0

# ratio of the energy density deposited in the bubbles to the radiation energy density
alpha = 1  # http://odahub.io/ontology#Float

# fraction of turbulent energy that goes to gw N.B. arXiv:1004.4187 claims that epsilon_turb=0.05, but checks below show that it is rather 0.01
epsilon_turb = 1.0  # http://odahub.io/ontology#Float

# terminal velocity of bubbles
v_w = 0.99  # http://odahub.io/ontology#Float

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

# Theo, please specify the units of $H_0$ in the cell below

T_star = T_star * u.GeV
T0 = 2.73  # Temperature of the universe now in Kelvin
conversion = 8.617 * 10 ** (-9 - 5)  # conversion factor from kelvin to GeV
H0 = 3.086 * 10 ** (22)  # current Hubble rate
h = 0.7  # dimensionless Hubble Constant
c_s = 3 ** (-0.5)  # speed of sound

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
print("Bubble wall limiting velocity:", v_w)
print(
    "Fraction of released energy in bulk motion:",
    kappa_v,
    ", Thermal energy fraction:",
    kappa_therm,
)
Omega_star = (kappa_v) * alpha / (1 + alpha)
print("Omega in bulk motion:", Omega_star)

# Comoving Hubble rate at the phase transition Eq. 11 of arXiv:1512.06239

def H_star_bis(T_star):
    return (
        16.5e-6
        * (T_star / (100.0 * u.GeV))
        * (g_star / 100) ** (1.0 / 6.0)
        * u.Hz
    )

# With factor of 2pi with definition of scale factor and Hubble rate
def H_star(T_star):
    scale = T0 * conversion * u.GeV / (T_star * (g_star) ** (1 / 3.0))
    Hexp = (
        (4.7 * 10 ** (-42) * g_star) ** (1 / 2)
        * (T_star / (conversion * u.GeV)) ** 2
    ) * u.Hz
    return scale * Hexp / (2 * pi)

R_H_star = (const.c / H_star(T_star)).to(u.Mpc)

print("Reference frequency:", H_star(T_star), ", Horizon size:", R_H_star)

lambda_star = (
    (8 * pi) ** (1 / 3) * max([v_w, c_s]) / (beta_H * H_star(T_star)) * c
)  # size of bubbles at percolation in units of Hubble scale
# there are different distance scales that we can define

R_star = (
    lambda_star  # size of bubbles at percolation (?) in units of Hubble scale
)
Delta_w = (
    sqrt((v_w - c_s) ** 2) / v_w
)  # <1, the width of the bubble wall in units of R_star
R_peak = (
    Delta_w * R_star
)  # the GW spectrum peaks ad the frequency correspoding to the bubble wall width
print(
    "light crossing distance on for the duration of phase transition:",
    lambda_star,
    ", bubble size at percolation:",
    R_star,
)

# Ellis definition of GW_spectrum with broken power law

def GW_old(f, T_star, alpha, beta_H, v_w):
    A = 5.1 * 10 ** (-2)
    a = 2.4
    b = 2.4
    c = 4.0
    delta = 1.0
    scale = T0 * conversion * u.GeV / (T_star * (g_star) ** (1 / 3.0))
    # Hubble rate from Friedman
    Hexp = (
        (4.7 * 10 ** (-42) * g_star) ** (1 / 2)
        * (T_star / (conversion * u.GeV)) ** 2
    ) * u.Hz
    fH = Hexp * scale / (2 * pi)
    Sh = (1 + (f / fH) ** ((-3 + a) / delta)) ** (-delta)
    fp = 0.7 * fH * beta_H
    factor = scale**4 * (Hexp / (70 * 10 ** (3) / (H0) * u.Hz)) ** 2
    return (
        (beta_H) ** 2
        * A
        * (a + b) ** c
        * Sh
        / (b * (f / fp) ** (-a / c) + a * (f / fp) ** (b / c)) ** c
        * factor
    )

# HL model formula

def GW_sound2(f, T_star, alpha, beta_H, v_w):
    # f=f*u.Hz
    (const.c / H_star(T_star)).to(u.Mpc)
    kappa_v = ka_v(alpha, v_w)
    Omega_star = (kappa_v) * alpha / (1 + alpha)
    Omega_gw_tilde = 10 ** (
        -2
    )  # new parameter (that comes from simulations apparently)
    lambda_star = (
        (8 * pi) ** (1 / 3) * max([v_w, c_s]) / (beta_H * H_star(T_star)) * c
    )  # characteristic light-crossing distance scale
    Delta_w = (
        sqrt((v_w - c_s) ** 2) / v_w
    )  # <1, the width of the bubble wall in units of R_star
    s = (
        Delta_w * lambda_star * f / c
    )  # the GW spectrum peaks ad the frequency correspoding to the bubble wall width
    # frequencies are measured in units of the peak frequency
    s1 = lambda_star * f / c
    r_b = Delta_w  # apparently, there are two different notations for the same
    (9 * r_b**4 + 1) / (r_b**4 + 1)
    M = (
        16
        * (1 + r_b ** (-3)) ** (2 / 3.0)
        * (r_b * s1) ** 3
        / ((1 + s1**3) ** (2.0 / 3.0) * (3 + s**2) ** 2)
    )
    mu = 4.78 - 6.27 * r_b + 3.34 * r_b**2
    B_cursive = Omega_gw_tilde / mu  # new parameter
    factor = (Omega_star * lambda_star * H_star(T_star) / c) ** 2 / (
        Omega_star ** (1 / 2.0) + lambda_star * H_star(T_star) / c
    )
    old_factor = (
        3
        * B_cursive
        * factor
        * 1.64
        * 10 ** (-5)
        * (100 / g_star) ** (1 / 3.0)
        * v_w
        / h**2
    )
    return (old_factor * M).cgs

#   SSM model formula

def GW_sound1(f, T_star, alpha, beta_H, v_w):
    # f=f*u.Hz
    (const.c / H_star(T_star)).to(u.Mpc)
    kappa_v = ka_v(alpha, v_w)
    Omega_star = (kappa_v) * alpha / (1 + alpha)
    Omega_gw_tilde = 10 ** (
        -2
    )  # new parameter (that comes from simulations apparently)
    lambda_star = (
        (8 * pi) ** (1 / 3) * max([v_w, c_s]) / (beta_H * H_star(T_star)) * c
    )  # characteristic light-crossing distance scale
    Delta_w = (
        sqrt((v_w - c_s) ** 2) / v_w
    )  # <1, the width of the bubble wall in units of R_star
    s = (
        Delta_w * lambda_star * f / c
    )  # the GW spectrum peaks ad the frequency correspoding to the bubble wall width                    # frequencies are measured in units of the peak frequency
    m = (9 * Delta_w**4 + 1) / (Delta_w**4 + 1)
    M = (
        s**9
        * ((Delta_w**4 + 1) / (Delta_w**4 + ((s) ** 4))) ** 2
        * (5 / (5 - m + m * (s) ** 2)) ** (5.0 / 2.0)
    )
    mu = 4.78 - 6.27 * Delta_w + 3.34 * Delta_w**2
    B_cursive = Omega_gw_tilde / mu  # new parameter
    factor = (Omega_star * lambda_star * H_star(T_star) / c) ** 2 / (
        Omega_star ** (1 / 2.0) + lambda_star * H_star(T_star) / c
    )
    old_factor = (
        3
        * B_cursive
        * factor
        * 1.64
        * 10 ** (-5)
        * (100 / g_star) ** (1 / 3.0)
        * v_w
        / h**2
    )
    return (old_factor * M).cgs

# Eq 1 of the new Overleaf, from Alberto

def GW_turb1(f, T_star, alpha, beta_H, v_w, epsilon_turb):
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

# Range of frequencies used

f = np.logspace(-10, -2, 50) * u.Hz

lambda_star = (
    (8 * pi) ** (1 / 3) * max([v_w, c_s]) / (beta_H * H_star(T_star)) * c
)
kappa_v = ka_v(alpha, v_w)  # kappa_therm=(1-kappa_v)
Omega_star = (kappa_v) * alpha / (1 + alpha)
u_star = sqrt(0.75 * epsilon_turb * Omega_star)  # alfven velocity
dt_fin = 2 * lambda_star / u_star / c
Delta_w = sqrt((v_w - c_s) ** 2) / v_w
delta = 1.0
scale = T0 * conversion * u.GeV / (T_star * (g_star ** (1 / 3.0)))
Hexp = (
    (4.7 * 10 ** (-42) * g_star) ** (1 / 2)
    * (T_star / (conversion * u.GeV)) ** 2
) * u.Hz
fH = Hexp * scale / (2 * pi)
fp = 0.7 * fH * beta_H

# Printing the spectrums

spectrum_turb = GW_turb1(f, T_star, alpha, beta_H, v_w, epsilon_turb)
plt.plot(
    f,
    spectrum_turb,
    label="turbulence",
    color="blue",
    linewidth=4,
    linestyle="dotted",
)
spectrum_sound = GW_sound2(f, T_star, alpha, beta_H, v_w)

spectrum_old = GW_old(f, T_star, alpha, beta_H, v_w)
plt.plot(
    f,
    spectrum_sound,
    label="sound waves",
    color="red",
    linewidth=4,
    linestyle="dashed",
)
plt.plot(
    f, spectrum_sound + spectrum_turb, color="grey", linewidth=4, label="total"
)
# plt.plot(f,spectrum_old,color='purple',linewidth=4,label='total')

# plt.axvline(1/(dt_fin*u.Hz).cgs,linestyle='dashed',color='black',label=r'$dt_{fin}**(-1)$')
# plt.axvline((c/lambda_star/u.Hz).cgs,linestyle='dotted',color='black',label=r'$c/\lambda_*$')
# plt.axvline((c/(Delta_w*lambda_star)/u.Hz).cgs,linestyle='dashdot',color='black',label=r'$c/(\Delta_w \lambda_*)$')
# plt.axvline(c/(R_peak*u.Hz),linestyle='dashed',color='blue',label='peak')

plt.legend()

# Butterfly of compatibility with the NANOgrav data

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
plt.xlabel("$f$, Hz")
plt.ylabel("$h_0^2\Omega_{gw}(f)$")
plt.savefig("Spectrum.png", format="png", bbox_inches="tight")

bin_image = PictureProduct.from_file("Spectrum.png")
from astropy.table import Table

data = [f, spectrum_sound, spectrum_turb]
names = ("f[Hz]", "Omega_sound_waves[MJD]", "Omega_turbulence[MJD]")
spectrum = ODAAstropyTable(Table(data, names=names))

picture = bin_image  # http://odahub.io/ontology#ODAPictureProduct
spectrum_astropy_table = spectrum  # http://odahub.io/ontology#ODAAstropyTable

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(("out_GW_spectrum_picture", "picture_galaxy.output", picture))
_oda_outs.append(
    (
        "out_GW_spectrum_spectrum_astropy_table",
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
