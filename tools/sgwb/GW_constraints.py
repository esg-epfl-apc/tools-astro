#!/usr/bin/env python
# coding: utf-8

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

# flake8: noqa

import json
import os

import matplotlib.pyplot as plt

# Loading necessary packages
import numpy as np

get_ipython().run_line_magic("matplotlib", "inline")   # noqa: F821
from math import *

import astropy.constants as const
import astropy.units as u
from astropy.constants import c
from numpy import log, pi, sqrt
from scipy.integrate import quad

# Parameters of the phase transition
T_star = 3.0 * u.GeV
g_star = 20
g_star_s = 20

T0 = 2.73  # Temperature of the universe now
conversion = 8.617 * 10 ** (-9 - 5)  # conversion factor from kelvin to GeV
H0 = 70 * 10**3 / (3.086 * 10 ** (22))  # current Hubble rate

beta_H = (
    1000.0  # beta/H : rate of the phase transition compared to Hubble rate
)
alpha = 2  # ratio of the energy density deposited in the bubbles to the radiation energy density
epsilon_turb = 0.1  # fraction of turbulent energy that goes to gw N.B. arXiv:1004.4187 claims that epsilon_turb=0.05, but checks below show that it is rather 0.01

# --------------------
v_w = 0.99  # terminal velocity of bubbles
h = 0.7  # dimensionless Hubble Constant
c_s = 3 ** (-0.5)  # speed of sound

taille = [1, 30, 30, 30]
taille_old = taille

vlist = np.linspace(0.99, 1, taille[0])
alphalist = np.logspace(-1, 2, taille[2])
betalist = np.logspace(0, 2, taille[3])
Templist = np.logspace(-3, 1, taille[1])

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

def H_star(T_star):
    return (
        16.5e-6
        * (T_star / (100.0 * u.GeV))
        * (g_star / 100) ** (1.0 / 6.0)
        * u.Hz
    )

# With factor of 2pi with definition of scale factor and Hubble rate
def H_star_bis(T_star):
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
    factor = scale**4 * (Hexp / (H0 * u.Hz)) ** 2
    return (
        1
        / (beta_H) ** 2
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

# Just a small modification on units of G_old to be able to use the quad module to do the integral of the spectrum

def GW_tot(f, T_star, alpha, beta_H, v_w, epsilon_turb):
    return GW_old(f * u.Hz, T_star, alpha, beta_H, v_w)
    # return (GW_turb1(f*u.Hz,T_star,alpha,beta_H,v_w,epsilon_turb)+GW_sound2(f*u.Hz,T_star,alpha,beta_H,v_w))

print(GW_tot(10 ** (-8), 0.3 * u.GeV, 6, 6, 0.9, 1))
quad(GW_tot, 0, 10 ** (-5), args=(0.34 * u.GeV, 6, 10, 0.9, 0.1))

# PBH constraints
def PBH_bis(alpha):
    return 3.8

def PBH(alpha):
    if alpha <= 1.0:
        return 0
    # if alpha<1.5 and alpha>1.0:
    #   return 2.8/0.5*alpha-2.3/0.5
    if alpha > 1.0:
        return 5.2 * (1.0 - exp(-1.1 * (alpha - 1.0) ** (1 / 3))) + 1.0

xx = np.arange(0, 100, 0.01)
yy = []
for k in range(len(xx)):
    yy.append(PBH(xx[k]))
plt.plot(xx, yy)

def Bturb(T_star, alpha, beta_H, v_w, epsilon_turb):
    kappa_v = ka_v(alpha, v_w)
    (1 - kappa_v)
    # Here are the prescriptions for the magnetic field from 1907.04315
    R_H_star = (const.c / H_star(T_star)).to(u.Mpc)
    Omega_B_star = epsilon_turb * (kappa_v) * alpha / (1 + alpha)
    lambda_B_star = (
        (8 * pi) ** (1.0 / 3.0) * max(v_w, c_s) / beta_H * R_H_star
    )  # size of bubbles at percolation (?) in units of Hubble scale
    # 3e-6 G would be in equipartition with CMB. We estimate comoving magnetic field at production based on this:
    Omega_ph = 2 / g_star
    B_star = 3e-6 * u.G * np.sqrt(Omega_B_star / Omega_ph)
    # print('Magnetic field at production: lambda=',lambda_B_star,', B=', B_star)
    # Helical magnetic field evolves with gamma=0.5:
    gamma = 0.5
    lambda_B_h = (
        u.Mpc
        * (B_star / (10 ** (-8.5) * u.G)) ** (1 / (1 + gamma))
        * (lambda_B_star / u.Mpc) ** (gamma / (1 + gamma))
    )
    B_h = 10 ** (-8.5) * u.G * (lambda_B_h / u.Mpc)
    # print('Helical magnetic field today: lambda=',lambda_B_h,' B=',B_h)
    # Non-Helical magnetic field evolves with gamma=1:
    gamma = 1.25
    lambda_B_nh = (
        u.Mpc
        * (B_star / (10.0 ** (-8.5) * u.G)) ** (1 / (1 + gamma))
        * (lambda_B_star / u.Mpc) ** (gamma / (1 + gamma))
    )
    B_nh = 10 ** (-8.5) * u.G * (lambda_B_nh / u.Mpc)
    # print('Non-helical magnetic field today: lambda=',lambda_B_nh,', B=',B_nh)
    # print(Omega_B_star,lambda_B_star)
    return (B_star, lambda_B_star, B_h, lambda_B_h, B_nh, lambda_B_nh)

# PTA data loading (approximation of constraints as ellipses)

DATA = np.loadtxt("Documents/PPTA_15.txt")
print(DATA)

u_ell3_PPTA = DATA[0][0]
v_ell3_PPTA = DATA[1][1]
a_ell3_PPTA = (
    (DATA[4][0] - u_ell3_PPTA) ** 2 + (DATA[4][1] - v_ell3_PPTA) ** 2
) ** (1 / 2)
print(a_ell3_PPTA)
b_ell3_PPTA = (
    (DATA[10][0] - u_ell3_PPTA) ** 2 + (DATA[10][1] - v_ell3_PPTA) ** 2
) ** (1 / 2)
print(b_ell3_PPTA)
u_ell2_PPTA = DATA[0][0]
v_ell2_PPTA = DATA[1][1]
a_ell2_PPTA = (
    (DATA[3][0] - u_ell3_PPTA) ** 2 + (DATA[3][1] - v_ell3_PPTA) ** 2
) ** (1 / 2)
print(a_ell2_PPTA)
b_ell2_PPTA = (
    (DATA[9][0] - u_ell3_PPTA) ** 2 + (DATA[9][1] - v_ell3_PPTA) ** 2
) ** (1 / 2)
print(b_ell2_PPTA)
u_ell_PPTA = DATA[0][0]
v_ell_PPTA = DATA[1][1]
a_ell_PPTA = (
    (DATA[2][0] - u_ell3_PPTA) ** 2 + (DATA[2][1] - v_ell3_PPTA) ** 2
) ** (1 / 2)
print(a_ell_PPTA)
b_ell_PPTA = (
    (DATA[8][0] - u_ell3_PPTA) ** 2 + (DATA[8][1] - v_ell3_PPTA) ** 2
) ** (1 / 2)
print(b_ell_PPTA)
theta_PPTA = -0.32

# output gathering
_galaxy_meta_data = {}

with open(os.path.join(_galaxy_wd, "galaxy.json"), "w") as fd:
    json.dump(_galaxy_meta_data, fd)
print("*** Job finished successfully ***")
