#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
from numpy import exp, log, log10, pi, sqrt
from oda_api.json import CustomJSONEncoder

# ////////////////////////////////////////////////////////////////////////
# // Calculation of synchrotron, inverse Compton (IC) and bremsstrahlung
# // emission from a distribution of electrons with given injection
# // spectrum, with account of synchrotron, IC, bremsstrahlung and Coulomb
# // energy loss.
# //
# // Based on Blumenthal&Gould Rev. Mod. Phys. 42, 237 (1970) [BG]
# ////////////////////////////////////////////////////////////////////////

Gamma_inj = 2.0  # http://odahub.io/ontology#Float ; oda:label "electron injection spectrum powerlaw index"
Ecut = 3e13  # http://odahub.io/ontology#Float ; oda:label "electron injeciton spectrum cut-off"

B = 3e-6  # http://odahub.io/ontology#Float ; oda:label "magnetic field [G]"
n = 1.0e-1  # http://odahub.io/ontology#Float ; oda:label "density of the medium [1/cm3]"
Z = 1.4  # http://odahub.io/ontology#Float ; oda:label "average atomic charge of the medium"
T = 2.73  # http://odahub.io/ontology#Float ; oda:label "Temperature of black body photon background [K]"
Back_norm = 1.0  # http://odahub.io/ontology#Float ; oda:label "Normalization of blackbody photon background (<=1)"

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in ["Gamma_inj", "Ecut", "B", "n", "Z", "T", "Back_norm"]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

# constants
m_e = 5.1e5  # electron mass [eV]
BeV = 7.0e-2  # Gauss-to-eV2 conversion factor
c = 3.0e10  # speed of light
e = 0.085  # electron charge in natural system
hbar = 6.6e-16  # Planck constant in eV-s
eV_cm = hbar * c  # hbar*c
sigma_T = 6.6e-25  # Thomson cross-seciton
alpha = 0.007299  # fine structure constant

Norm = 1.0  # normalization of the electron spectrum
# in principle, dN_e/dE is in [1/cm3/eV]

# Energy binning parameters
Emin = 1.0e-10  # maximum of the energy range
Emax = 10 ** (int(log10(Ecut)) + 2)  # minimum of the energy range
print(int(log10(Ecut)), log10(Emax))
Ndec = 10  # number of energy bins per decade;
NEbins = int((log10(Emax) - log10(Emin)) * Ndec)
print(NEbins)

# Determine the energy binning
energy = np.logspace(log10(Emin), log10(Emax), NEbins)
# bin in which E=m_e and bin in which
# Coulomb los ceases to dominate over bremsstrahlung
for i in range(NEbins):
    if (energy[i] >= m_e) & (energy[i - 1] < m_e):
        j_me = i
    if (energy[i] >= 3.5e8) & (energy[i - 1] < 3.5e8):
        j_brems = i
print(j_me, j_brems)

# electron injection spectrum, set by hands
electrons = (
    Norm
    * (energy / m_e) ** (-Gamma_inj)
    * exp(-energy / Ecut)
    * (energy > m_e)
)

# Synchrotron emissivity
def Fsynch(E, E_e, B):
    if E < E_e:
        # critical energy of synchrotron radiation, BG (4.32)
        E_crit = 3.0 / 2.0 * e * BeV * B / m_e**3 * E_e**2
        # Approximation BG (4.33), (4.34)
        Fx = (
            4 * pi / 2.679 * (E / 2 / E_crit) ** (1.0 / 3.0) * exp(-E / E_crit)
        )
        return e**3 * BeV / m_e / hbar * B * Fx
    else:
        return 0.0

d_energy = energy[1:] - energy[:-1]

def integrate(arr, ind1, ind2):
    arr_mean = sqrt(arr[1:] * arr[:-1])
    return sum(arr_mean[ind1:ind2] * d_energy[ind1:ind2])

# Spectrum of "direct" emission from the freshly injected electrons
# synchrotron emission

synch = np.zeros(NEbins)
# at given emission energy
for i in range(NEbins):
    synch1 = np.zeros(NEbins)
    # from electrons of all energies
    for j in range(j_me, NEbins):
        # Fsynch is in 1/s
        synch1[j] = electrons[j] * Fsynch(energy[i], energy[j], B)
    # synch[i] in 1/cm3/s
    synch[i] = integrate(synch1, 1, -1)

# thermal emission spectrum
def Thermal(ee, T):
    T_eV = 8.6e-5 * T  # temperature in eV
    # thermal spectrum BG (2.58) [1/eV/cm3]
    return ee**2 / (exp(ee / T_eV) - 1) / pi**2 / (hbar * c) ** 3

backgr = energy * Back_norm * Thermal(energy, T)
plt.plot(energy, backgr)
plt.ylim(max(backgr) / 1e3, max(backgr) * 2)
plt.xlim(energy[np.argmax(backgr)] / 100, energy[np.argmax(backgr)] * 100)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("$E$, eV")
plt.ylabel("$E\cdot dn/dE$, 1/cm$^3$")
print(integrate(backgr, 1, -1))

print("    Photon background ")
backgr = energy * Back_norm * Thermal(energy, T)
U_B = (B * BeV) ** 2 / (hbar * c) ** 3 / 8 / pi
U_back = integrate(backgr, 1, -1)
print("    .... Magnetic field energy density:", U_B, "eV/cm3\n")
print("    .... Background energy density:", U_back, "eV/cm3\n")

# inverse Compton emissivity
def Fics(k, l, m):
    # Gamma_\epsilon in BG (2.49)
    b = 4.0 * energy[k] * energy[m] / m_e**2
    emax = b * energy[m] / (1 + b)
    if (k <= l) & (energy[l] < emax):
        # E_1 in BG (2.47)
        z = energy[l] / energy[m]
        # q in BG (2.49)
        q = z / (b * (1.0 - z))
        0.5 * z * z / (1.0 - z)
        # square bracket from (2.48) of BG
        s = (
            2.0 * q * log(q)
            + (1.0 + 2.0 * q) * (1.0 - q)
            + 0.5 * b * b * q * q * (1.0 - q) / (1 + b * q)
        )
        # BG (2.48) in units of cm^2
        return 2.0 * pi * 3 / 4 * sigma_T * m_e**2 * s / energy[k] / energy[m]
    else:
        return 0.0

# Inverse Compton emission
ics = np.zeros(NEbins)
for i in range(NEbins):
    compton = np.zeros(NEbins)
    for j in range(i, NEbins):
        compton1 = np.zeros(NEbins)
        # from background photons of all energies
        for k in range(NEbins):
            # Fics in cm^2; compton1 in 1/s
            compton1[k] = c * backgr[k] / energy[k] * energy[i] * Fics(k, i, j)
        # BG (2.61) [1/cm3/s/eV]
        compton[j] = electrons[j] * integrate(compton1, 1, -1) / energy[j]
    # integration over gamma in BG (2.61) [1/cm3/s]
    ics[i] = integrate(compton, 1, -1)

# bremsstrahlung emissivity
def Fbrems(E, E_e):
    if E < E_e:
        Z = 1.4
        E_fin = E_e - E
        # Blumenthal&Gould (3.31) [cm2/eV]
        return (
            alpha
            * 3.0
            * sigma_T
            / 8.0
            / pi
            / E
            / E_e**2
            * (E_e**2 + E_fin**2 - 2.0 * E_e * E_fin / 3.0)
            * 4
            * (log(2.0 * E_e * E_fin / E / m_e) - 0.5)
            * Z**2
        )
    else:
        return 0.0

# Bremsstrahlung emission
brems = np.zeros(NEbins)
for i in range(NEbins):
    brems1 = np.zeros(NEbins)
    for j in range(i, NEbins):
        # // [1/s/eV2/cm3]
        brems1[j] = c * electrons[j] * n * Fbrems(energy[i], energy[j])
    # integration over electron energies [1/s/cm3]
    brems[i] = energy[i] * integrate(brems1 * (brems1 > 0), 1, -1)

factor = max(max(synch * energy), max(ics * energy), max(brems * energy))
plt.plot(energy, synch * energy / factor, color="black", label="synchrotron")
plt.plot(
    energy,
    ics * energy / factor,
    color="black",
    linestyle="dashed",
    label=r"inverse Compton",
)
plt.plot(
    energy,
    brems * energy / factor,
    color="black",
    linestyle="dotted",
    label=r"Bremsstrahlung",
)
plt.xscale("log")
plt.yscale("log")
plt.ylim(1e-10, 10)
m = (
    synch * energy / factor + ics * energy / factor + brems * energy / factor
) > 1e-5
# plt.xlim(energy[m][0],energy[m][-1])
plt.legend(loc="lower right")
plt.xlabel("Energy, eV")
plt.ylabel("Flux, erg/(cm$^2$s)")
plt.savefig("Spectrum.png", format="png", bbox_inches="tight")

from oda_api.data_products import ODAAstropyTable, PictureProduct

bin_image = PictureProduct.from_file("Spectrum.png")

from astropy.table import Table

data = [
    energy,
    synch * energy / factor,
    ics * energy / factor,
    brems * energy / factor,
]
names = ("E[eV]", "E2dNdE_synch", "E2dNdE_ics", "E2dNdE_brems")
spec = ODAAstropyTable(Table(data, names=names))

spectrum_png = bin_image  # http://odahub.io/ontology#ODAPictureProduct
spectrum_table = spec  # http://odahub.io/ontology#ODAAstropyTable

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_Synchrotron_Inverse_Compton_Bremsstrahlung_spectrum_png",
        "spectrum_png_galaxy.output",
        spectrum_png,
    )
)
_oda_outs.append(
    (
        "out_Synchrotron_Inverse_Compton_Bremsstrahlung_spectrum_table",
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
