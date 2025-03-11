#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

# oda:usesResource oda:CRBeamS3 .
# oda:CRBeamS3 a oda:S3 .
# oda:CRBeamS3 oda:resourceBindingEnvVarName "S3_CREDENTIALS" .

src_name = "NGC 1365"  # http://odahub.io/ontology#AstrophysicalObject

z_start = 0  # http://odahub.io/ontology#Float
Npart = 2000  # http://odahub.io/ontology#Integer ; oda:lower_limit 1 ; oda:upper_limit 100000
particle_type = "gamma"  # http://odahub.io/ontology#String ; oda:allowed_value "gamma","electron","proton"
Emax = 30  # http://odahub.io/ontology#Energy_TeV
Emin = 0.01  # http://odahub.io/ontology#Energy_TeV
EminSource = 1.0  # http://odahub.io/ontology#Energy_TeV
Gamma = 2.0  # http://odahub.io/ontology#Float
EGMF_fG = 10  # http://odahub.io/ontology#Float
lmaxEGMF_Mpc = 5  # http://odahub.io/ontology#Float
jet_half_size = 5.0  # http://odahub.io/ontology#degree
jet_direction = 0.0  # http://odahub.io/ontology#degree
psf = 1.0  # http://odahub.io/ontology#degree
window_size_RA = 2.0  # http://odahub.io/ontology#degree
window_size_DEC = 1.0  # http://odahub.io/ontology#degree
EBL = "Franceschini 2017"  # http://odahub.io/ontology#String ; oda:allowed_value "Franceschini 2017","Stecker 2016 lower limit","Stecker 2016 upper limit","Inoue 2012 Baseline","Inoue 2012 lower limit","Inoue 2012 upper limit"

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in [
    "src_name",
    "z_start",
    "Npart",
    "particle_type",
    "Emax",
    "Emin",
    "EminSource",
    "Gamma",
    "EGMF_fG",
    "lmaxEGMF_Mpc",
    "jet_half_size",
    "jet_direction",
    "psf",
    "window_size_RA",
    "window_size_DEC",
    "EBL",
]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

get_ipython().run_cell_magic(   # noqa: F821
    "bash",
    "",
    "if [ ! -f utils.py ]\nthen\n    git clone https://gitlab.renkulab.io/astronomy/mmoda/crbeam.git tmp_src\n    cp tmp_src/*.sh tmp_src/*.py ./\nfi\n",
)

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.nddata import StdDevUncertainty
from astropy.table import Table
from gammapy.maps import Map, MapAxis
from oda_api.api import ProgressReporter
from oda_api.data_products import (
    LightCurveDataProduct,
    NumpyDataProduct,
    ODAAstropyTable,
    PictureProduct,
)
from specutils import Spectrum1D
from utils import find_redshift

if z_start <= 0:
    z_start = find_redshift(src_name)

source = SkyCoord.from_name(src_name, frame="icrs", parse=False, cache=True)

pr = ProgressReporter()
pr.report_progress(
    stage="Initialization",
    progress=0,
    substage="spectra",
    subprogress=30,
    message="some message",
)

import subprocess

import light_curve as lc
import numpy as np
from crbeam import CRbeam
from matplotlib import pyplot as plt
from spec_plot import plot_spec

EGMF = EGMF_fG * 1e-15

# internal units are eV
Emax *= 1e12
Emin *= 1e12
EminSource *= 1e12

background = {
    "Franceschini 2017": 12,
    "Stecker 2016 lower limit": 10,
    "Stecker 2016 upper limit": 11,
    "Inoue 2012 Baseline": 3,
    "Inoue 2012 lower limit": 4,
    "Inoue 2012 upper limit": 5,
}[EBL]

prog = CRbeam(
    z=z_start,
    nparticles=Npart,
    primary=particle_type,
    emax=Emax,
    emin=Emin,
    emin_source=EminSource,
    EGMF=EGMF,
    lmaxEGMF=lmaxEGMF_Mpc,
    background=background,
)
cmd = prog.command
cmd

# Here is one way to call CRbeam without global cache
# prog.run(force_overwrite=False)
# Here we call CRbeam with global cache
# prog.run_cached(overwrite_local_cache=True)
# to clear s3 cache run the following command
# prog.remove_cache()

# To report the progress we will split running CRbeam into 10 parts
n_steps = 10
# Initialize multistep simulation
data_exists = not prog.start_multistage_run(
    overwrite_local_cache=True, n_steps=n_steps
)
proceed = not data_exists

if proceed:
    for step in range(n_steps):
        pr.report_progress(
            stage="Running Simulation", progress=100 * step // n_steps
        )
        proceed = prog.simulation_step()
        # todo: report progress using rest API
    pr.report_progress(stage="Running Simulation", progress=100)

assert not proceed, "must be completed before this cell"
if not data_exists:
    prog.end_multistep_run()

def adjust_weights(mc_file, power):
    # converting weights to mimic required injection spectrum power
    header = ""
    with open(mc_file, "rt") as lines:
        for line in lines:
            if len(line) > 0 and line[0] == "#":
                header += line[1:].strip() + "\n"
            else:
                break
    weight_col = 2
    E_src_col = 12
    data = np.loadtxt(mc_file)
    weight = data[:, weight_col - 1]
    E_src = data[:, E_src_col - 1]
    orig_power = (
        prog.power_law
    )  # CRBeam is always called with fixed power=1 to optimize cache usage
    weight *= np.power(E_src / Emax, -(power - orig_power))
    output_file = f"{mc_file}_p{power}"
    np.savetxt(output_file, data, header=header.strip(), fmt="%.6g")
    return output_file

# this code will not execute program since files already exist on server
# prog.run(force_overwrite=False)

# one can upload cache explicitely by call
# prog.upload_cache()

pr.report_progress(stage="Building Plots and Tables", progress=0)

print(prog.output_path)
get_ipython().system("ls {prog.output_path}")   # noqa: F821

mc_file = prog.output_path + "/photon"

# how the data looks like
get_ipython().system("cat {mc_file} | awk 'NR<=5'")   # noqa: F821

if Emax != EminSource:
    mc_file = adjust_weights(mc_file, Gamma)

# how the data looks like
get_ipython().system("cat {mc_file} | awk 'NR<=5'")   # noqa: F821

#! . ./makeSpecE2.sh {mc_file}
# the code above will not work on MMODA as of Sep 30 2023
# here is workaround
subprocess.run(["bash", "makeSpecE2.sh", mc_file])

# rotating the beam

if EGMF > 0:
    from mc_rotate import mc_rotate

    mc_rotated_file = mc_rotate(
        mc_file, jet_half_size, jet_direction, psf_deg=psf
    )
else:
    mc_rotated_file = mc_file
mc_rotated_file

# calculating the energy spectrum
#! . ./makeSpecE2.sh {mc_rotated_file}
subprocess.run(["bash", "makeSpecE2.sh", mc_rotated_file])

# how the rotated data looks like
get_ipython().system("cat {mc_rotated_file} | awk 'NR<=5'")   # noqa: F821

# reading the source distance in Mpc from the data file
T_Mpc = lc.get_distance_Mpc(mc_file)
T_Mpc

# building the spectrum figure for total flux
spec_file = mc_file + ".spec"
spec_fig = plot_spec(
    spec_file,
    title="spectrum",
    Emin=Emin,
    Emax=Emax,
    ext="png",
    show=False,
    logscale=True,
)
spec_image = PictureProduct.from_file(spec_fig)

# building the spectrum figure for the flux within PSF
spec_rotated_file = mc_rotated_file + ".spec"
spec_rotated_fig = plot_spec(
    spec_rotated_file,
    title="spectrum",
    Emin=Emin,
    Emax=Emax,
    ext="png",
    show=False,
    logscale=True,
)
spec_rotated_image = PictureProduct.from_file(spec_rotated_fig)

spec_rotated_file

lc_params = dict(
    logscale=True,
    max_t=-99,  # 8760000,
    n_points=30,
    psf=psf,
    min_n_particles=10,
    cut_0=True,
)

# building the light curve figure
if EGMF > 0:
    light_curve_fig = lc.make_plot(mc_rotated_file, **lc_params)
    light_curve_image = PictureProduct.from_file(light_curve_fig)
else:
    # to avoid possible problems with absent figure
    light_curve_image = PictureProduct.from_file(spec_fig)

l_curve = None
if EGMF > 0:
    delay, weights = lc.get_counts(mc_rotated_file, **lc_params)
    t, f, N = lc.light_curve(delay, weights, **lc_params)
    l_curve = LightCurveDataProduct.from_arrays(
        times=t * 3600.0,  # t is in hours
        fluxes=f,
        errors=f / np.sqrt(N),
        time_format="unix",
    )

def convert_to_ICRS(phi: np.array, theta: np.array, source: SkyCoord):
    # prime system has z-axes pointing the source centered at observer
    # it is obtained by rotation of source system by 180 degrees about x-axis
    # prime system coords have suffix "prime" (')
    # icrs system coords has empty suffix
    # TODO: add param ic_jet_plane_direction: SkyCoord - gives plane which contains jet
    # by definition in prime system the jet is in y'-z' plane and z' axes points from the source towards the observer
    # for simplicity we will first assume that prime frame is oriented as ICRS frame (use SkyOffsetFrame with rotation=0)

    # coordinates of unit direction vector in prime system

    # direction of event arrival at prime system
    z_prime = np.cos(
        theta / 180.0 * np.pi
    )  # (-1)*(-1) = 1 (rotating system and tacking oposite vector)
    r_xy_prime = np.sin(theta / 180.0 * np.pi)
    x_prime = -r_xy_prime * np.cos(
        phi / 180.0 * np.pi
    )  # (1)*(-1) = -1 (rotating system and tacking oposite vector)
    y_prime = r_xy_prime * np.sin(
        phi / 180.0 * np.pi
    )  # (-1)*(-1) = 1 (rotating system and tacking oposite vector)

    print("x',y',z' =", x_prime, y_prime, z_prime)

    # angle between z and z' axes
    delta1 = 90 * u.deg - source.dec

    print("source:", source.cartesian)
    print("delta1 = ", delta1)

    # rotating about y-axes

    x1 = x_prime * np.cos(delta1) + z_prime * np.sin(delta1)
    y1 = y_prime
    z1 = -x_prime * np.sin(delta1) + z_prime * np.cos(delta1)

    print("step 1: x,y,z =", x1, y1, z1)

    # rotation to -RA about z axis
    delta2 = source.ra
    x = x1 * np.cos(delta2) - y1 * np.sin(delta2)
    y = x1 * np.sin(delta2) + y1 * np.cos(delta2)
    z = z1

    print("x,y,z =", x, y, z)

    event_coords = SkyCoord(
        x=x, y=y, z=z, frame="icrs", representation_type="cartesian"
    ).spherical
    # print(event_coords.spherical)
    return event_coords

def LoadCubeTemplate(
    mc_file: str,
    source: SkyCoord,
    Emax=1e15,
    Emin=1e9,
    bins_per_decade=20,
    binsz=0.02,
    theta_mult=1,
    sec_component_mult=1,
    remove_straight=False,
    symmetric_split=0,
):
    pass

    mc_data = np.loadtxt(mc_file)
    E = mc_data[:, 0]
    print(f"dataset energy: {np.min(E)} <= E/eV <= {np.max(E)}")
    Theta = mc_data[:, 2] * theta_mult  # theta_mult is used for debug only
    print(f"dataset Theta: {np.min(Theta)} <= Theta <= {np.max(Theta)}")
    # idx = np.where((E >= Emin) & (E<=Emax) & (Theta<theta_max))[0]
    print(f"Filters: {Emin} <= E/eV <= {Emax}")
    condition = (E >= Emin) & (E <= Emax)
    if remove_straight:
        condition &= Theta > 0

    idx = np.where(condition)[0]

    print("filtered dataset length:", len(idx))

    E = E[idx]
    Theta = Theta[idx]
    w = mc_data[:, 1][idx]
    Phi = mc_data[:, 3][idx]
    # t = mc_data[:,5][idx]
    # t *= time_scale_hours
    if sec_component_mult != 1:  # sec_component_mult is used for debug only
        print("WARNING: Sec component multiplicity:", sec_component_mult)
        w[Theta > 0] *= sec_component_mult

    if len(idx) > 0:
        #     print(f'{np.min(t)} <= t/h <= {np.max(t)}')
        print(f"{np.min(E)} <= E <= {np.max(E)}")

    energy_axis_name = "energy_true"

    energy_axis = MapAxis.from_energy_bounds(
        Emin * u.eV,
        Emax * u.eV,
        nbin=bins_per_decade,
        name=energy_axis_name,
        per_decade=True,
    )

    m_cube = Map.create(
        binsz=binsz,
        width=(window_size_RA * u.deg, window_size_DEC * u.deg),
        frame="icrs",
        axes=[energy_axis],
        skydir=SkyCoord(source),
    )

    print(m_cube.geom)

    if len(idx) > 0:
        if symmetric_split > 1:
            Theta = np.repeat(Theta, symmetric_split)
            E = np.repeat(E, symmetric_split)
            w = np.repeat(w, symmetric_split) / symmetric_split
            Phi = np.random.random(len(Theta)) * 360

        sc = convert_to_ICRS(Phi, Theta, source)
        m_cube.fill_by_coord(
            {"lat": sc.lat, "lon": sc.lon, energy_axis_name: E * u.eV},
            weights=w,
        )

    return m_cube

def LoadTemplate4D(
    mc_file: str,
    source: SkyCoord,
    redshift,
    Emax=1e15,
    Emin=1e9,
    bins_per_decade=20,
    binsz=0.02,
    theta_mult=1,
    max_time_quantile=1.0,
    n_time_bins=100,
    symmetric_split=0,
):
    from astropy.cosmology import FlatLambdaCDM

    cosmo = FlatLambdaCDM(
        H0=67.8 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=(1 - 0.692)
    )
    time_scale_hours = float(
        ((cosmo.age(0) - cosmo.age(z_start)) / u.h).decompose()
    )

    mc_data = np.loadtxt(mc_file)
    E = mc_data[:, 0]

    print("dataset length:", len(mc_data))
    print(f"dataset energy: {np.min(E)} <= E/eV <= {np.max(E)}")
    Theta = mc_data[:, 2] * theta_mult  # theta_mult is used for debug only
    print(f"dataset Theta: {np.min(Theta)} <= Theta <= {np.max(Theta)}")
    # idx = np.where((E >= Emin) & (E<=Emax) & (Theta<theta_max))[0]

    idx = np.where((E >= Emin) & (E <= Emax))[0]

    print(f"Filters: {Emin} <= E/eV <= {Emax}")
    print("filtered dataset length:", len(idx))

    E = E[idx]
    Theta = Theta[idx]
    w = mc_data[:, 1][idx]
    Phi = mc_data[:, 3][idx]
    t = mc_data[:, 5][idx]
    t *= time_scale_hours
    min_time = np.min(t)
    if max_time_quantile >= 1.0:
        max_time = np.max(t)
    else:
        max_time = np.quantile(t, max_time_quantile)

    if max_time == min_time:
        max_time = min_time + 1e-10

    idx = np.where(t <= max_time)[0]
    print(f"Filters: {Emin} <= E/eV <= {Emax}, t/h <= {max_time}")
    print("filtered dataset length:", len(idx))

    E = E[idx]
    Theta = Theta[idx]
    w = w[idx]
    Phi = Phi[idx]
    t = t[idx]

    energy_axis_name = "energy_true"

    energy_axis = MapAxis.from_energy_bounds(
        Emin * u.eV,
        Emax * u.eV,
        nbin=bins_per_decade,
        name=energy_axis_name,
        per_decade=True,
    )

    from astropy.time import Time

    reference_time = Time(0, format="unix")

    delta = max_time - min_time
    (min_time + delta * np.linspace(0, 1, n_time_bins + 1)) * u.h
    # time_axis = TimeMapAxis.from_time_edges(
    #     time_min=time_edges[:-1],
    #     time_max=time_edges[1:],
    #     interp="lin",
    #     unit=u.h,
    #     name='time',
    #     reference_time=reference_time
    #  )

    time_axis = MapAxis.from_bounds(
        min_time, max_time, n_time_bins, name="time", unit=u.h
    )

    # time_axis = TimeMapAxis(time_edges[:-1], time_edges[1:], reference_time, name='time', interp='lin')

    # time_axis = TimeMapAxis.from_time_bounds(min_time, max_time, n_time_bins, unit=u.h, name='time')

    map4d = Map.create(
        binsz=binsz,
        width=(window_size_RA * u.deg, window_size_DEC * u.deg),
        frame="icrs",
        axes=[energy_axis, time_axis],
        skydir=SkyCoord(source),
    )

    print(map4d.geom)

    if len(idx) > 0:
        if symmetric_split > 1:
            Theta = np.repeat(Theta, symmetric_split)
            E = np.repeat(E, symmetric_split)
            t = np.repeat(t, symmetric_split)
            w = np.repeat(w, symmetric_split) / symmetric_split
            Phi = np.random.random(len(Theta)) * 360

        sc = convert_to_ICRS(Phi, Theta, source)
        map4d.fill_by_coord(
            {
                "lat": sc.lat,
                "lon": sc.lon,
                energy_axis_name: E * u.eV,
                "time": t * u.h,
            },
            weights=w,
        )

    return map4d

bins_per_decade = 20

symmetric_split = 0 if jet_direction != 0 else 100

map3d = LoadCubeTemplate(
    mc_rotated_file,
    source=source,
    Emin=Emin,
    Emax=Emax,
    bins_per_decade=bins_per_decade,
    symmetric_split=symmetric_split,
)

map3d.write("map3d.fits", overwrite=True)

map4d = LoadTemplate4D(
    mc_rotated_file,
    source=source,
    redshift=z_start,
    Emin=Emin,
    Emax=Emax,
    bins_per_decade=bins_per_decade,
    symmetric_split=symmetric_split,
)
map4d.write("map4d.fits", overwrite=True)

d = np.genfromtxt(spec_file)
ee = d[:, 0]
ff = d[:, 1]
ff_err = ff / np.sqrt(d[:, 2])
plt.errorbar(ee, ff, yerr=ff_err, label="total flux")
d = np.genfromtxt(spec_rotated_file)
ee1 = d[:, 0]
ff1 = d[:, 1]
ff_err1 = ff1 / np.sqrt(d[:, 2])
plt.errorbar(ee1, ff1, yerr=ff_err1, label="PSF flux")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Energy, eV")
plt.ylabel("$E^2dN/dE$, arbitrary units")
plt.legend(loc="lower right")
plt.savefig("Spectrum.png", format="png", bbox_inches="tight")

bin_image = PictureProduct.from_file("Spectrum.png")

data = [ee, ff, ff_err]
names = ("Energy", "Total_flux", "Total_flux_err")
table = ODAAstropyTable(Table(data, names=names))
data = [ee, ff, ff_err]
names = ("Energy", "PSF_flux", "PSF_flux_err")
table1 = ODAAstropyTable(Table(data, names=names))

flux_unit = u.eV / u.cm / u.cm / u.s / u.sr
spec = Spectrum1D(
    spectral_axis=ee * u.eV,
    flux=ff * flux_unit,
    uncertainty=StdDevUncertainty(ff_err),
)
spec.write("spec.fits", overwrite=True)
dp = NumpyDataProduct.from_fits_file("spec.fits")
spec_rotated = Spectrum1D(
    spectral_axis=ee1 * u.eV,
    flux=ff1 * flux_unit,
    uncertainty=StdDevUncertainty(ff_err1),
)
spec_rotated.write("spec_rotated.fits", overwrite=True)
dp_rotated = NumpyDataProduct.from_fits_file("spec_rotated.fits")
map3d = NumpyDataProduct.from_fits_file("map3d.fits")
map4d = NumpyDataProduct.from_fits_file("map4d.fits")

pr.report_progress(stage="Building Plots and Tables", progress=100)

spectrum_png = bin_image  # http://odahub.io/ontology#ODAPictureProduct
light_curve_png = (
    light_curve_image  # http://odahub.io/ontology#ODAPictureProduct
)
total_spectrum_table = table  # http://odahub.io/ontology#ODAAstropyTable
psf_spectrum_table = table1  # http://odahub.io/ontology#ODAAstropyTable
lc_result = l_curve  # http://odahub.io/ontology#LightCurve
spectrum = dp  # https://odahub.io/ontology/#Spectrum
spectrum_rotated = dp_rotated  # https://odahub.io/ontology/#Spectrum
map3d = map3d  # https://odahub.io/ontology/#Spectrum
map4d = map4d  # https://odahub.io/ontology/#Spectrum

# spec_png = spec_image # http://odahub.io/ontology#ODAPictureProduct
# spec_rotated_png = spec_rotated_image # http://odahub.io/ontology#ODAPictureProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    ("out_CRbeam_spectrum_png", "spectrum_png_galaxy.output", spectrum_png)
)
_oda_outs.append(
    (
        "out_CRbeam_light_curve_png",
        "light_curve_png_galaxy.output",
        light_curve_png,
    )
)
_oda_outs.append(
    (
        "out_CRbeam_total_spectrum_table",
        "total_spectrum_table_galaxy.output",
        total_spectrum_table,
    )
)
_oda_outs.append(
    (
        "out_CRbeam_psf_spectrum_table",
        "psf_spectrum_table_galaxy.output",
        psf_spectrum_table,
    )
)
_oda_outs.append(
    ("out_CRbeam_lc_result", "lc_result_galaxy.output", lc_result)
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
_simple_outs = []
_simple_outs.append(
    ("out_CRbeam_spectrum", "spectrum_galaxy.output", spectrum)
)
_simple_outs.append(
    (
        "out_CRbeam_spectrum_rotated",
        "spectrum_rotated_galaxy.output",
        spectrum_rotated,
    )
)
_simple_outs.append(("out_CRbeam_map3d", "map3d_galaxy.output", map3d))
_simple_outs.append(("out_CRbeam_map4d", "map4d_galaxy.output", map4d))
_numpy_available = True

for _outn, _outfn, _outv in _simple_outs:
    _galaxy_outfile_name = os.path.join(_galaxy_wd, _outfn)
    if isinstance(_outv, str) and os.path.isfile(_outv):
        shutil.move(_outv, _galaxy_outfile_name)
        _galaxy_meta_data[_outn] = {"ext": "_sniff_"}
    elif _numpy_available and isinstance(_outv, np.ndarray):
        with open(_galaxy_outfile_name, "wb") as fd:
            np.savez(fd, _outv)
        _galaxy_meta_data[_outn] = {"ext": "npz"}
    else:
        with open(_galaxy_outfile_name, "w") as fd:
            json.dump(_outv, fd)
        _galaxy_meta_data[_outn] = {"ext": "expression.json"}

with open(os.path.join(_galaxy_wd, "galaxy.json"), "w") as fd:
    json.dump(_galaxy_meta_data, fd)
print("*** Job finished successfully ***")
