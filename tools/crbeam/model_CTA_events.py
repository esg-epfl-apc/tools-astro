#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from gammapy.data import (
    FixedPointingInfo,
    Observation,
    PointingMode,
    observatory_locations,
)
from gammapy.datasets import MapDataset, MapDatasetEventSampler
from gammapy.irf import load_irf_dict_from_file
from gammapy.makers import MapDatasetMaker
from gammapy.maps import Map, MapAxis, WcsGeom
from gammapy.modeling.models import (
    FoVBackgroundModel,
    Models,
    SkyModel,
    TemplateSpatialModel,
)
from oda_api.data_products import NumpyDataProduct, PictureProduct
from oda_api.json import CustomJSONEncoder
from utils import find_redshift

src_name = "Mrk 501"  # http://odahub.io/ontology#AstrophysicalObject
z_start = 0  # http://odahub.io/ontology#Float
Npart = 5000  # http://odahub.io/ontology#Integer ; oda:lower_limit 1 ; oda:upper_limit 100000
particle_type = "gamma"  # http://odahub.io/ontology#String ; oda:allowed_value "gamma","electron","proton"
Emax = 30  # http://odahub.io/ontology#Energy_TeV
Emin = 0.1  # http://odahub.io/ontology#Energy_TeV
EminSource = 1.0  # http://odahub.io/ontology#Energy_TeV
Gamma = 2.0  # http://odahub.io/ontology#Float
EGMF_fG = 10  # http://odahub.io/ontology#Float
lmaxEGMF_Mpc = 5  # http://odahub.io/ontology#Float
jet_half_size = 5.0  # http://odahub.io/ontology#degree
jet_direction = 0.0  # http://odahub.io/ontology#degree
window_size_RA = 2.0  # http://odahub.io/ontology#degree
window_size_DEC = 1.0  # http://odahub.io/ontology#degree
livetime = 0.1  # http://odahub.io/ontology#TimeIntervalDays
EBL = "Franceschini 2017"  # http://odahub.io/ontology#String ; oda:allowed_value "Franceschini 2017","Stecker 2016 lower limit","Stecker 2016 upper limit","Inoue 2012 Baseline","Inoue 2012 lower limit","Inoue 2012 upper limit"

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

livetime = livetime * u.hr * 24

if z_start <= 0:
    z_start = find_redshift(src_name)

source = SkyCoord.from_name(src_name, frame="icrs", parse=False, cache=True)

z_start, source

import subprocess

import numpy as np
from crbeam import CRbeam
from matplotlib import pyplot as plt

EGMF = EGMF_fG * 1e-15

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
    emax=Emax * 1e12,
    emin=Emin * 1e12,
    emin_source=EminSource * 1e12,
    EGMF=EGMF,
    lmaxEGMF=lmaxEGMF_Mpc,
    background=background,
)
cmd = prog.command
cmd

n_steps = 10
# Initialize multistep simulation
data_exists = not prog.start_multistage_run(
    overwrite_local_cache=True, n_steps=n_steps
)
proceed = not data_exists

if proceed:
    for step in range(n_steps):
        print(f"running simulation {100 * step // n_steps}%", prog.output_dir)
        proceed = prog.simulation_step()
        # todo: report progress using rest API
    print(f"running simulation 100%", prog.output_dir)

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

mc_file = prog.output_path + "/photon"
if Emax != EminSource:
    mc_file = adjust_weights(mc_file, Gamma)

# rotating the beam

if EGMF > 0:
    from mc_rotate import mc_rotate

    mc_rotated_file = mc_rotate(mc_file, jet_half_size, jet_direction)
else:
    mc_rotated_file = mc_file
mc_rotated_file

selected_n_bins_per_decade = 20  # n bins per decade
n_events_reduction_factor = 1  # suppress flux factor
n_models_suppress_factor = (
    1  # exclude only every n_models_suppress_factor-th energy bin
)
theta_mult = 1  # manually increase deflection angle
max_rel_energy_error = 3

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
    redshift,
    Emax=1e3,
    Emin=1e-3,
    bins_per_decade=20,
    binsz=0.02,
    theta_mult=1,
):
    from astropy.cosmology import FlatLambdaCDM

    cosmo = FlatLambdaCDM(
        H0=67.8 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=(1 - 0.692)
    )
    time_scale_hours = float(
        ((cosmo.age(0) - cosmo.age(z_start)) / u.h).decompose()
    )

    mc_data = np.loadtxt(mc_file)
    E = mc_data[:, 0] * 1e-12
    print(f"dataset energy: {np.min(E)} <= E <= {np.max(E)}")
    Theta = mc_data[:, 2] * theta_mult  # theta_mult is used for debug only
    print(f"dataset Theta: {np.min(Theta)} <= Theta <= {np.max(Theta)}")
    # idx = np.where((E >= Emin) & (E<=Emax) & (Theta<theta_max))[0]
    print(f"Filters: {Emin} <= E <= {Emax}")
    idx = np.where((E >= Emin) & (E <= Emax))[0]

    print("filtered dataset length:", len(idx))

    E = E[idx]
    Theta = Theta[idx]
    w = mc_data[:, 1][idx]
    Phi = mc_data[:, 3][idx]
    t = mc_data[:, 5][idx]
    t *= time_scale_hours

    if len(idx) > 0:
        print(f"{np.min(t)} <= t/h <= {np.max(t)}")
        print(f"{np.min(E)} <= E/TeV <= {np.max(E)}")

    energy_axis = MapAxis.from_energy_bounds(
        Emin * u.TeV,
        Emax * u.TeV,
        nbin=bins_per_decade,
        name="energy",
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
        sc = convert_to_ICRS(Phi, Theta, source)
        m_cube.fill_by_coord(
            {"lat": sc.lat, "lon": sc.lon, "energy": E * u.TeV}, weights=w
        )

    return m_cube

source

# telescope is pointing at a fixed position in ICRS for the observation
pointing = FixedPointingInfo(fixed_icrs=source, mode=PointingMode.POINTING)

location = observatory_locations["cta_south"]

# irfs = load_irf_dict_from_file(path / irf_filename)
filename = "data/Prod5-North-20deg-AverageAz-4LSTs09MSTs.180000s-v0.1.fits.gz"
irfs = load_irf_dict_from_file(filename)

observation = Observation.create(
    obs_id=1001,
    pointing=pointing,
    livetime=livetime,
    irfs=irfs,
    location=location,
)
print(observation)

energy_axis = MapAxis.from_energy_bounds(
    max_rel_energy_error * Emin * u.TeV,
    Emax * u.TeV,
    nbin=selected_n_bins_per_decade,
    per_decade=True,
)
energy_axis_true = MapAxis.from_energy_bounds(
    Emin * u.TeV,
    max_rel_energy_error * Emax * u.TeV,
    nbin=selected_n_bins_per_decade,
    per_decade=True,
    name="energy_true",
)
migra_axis = MapAxis.from_bounds(
    0.5, 2, nbin=150, node_type="edges", name="migra"
)

geom = WcsGeom.create(
    skydir=pointing.fixed_icrs,
    width=(2, 2),
    binsz=0.02,
    frame="icrs",
    axes=[energy_axis],
)

cube_map = LoadCubeTemplate(
    mc_rotated_file,
    source=source,
    redshift=z_start,
    theta_mult=theta_mult,
    Emin=Emin,
    Emax=Emax,
    bins_per_decade=selected_n_bins_per_decade,
)
cube_map.write("3d.fits", overwrite=True)

mid_energy = cube_map.data.shape[0] // 2
map0 = cube_map.slice_by_idx({"energy": mid_energy})
map0.plot()
plt.show()
# print(map0)

empty = MapDataset.create(
    geom,
    energy_axis_true=energy_axis_true,
    migra_axis=migra_axis,
    name="my-dataset",
)
maker = MapDatasetMaker(selection=["exposure", "background", "psf", "edisp"])
dataset = maker.run(empty, observation)

Path("event_sampling").mkdir(exist_ok=True)
dataset.write("./event_sampling/dataset.fits", overwrite=True)

def GetBinSpectralModel(E, bins_per_decade=20, norm=1):
    amplitude = 1e-12 * u.Unit("cm-2 s-1") * norm
    from gammapy.modeling.models import GaussianSpectralModel

    sigma = (10 ** (1 / bins_per_decade) - 1) * E
    return GaussianSpectralModel(mean=E, sigma=sigma, amplitude=amplitude)

spec = cube_map.get_spectrum()
spec

spec.plot()  # this plot shows dN/dE * E (below we check this)

spec.data.shape, spec.data[spec.data.shape[0] // 2, 0, 0]

# ### Let us build spectrum manually and compare with output of cube_map.get_spectrum()

subprocess.run(["bash", "makeSpecE2.sh", mc_rotated_file])

spec_file = mc_rotated_file + ".spec"
spec_data = np.loadtxt(spec_file)
E = spec_data[:, 0]
fluxE = spec_data[:, 1] / E  # convert dN/dE * E^2 to dN/dE * E
plt.scatter(np.log10(E) - 12, np.log10(fluxE))
plt.xlabel("Log(E/TeV)")
plt.ylabel("Log(dN/dE * E)")

energy_bins = cube_map.geom.axes["energy"].center
len(energy_bins), float(np.max(energy_bins) / u.TeV)

int_bin_flux = (
    spec.data.flatten()
)  # we don't have to multiply by energy_bins /u.TeV since spectrum is in already multiplied by E (see above)
int_bin_flux /= (
    Npart
    / 200000
    * np.max(int_bin_flux)
    * n_events_reduction_factor
    * 20
    / len(energy_bins)
)  # roughly 100 events
int_bin_flux

bin_models = []
for i, (flux, E) in enumerate(zip(int_bin_flux, energy_bins)):
    if n_models_suppress_factor > 1 and i % n_models_suppress_factor != 0:
        continue
    if flux == 0:
        continue
    spectral_model_delta = GetBinSpectralModel(
        E, norm=flux
    )  # normalizing here
    spacial_template_model = TemplateSpatialModel(
        cube_map.slice_by_idx({"energy": i}),
        filename=f"cube_bin{i}.fit",
        normalize=True,
    )
    sky_bin_model = SkyModel(
        spectral_model=spectral_model_delta,
        spatial_model=spacial_template_model,
        name=f"bin_{i}",
    )
    bin_models.append(sky_bin_model)

bkg_model = FoVBackgroundModel(dataset_name="my-dataset")
models = Models(bin_models + [bkg_model])
file_model = "./event_sampling/cube.yaml"
models.write(file_model, overwrite=True)

dataset.models = models
print(dataset.models)

sampler = MapDatasetEventSampler(random_state=0)
events = sampler.run(dataset, observation)

print(f"Source events: {(events.table['MC_ID'] > 0).sum()}")
print(f"Background events: {(events.table['MC_ID'] == 0).sum()}")

for i in range(1, len(bin_models) + 1):
    n = (events.table["MC_ID"] == i).sum()
    if n > 1:
        print(f"\tmodel {i}: {n} events")

events.select_offset([0, 1] * u.deg).peek()
plt.savefig("events.png", format="png", bbox_inches="tight")
plt.show()

events.table

print(f"Save events ...")
primary_hdu = fits.PrimaryHDU()
hdu_evt = fits.BinTableHDU(events.table)
hdu_gti = fits.BinTableHDU(dataset.gti.table, name="GTI")
hdu_all = fits.HDUList([primary_hdu, hdu_evt, hdu_gti])
hdu_all.writeto(f"./events.fits", overwrite=True)
####################
hdul = fits.open("events.fits")

events.table.meta

events_fits = NumpyDataProduct.from_fits_file("events.fits")
events_summary = PictureProduct.from_file("events.png")

events_summary_png = (
    events_summary  # http://odahub.io/ontology#ODAPictureProduct
)
events_fits_fits = events_fits  # https://odahub.io/ontology/#Spectrum

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_model_CTA_events_events_summary_png",
        "events_summary_png_galaxy.output",
        events_summary_png,
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
_simple_outs = []
_simple_outs.append(
    (
        "out_model_CTA_events_events_fits_fits",
        "events_fits_fits_galaxy.output",
        events_fits_fits,
    )
)
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
