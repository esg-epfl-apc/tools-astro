#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python

# This script is generated with nb2galaxy

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

# oda:usesResource oda:CRBeamS3 .
# oda:CRBeamS3 a oda:S3 .
# oda:CRBeamS3 oda:resourceBindingEnvVarName "S3_CREDENTIALS" .

src_name = "1ES 1215+303"  # http://odahub.io/ontology#AstrophysicalObject

z_start = 0.13  # http://odahub.io/ontology#Float
Npart = 10000  # http://odahub.io/ontology#Integer ; oda:lower_limit 1 ; oda:upper_limit 100000
particle_type = "gamma"  # http://odahub.io/ontology#String ; oda:allowed_value "gamma","electron","proton"
Emax = 50  # http://odahub.io/ontology#Energy_TeV
Emin = 0.01  # http://odahub.io/ontology#Energy_TeV
EminSource = 0.01  # http://odahub.io/ontology#Energy_TeV
Gamma = 2.0  # http://odahub.io/ontology#Float
EGMF_fG = 100  # http://odahub.io/ontology#Float
lmaxEGMF_Mpc = 5  # http://odahub.io/ontology#Float
jet_half_size = 180.0  # http://odahub.io/ontology#AngleDegrees
jet_direction = 0.0  # http://odahub.io/ontology#AngleDegrees
psf = 180.0  # http://odahub.io/ontology#AngleDegrees
window_size_RA = 4.0  # http://odahub.io/ontology#AngleDegrees
window_size_DEC = 4.0  # http://odahub.io/ontology#AngleDegrees
EBL = "Franceschini 2017"  # http://odahub.io/ontology#String ; oda:allowed_value "Franceschini 2017","Stecker 2016 lower limit","Stecker 2016 upper limit","Inoue 2012 Baseline","Inoue 2012 lower limit","Inoue 2012 upper limit"

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "C_data_product_" in inp_dic.keys():
    inp_pdic = inp_dic["C_data_product_"]
else:
    inp_pdic = inp_dic
src_name = str(inp_pdic["src_name"])
z_start = float(inp_pdic["z_start"])
Npart = int(inp_pdic["Npart"])
particle_type = str(inp_pdic["particle_type"])
Emax = float(inp_pdic["Emax"])
Emin = float(inp_pdic["Emin"])
EminSource = float(inp_pdic["EminSource"])
Gamma = float(inp_pdic["Gamma"])
EGMF_fG = float(inp_pdic["EGMF_fG"])
lmaxEGMF_Mpc = float(inp_pdic["lmaxEGMF_Mpc"])
jet_half_size = float(inp_pdic["jet_half_size"])
jet_direction = float(inp_pdic["jet_direction"])
psf = float(inp_pdic["psf"])
window_size_RA = float(inp_pdic["window_size_RA"])
window_size_DEC = float(inp_pdic["window_size_DEC"])
EBL = str(inp_pdic["EBL"])

import os
import sys

sys.path.append(os.environ.get("BASEDIR", os.getcwd()))

from astropy.coordinates import SkyCoord
from astropy.table import Table
from oda_api.api import ProgressReporter
from oda_api.data_products import ODAAstropyTable
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

import light_curve as lc
import numpy as np
from crbeam import CRbeam

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

if Emax != EminSource:
    mc_file = adjust_weights(mc_file, Gamma)

# rotating the beam

if EGMF > 0:
    from mc_rotate import mc_rotate

    mc_rotated_file = mc_rotate(
        mc_file, jet_half_size, jet_direction, psf_deg=psf
    )
else:
    mc_rotated_file = mc_file
mc_rotated_file

# reading the source distance in Mpc from the data file
T_Mpc = lc.get_distance_Mpc(mc_file)
T_Mpc

from oda_api.data_products import ODAAstropyTable

d = np.genfromtxt(mc_rotated_file, skip_header=3)
# E/eV,weight,Theta,Phi,dT_raw/T	dT_calc/T,z_cascade_production,N_interactions,z_src,E_src/eV,z_prod
names = (
    "E/eV",
    "weight",
    "Theta",
    "Phi",
    "dT_raw/T",
    "dT_calc/T",
    "z_cascade_production",
    "N_interactions",
    "z_src",
    "E_src/eV",
    "z_prod",
)
res = ODAAstropyTable(Table(d, names=names))

pr.report_progress(stage="Building Plots and Tables", progress=100)

Event_file = res  # http://odahub.io/ontology#ODAAstropyTable

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    ("out_Generate_events_Event_file", "Event_file_galaxy.output", Event_file)
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
