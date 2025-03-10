#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

token = ""  # http://odahub.io/ontology#LongString  ; oda:label "token used for data access"
src_name = "Mrk 501"  # http://odahub.io/ontology#AstrophysicalObject  ; oda:label "source name"
T1 = "2028-01-27T00:00:00.0"  # http://odahub.io/ontology#StartTime
T2 = "2028-01-30T23:59:59.0"  # http://odahub.io/ontology#EndTime
RA = 0  # http://odahub.io/ontology#PointOfInterestRA
DEC = 0  # http://odahub.io/ontology#PointOfInterestDEC
Emin = 0.1  # http://odahub.io/ontology#Energy_TeV ; oda:label "minimal energy"
Emax = (
    10.0  # http://odahub.io/ontology#Energy_TeV ; oda:label "maximal energy"
)
radius = 2.0  # http://odahub.io/ontology#AngleDegrees ; oda:label "Size of the Region-Of-Interest (ROI)"
max_observations = 500  # http://odahub.io/ontology#Integer ; oda:label "limit total amount of observations to use"
on_region_radius = 0.1  # http://odahub.io/ontology#AngleDegrees ; oda:label "Size of ON-Region"

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in [
    "token",
    "src_name",
    "T1",
    "T2",
    "RA",
    "DEC",
    "Emin",
    "Emax",
    "radius",
    "max_observations",
    "on_region_radius",
]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

get_ipython().run_cell_magic(   # noqa: F821
    "bash",
    "",
    "if [ ! -f sdc_setup.py ]\nthen\n    git clone https://gitlab.renkulab.io/astronomy/mmoda/ctao.git tmp_src\n    cp tmp_src/*.sh tmp_src/*.py ./\nfi\npip install git+https://github.com/cta-epfl/ctadata.git@54-temp-dir-creation-failure\n",
)

# ## SDC data access setup

from sdc_setup import event_files, load_observation, setup

setup(token)

get_ipython().run_line_magic("matplotlib", "inline")   # noqa: F821
import logging
import os

# Check package versions
import astropy.units as u

# %matplotlib inline
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from gammapy.analysis import Analysis, AnalysisConfig
from gammapy.data import DataStore
from gammapy.estimators import LightCurveEstimator
from gammapy.modeling.models import (
    ExpCutoffPowerLawSpectralModel,
    Models,
    PointSpatialModel,
    SkyModel,
)
from IPython.display import display

log = logging.getLogger(__name__)

GAMMAPY_DATA = os.path.join(os.getcwd(), ".")
os.environ["GAMMAPY_DATA"] = GAMMAPY_DATA
CALDB = os.path.join(os.getcwd(), "IRFS")
os.environ["CALDB"] = "IRFS"

if src_name:
    source = SkyCoord.from_name(
        src_name, frame="icrs", parse=False, cache=True
    )
    RA = source.ra
    DEC = source.dec
else:
    RA = RA * u.deg
    DEC = DEC * u.deg
    source = SkyCoord(ra=RA, dec=DEC)
    src_name = "coord_based"
source

# ### List observations available for the source selected
# We select observations within 2 degrees of the source

T1 = Time(T1, format="isot", scale="utc")
T2 = Time(T2, format="isot", scale="utc")
data_store = DataStore.from_dir(".")
selection = dict(
    type="sky_circle",
    frame="icrs",
    lon=RA,
    lat=DEC,
    radius=radius * u.deg,
)
selected_obs_table = data_store.obs_table.select_observations(selection)
print(f"Number of observations in selected region: {len(selected_obs_table)}")
selected_obs_table = selected_obs_table.select_time_range((T1, T2))[
    :max_observations
]
obs_ids = selected_obs_table["OBS_ID"]
observations = data_store.get_observations(obs_ids)
print(f"Number of selected observations : {len(observations)}")

from oda_api.data_products import ODAAstropyTable

output_observations_table = ODAAstropyTable(selected_obs_table)

selected_obs_table

# ### Loading observations

if len(event_files) == 0:
    print("offline mode: skipping event loading")
else:
    for obs_id in observations.ids:
        load_observation(obs_id)

selected_observations = obs_ids

# ## Running the light curve extraction in 1D

conf_1d = AnalysisConfig()

conf_1d.observations.obs_ids = selected_observations

# We define the datastore containing the data
conf_1d.observations.datastore = "."

# We want a 1D analysis
conf_1d.datasets.type = "1d"

# We want to extract the data by observation and therefore to not stack them
conf_1d.datasets.stack = False

# Here we define the ON region and make sure that PSF leakage is corrected
conf_1d.datasets.on_region = dict(
    frame="icrs",
    lon=RA,
    lat=DEC,
    radius=on_region_radius * u.deg,
)
conf_1d.datasets.containment_correction = True

# Finally we define the energy binning for the spectra
conf_1d.datasets.geom.axes.energy = dict(
    min=Emin * u.TeV, max=Emax * u.TeV, nbins=5
)
conf_1d.datasets.geom.axes.energy_true = dict(
    min=0.5 * Emin * u.TeV, max=2 * Emax * u.TeV, nbins=20
)

# ### Run the 1D data reduction

analysis_1d = Analysis(conf_1d)
analysis_1d.get_observations()
datasets = analysis_1d.get_datasets()

# ### Define the model to be used
# Here we use spectral model, which was fitted before

spatial_model = PointSpatialModel(lon_0=RA, lat_0=DEC, frame="icrs")

source_model_name = src_name.replace(" ", "")
spectral_model = ExpCutoffPowerLawSpectralModel(
    amplitude=1e-12 * u.Unit("cm-2 s-1 TeV-1"),
    index=2,
    lambda_=0.1 * u.Unit("TeV-1"),
    reference=1 * u.TeV,
)

sky_model = SkyModel(
    spatial_model=spatial_model,
    spectral_model=spectral_model,
    name=source_model_name,
)

models = Models([sky_model])

analysis_1d.set_models(models)

lc_maker_1d = LightCurveEstimator(
    energy_edges=[1, 10] * u.TeV, source=source_model_name, reoptimize=False
)
lc_1d = lc_maker_1d.run(analysis_1d.datasets)

print(lc_1d.geom.axes.names)

output_lc = lc_1d.to_table(sed_type="flux", format="lightcurve")

display(output_lc)   # noqa: F821

output_lc_table = ODAAstropyTable(output_lc)

fig, ax = plt.subplots(
    figsize=(8, 6),
    gridspec_kw={"left": 0.16, "bottom": 0.2, "top": 0.98, "right": 0.98},
)
lc_1d.plot(ax=ax, marker="o", label="1D")
# lc_3d.plot(ax=ax, marker="o", label="3D")
plt.legend()
plt.savefig("lightcurve.png")
from oda_api.data_products import PictureProduct

output_lightcurve_image = PictureProduct.from_file("lightcurve.png")

lc_image = output_lightcurve_image  # oda:ODAPictureProduct
lc_table = output_lc_table  # oda:ODAAstropyTable
observations = output_observations_table  # oda:ODAAstropyTable

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    ("out_lightcurve_lc_image", "lc_image_galaxy.output", lc_image)
)
_oda_outs.append(
    ("out_lightcurve_lc_table", "lc_table_galaxy.output", lc_table)
)
_oda_outs.append(
    ("out_lightcurve_observations", "observations_galaxy.output", observations)
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
