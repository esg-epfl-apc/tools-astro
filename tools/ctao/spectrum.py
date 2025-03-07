#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

token = ""  # http://odahub.io/ontology#LongString  ; oda:label "used for data access"
src_name = "Mrk 501"  # http://odahub.io/ontology#AstrophysicalObject
T1 = "2028-01-01T00:00:00.0"  # http://odahub.io/ontology#StartTime
T2 = "2028-12-31T23:59:59.0"  # http://odahub.io/ontology#EndTime
RA = 0  # http://odahub.io/ontology#PointOfInterestRA
DEC = 0  # http://odahub.io/ontology#PointOfInterestDEC
Emin = 0.1  # http://odahub.io/ontology#Energy_TeV ; oda:label "minimal energy"
Emax = (
    10.0  # http://odahub.io/ontology#Energy_TeV ; oda:label "maximal energy"
)
radius = 2.0  # http://odahub.io/ontology#AngleDegrees ;  oda:label "Size of the Region-Of-Interest (ROI)"
on_region_radius = 0.2  # http://odahub.io/ontology#AngleDegrees ; oda:label "ON region radius"
nbin_per_decade = 5  # http://odahub.io/ontology#Integer
max_observations = 50  # http://odahub.io/ontology#Integer ; oda:label "limit total amount of observations to use"

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
    "on_region_radius",
    "nbin_per_decade",
    "max_observations",
]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

get_ipython().run_cell_magic(   # noqa: F821
    "bash",
    "",
    "if [ ! -f sdc_setup.py ]\nthen\n    git clone https://gitlab.renkulab.io/astronomy/mmoda/ctao.git tmp_src\n    cp tmp_src/*.sh tmp_src/*.py ./\nfi\n",
)

# ## SDC data access setup

from sdc_setup import load_observation, setup

setup(token)

get_ipython().run_line_magic("matplotlib", "inline")   # noqa: F821
import logging
import os

# Check package versions
import astropy.units as u

# %matplotlib inline
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.time import Time
from gammapy.data import DataStore
from gammapy.datasets import (
    Datasets,
    FluxPointsDataset,
    SpectrumDataset,
)
from gammapy.estimators import FluxPointsEstimator
from gammapy.estimators.utils import resample_energy_edges
from gammapy.makers import (
    ReflectedRegionsBackgroundMaker,
    SafeMaskMaker,
    SpectrumDatasetMaker,
)
from gammapy.maps import MapAxis, RegionGeom
from gammapy.modeling import Fit
from gammapy.modeling.models import (
    ExpCutoffPowerLawSpectralModel,
    SkyModel,
)
from gammapy.utils.regions import CircleSkyRegion
from IPython.display import display
from regions import CircleSkyRegion

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

selected_obs_table = selected_obs_table.select_time_range((T1, T2))
obs_ids = selected_obs_table["OBS_ID"]
observations = data_store.get_observations(obs_ids[:max_observations])
print(f"Number of selected observations : {len(observations)}")

if len(observations) == 0:
    raise Exception("No observations found")

from astropy.table import Table
from oda_api.data_products import ODAAstropyTable

data = [obs_ids]
names = ("Id",)
output_observations_table = ODAAstropyTable(Table(data, names=names))

# ### Loading observations

for obs_id in observations.ids:
    load_observation(obs_id)

# ### Define Target Region
#
# The next step is to define a signal extraction region, also known as on
# region. In the simplest case this is just a
# [CircleSkyRegion](http://astropy-regions.readthedocs.io/en/latest/api/regions.CircleSkyRegion.html)_.
#
#
#

on_region_radius = Angle(on_region_radius * u.deg)
on_region = CircleSkyRegion(center=source, radius=on_region_radius)

# Exclusion mask can be applied if needed

# ### Run data reduction chain
#
# We begin with the configuration of the maker classes:
#
#
#

energy_axis = MapAxis.from_energy_bounds(
    Emin,
    Emax,
    nbin=nbin_per_decade,
    per_decade=True,
    unit="TeV",
    name="energy",
)
energy_axis_true = MapAxis.from_energy_bounds(
    0.5 * Emin,
    2 * Emax,
    nbin=nbin_per_decade,
    per_decade=True,
    unit="TeV",
    name="energy_true",
)

geom = RegionGeom.create(region=on_region, axes=[energy_axis])
dataset_empty = SpectrumDataset.create(
    geom=geom, energy_axis_true=energy_axis_true
)

dataset_maker = SpectrumDatasetMaker(
    containment_correction=True, selection=["counts", "exposure", "edisp"]
)
# bkg_maker = ReflectedRegionsBackgroundMaker(exclusion_mask=exclusion_mask)
bkg_maker = ReflectedRegionsBackgroundMaker()
safe_mask_masker = SafeMaskMaker(methods=["aeff-max"], aeff_percent=10)

datasets = Datasets()

for obs_id, observation in zip(observations.ids, observations):
    dataset = dataset_maker.run(
        dataset_empty.copy(name=str(obs_id)), observation
    )
    dataset_on_off = bkg_maker.run(dataset, observation)
    dataset_on_off = safe_mask_masker.run(dataset_on_off, observation)
    datasets.append(dataset_on_off)

print(datasets)

# ### Source statistic
#
# Next we’re going to look at the overall source statistics in our signal
# region.
#
#
#

info_table = datasets.info_table(cumulative=True)

display(info_table)   # noqa: F821

from oda_api.data_products import ODAAstropyTable

output_info_table = ODAAstropyTable(info_table)
type(info_table)

# And make the corresponding plots
#
#

fig, (ax_excess, ax_sqrt_ts) = plt.subplots(figsize=(10, 4), ncols=2, nrows=1)
ax_excess.plot(
    info_table["livetime"].to("h"),
    info_table["excess"],
    marker="o",
    ls="none",
)

ax_excess.set_title("Excess")
ax_excess.set_xlabel("Livetime [h]")
ax_excess.set_ylabel("Excess events")

ax_sqrt_ts.plot(
    info_table["livetime"].to("h"),
    info_table["sqrt_ts"],
    marker="o",
    ls="none",
)

ax_sqrt_ts.set_title("Sqrt(TS)")
ax_sqrt_ts.set_xlabel("Livetime [h]")
ax_sqrt_ts.set_ylabel("Sqrt(TS)")
plt.savefig("excess_events.png")
from oda_api.data_products import PictureProduct

output_excess_events_image = PictureProduct.from_file("excess_events.png")

# Finally you can write the extracted datasets to disk using the OGIP
# format (PHA, ARF, RMF, BKG, see
# [here](https://gamma-astro-data-formats.readthedocs.io/en/latest/spectra/ogip/index.html)_
# for details):
#
#
#

# ### Fit spectrum
#
# Now we’ll fit a global model to the spectrum. First we do a joint
# likelihood fit to all observations. If you want to stack the
# observations see below. We will also produce a debug plot in order to
# show how the global fit matches one of the individual observations.
#
#
#

source_model_name = src_name.replace(" ", "")
spectral_model = ExpCutoffPowerLawSpectralModel(
    amplitude=1e-12 * u.Unit("cm-2 s-1 TeV-1"),
    index=2,
    lambda_=0.1 * u.Unit("TeV-1"),
    reference=1 * u.TeV,
)
model = SkyModel(spectral_model=spectral_model, name=source_model_name)

datasets.models = [model]

fit_joint = Fit()
result_joint = fit_joint.run(datasets=datasets)

# we make a copy here to compare it later
model_best_joint = model.copy()

# ### Fit quality and model residuals
#
#
#

# We can access the results dictionary to see if the fit converged:
#
#
#

print(result_joint)

# and check the best-fit parameters
#
#
#

joined_model_params_table = result_joint.models.to_parameters_table()
display(joined_model_params_table)   # noqa: F821
output_joined_model_params = ODAAstropyTable(joined_model_params_table)

# A simple way to inspect the model residuals is using the function
# `~SpectrumDataset.plot_fit()`
#
#
#

ax_spectrum, ax_residuals = datasets[0].plot_fit()
ax_spectrum.set_ylim(0.1, 40)
datasets[0].plot_masks(ax=ax_spectrum)
plt.savefig("counts_fit.png")
from oda_api.data_products import PictureProduct

output_counts_fit_image = PictureProduct.from_file("counts_fit.png")

dataset_stacked = Datasets(datasets).stack_reduce()
energy_edges = resample_energy_edges(
    dataset_stacked, conditions={"sqrt_ts_min": 2}
)

# ### Compute Flux Points
#
# To round up our analysis we can compute flux points by fitting the norm
# of the global model in energy bands.
# We can utilise the `~gammapy.estimators.utils.resample_energy_edges`
# for defining the energy bins in which the minimum number of `sqrt_ts` is 2.
# To do so we first stack the individual datasets, only for obtaining the energies:
#
#
#

# Now we create an instance of the
# `~gammapy.estimators.FluxPointsEstimator`, by passing the dataset and
# the energy binning:
#
#
#

fpe = FluxPointsEstimator(
    energy_edges=energy_edges,
    source=source_model_name,
    selection_optional="all",
)
flux_points = fpe.run(datasets=datasets)

# Here is a the table of the resulting flux points:
#
#
#

flux_points_table = flux_points.to_table(sed_type="dnde", formatted=True)
output_flux_points = ODAAstropyTable(flux_points_table)
display(flux_points_table)   # noqa: F821

# Now we plot the flux points and their likelihood profiles. For the
# plotting of upper limits we choose a threshold of TS < 4.
#
#
#

fig, ax = plt.subplots()
flux_points.plot(ax=ax, sed_type="e2dnde", color="darkorange")
flux_points.plot_ts_profiles(ax=ax, sed_type="e2dnde")
plt.savefig("flux_points.png")
output_flux_points_image = PictureProduct.from_file("flux_points.png")

# The final plot with the best fit model, flux points and residuals can be
# quickly made like this:
#
#
#

flux_points_dataset = FluxPointsDataset(
    data=flux_points, models=model_best_joint
)
flux_points_dataset.plot_fit()
plt.savefig("spec_fit.png")
output_spec_fit_image = PictureProduct.from_file("spec_fit.png")

# ### Stack observations
#
# An alternative approach to fitting the spectrum is stacking all
# observations first and then fitting a model. For this we first stack the
# individual datasets:
#
#
#

dataset_stacked = Datasets(datasets).stack_reduce()

# Again we set the model on the dataset we would like to fit (in this case
# it’s only a single one) and pass it to the `~gammapy.modeling.Fit`
# object:
#
#
#

dataset_stacked.models = model
stacked_fit = Fit()
result_stacked = stacked_fit.run([dataset_stacked])

# Make a copy to compare later
model_best_stacked = model.copy()

print(result_stacked)

# And display the parameter table
#
#

display(model_best_joint.parameters.to_table())   # noqa: F821
stacked_model_params_table = model_best_stacked.parameters.to_table()
display(stacked_model_params_table)   # noqa: F821
output_stacked_model_params = ODAAstropyTable(stacked_model_params_table)

stacked_model_params = (
    output_stacked_model_params  # http://odahub.io/ontology#ODAAstropyTable
)
joined_model_params = (
    output_joined_model_params  # http://odahub.io/ontology#ODAAstropyTable
)
info = output_info_table  # http://odahub.io/ontology#ODAAstropyTable
observations = (
    output_observations_table  # http://odahub.io/ontology#ODAAstropyTable
)
flux_points = output_flux_points  # http://odahub.io/ontology#ODAAstropyTable

counts_fit_image = output_counts_fit_image  # oda:ODAPictureProduct
excess_events_image = output_excess_events_image  # oda:ODAPictureProduct
flux_points_image = output_flux_points_image  # oda:ODAPictureProduct
spec_fit_image = output_spec_fit_image  # oda:ODAPictureProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_spectrum_stacked_model_params",
        "stacked_model_params_galaxy.output",
        stacked_model_params,
    )
)
_oda_outs.append(
    (
        "out_spectrum_joined_model_params",
        "joined_model_params_galaxy.output",
        joined_model_params,
    )
)
_oda_outs.append(("out_spectrum_info", "info_galaxy.output", info))
_oda_outs.append(
    ("out_spectrum_observations", "observations_galaxy.output", observations)
)
_oda_outs.append(
    ("out_spectrum_flux_points", "flux_points_galaxy.output", flux_points)
)
_oda_outs.append(
    (
        "out_spectrum_counts_fit_image",
        "counts_fit_image_galaxy.output",
        counts_fit_image,
    )
)
_oda_outs.append(
    (
        "out_spectrum_excess_events_image",
        "excess_events_image_galaxy.output",
        excess_events_image,
    )
)
_oda_outs.append(
    (
        "out_spectrum_flux_points_image",
        "flux_points_image_galaxy.output",
        flux_points_image,
    )
)
_oda_outs.append(
    (
        "out_spectrum_spec_fit_image",
        "spec_fit_image_galaxy.output",
        spec_fit_image,
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
