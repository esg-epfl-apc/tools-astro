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
radius = 2.0  # http://odahub.io/ontology#AngleDegrees ; oda:label "Size of the Region-Of-Interest (ROI)"
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
from astropy.coordinates import SkyCoord
from astropy.time import Time
from gammapy.data import DataStore
from gammapy.datasets import (
    MapDataset,
)
from gammapy.makers import FoVBackgroundMaker, MapDatasetMaker, SafeMaskMaker
from gammapy.maps import MapAxis, WcsGeom
from gammapy.utils.regions import CircleSkyRegion
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

# ## View source image

# ### Preparing reduced datasets geometry
#
# Now we define a reference geometry for our analysis, We choose a WCS
# based geometry with a binsize of 0.02 deg and also define an energy
# axis:

energy_axis = MapAxis.from_energy_bounds(1.0, 10.0, 4, unit="TeV")

geom = WcsGeom.create(
    skydir=source,
    binsz=0.01 * radius,
    width=(radius, radius),
    frame="icrs",
    proj="CAR",
    axes=[energy_axis],
)

# Reduced IRFs are defined in true energy (i.e. not measured energy).
energy_axis_true = MapAxis.from_energy_bounds(
    Emin, Emax, 10, unit="TeV", name="energy_true"
)

# Now we can define the target dataset with this geometry.

stacked = MapDataset.create(
    geom=geom, energy_axis_true=energy_axis_true, name="Mkn-501-stacked"
)

# ### Data reduction
#
# Create the maker classes to be used

offset_max = 1.0 * u.deg
maker = MapDatasetMaker()
maker_safe_mask = SafeMaskMaker(
    methods=["offset-max", "aeff-max"], offset_max=offset_max
)

circle = CircleSkyRegion(center=source, radius=0.2 * u.deg)

exclusion_mask = ~geom.region_mask(regions=[circle])
maker_fov = FoVBackgroundMaker(method="fit", exclusion_mask=exclusion_mask)

# Perform the data reduction loop

for obs in observations:
    # First a cutout of the target map is produced
    cutout = stacked.cutout(
        obs.get_pointing_icrs(obs.tmid),
        width=2 * offset_max,
        name=f"obs-{obs.obs_id}",
    )
    # A MapDataset is filled in this cutout geometry
    dataset = maker.run(cutout, obs)
    # The data quality cut is applied
    dataset = maker_safe_mask.run(dataset, obs)
    # fit background model
    dataset = maker_fov.run(dataset)
    print(
        f"Background norm obs {obs.obs_id}: {dataset.background_model.spectral_model.norm.value:.2f}"
    )
    # The resulting dataset cutout is stacked onto the final one
    stacked.stack(dataset)

print(stacked)

# ### Inspect image summed over energy

# stacked.counts.sum_over_axes().smooth(0.05 * u.deg).plot(stretch="sqrt", add_cbar=True)
stacked.counts.sum_over_axes().plot(stretch="sqrt", add_cbar=True)
plt.savefig("image.png")

from oda_api.data_products import PictureProduct

output_source_image = PictureProduct.from_file(os.getcwd() + "/image.png")

source_image = output_source_image  # oda:ODAPictureProduct
observations = (
    output_observations_table  # http://odahub.io/ontology#ODAAstropyTable
)

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    ("out_image_source_image", "source_image_galaxy.output", source_image)
)
_oda_outs.append(
    ("out_image_observations", "observations_galaxy.output", observations)
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
