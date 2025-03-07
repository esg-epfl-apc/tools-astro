#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

token = ""  # http://odahub.io/ontology#LongString  ; oda:label "token used for data access"
src_name = "Mrk 501"  # http://odahub.io/ontology#AstrophysicalObject
T1 = "2028-01-01T00:00:00.0"  # http://odahub.io/ontology#StartTime
T2 = "2028-12-31T23:59:59.0"  # http://odahub.io/ontology#EndTime
RA = 0  # http://odahub.io/ontology#PointOfInterestRA
DEC = 0  # http://odahub.io/ontology#PointOfInterestDEC
radius = 2.0  # http://odahub.io/ontology#AngleDegrees ; oda:label "Size of the Region-Of-Interest (ROI)"

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in ["token", "src_name", "T1", "T2", "RA", "DEC", "radius"]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

get_ipython().run_cell_magic(   # noqa: F821
    "bash",
    "",
    "if [ ! -f sdc_setup.py ]\nthen\n    git clone https://gitlab.renkulab.io/astronomy/mmoda/ctao.git tmp_src\n    cp tmp_src/*.sh tmp_src/*.py ./\nfi\n",
)

# ## SDC data access setup

from sdc_setup import setup

setup(token)

get_ipython().run_line_magic("matplotlib", "inline")   # noqa: F821
import logging
import os

# Check package versions
import astropy.units as u
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

# %matplotlib inline
from gammapy.data import DataStore

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
observations = data_store.get_observations(obs_ids)
print(f"Number of selected observations : {len(observations)}")

from astropy.table import Table
from oda_api.data_products import ODAAstropyTable

data = [obs_ids]
names = ("Id",)
output_observations_table = ODAAstropyTable(Table(data, names=names))

observations = (
    output_observations_table  # http://odahub.io/ontology#ODAAstropyTable
)

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_observations_observations",
        "observations_galaxy.output",
        observations,
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
