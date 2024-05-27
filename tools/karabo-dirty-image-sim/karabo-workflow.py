#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

from __future__ import annotations

import json
import os
import shutil
from typing import get_args

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from karabo.imaging.imager import Imager
from karabo.simulation.interferometer import InterferometerSimulation
from karabo.simulation.observation import Observation
from karabo.simulation.sky_model import SkyModel
from karabo.simulation.telescope import (
    OSKARTelescopesWithoutVersionType,
    Telescope,
)
from karabo.simulation.visibility import Visibility
from oda_api.data_products import PictureProduct
from oda_api.json import CustomJSONEncoder

# # Karabo Workflow
#
# This notebook acts as MMODA-backend to create a dirty-image. Karabo version: 0.24.0
#
# TODOs:
# - Use MMODA-object storage for remote-resources. Otherwise the surveys have to be downloaded through the internet for each call.
# - Probably provide telescope frequency-bands as info? And check user-input accordingly
# - Constrain params to avoid memory and time issues (duration-simulation, freq-range, field-of-view), through ontology and run-time-checks. Explore constraints using front-end testing. Runners should have about 20GB RAM and 2 CPU's?
# - Cleanup before going live by adding `live-workflow` and testing afterwards
#
# The most important information can be seen at `https://odahub.io/docs/guide-development/`. The ontology items can be seen at `https://odahub.io/ontology/`. An important thing to mention is, that MMODA doesn't consider the Dockerfile, and should get stuff done automatically like checkout on correct venv and install from environment.yaml and requirements.txt!

# Render plots inline
get_ipython().run_line_magic("matplotlib", "inline")   # noqa: F821

# - Next cell is tagged as `parameters`!
# - Do not create parameters which depend on inputs of other parameters because of `injected-parameters`!
# - Use the same python-types for type-hints as depicted from ontology (NoneTypes not supported yet).
# - Ontology-time-format: YYYY-MM-DDThh:mm:ss.s
# - MMODA-default-params:
#     - `RA`: http://odahub.io/ontology#PointOfInterestRA
#     - `DEC`: http://odahub.io/ontology#PointOfInterestDEC
#     - `start_time`: http://odahub.io/ontology#StartTime
#     - `end_time`: http://odahub.io/ontology#EndTime
#     - `src_name`: http://odahub.io/ontology#AstrophysicalObject
#         - ignoring `src_name` because RA-DEC are always provided and must align with `src_name` from the front-end
# - Default values for default-params are not relevant, because they're not shown in the front-end as defaults and have to be provided as input.
# - All parameters should have default-values according to MMODA-team.
# - Limits are also not relevant for default-params, because they're already applied, however not for non-default-params.
# - Survey-pattern: '\<survey\> \<low-hz\>-\<high-hz\> Hz, RA\[\<low\>\<high\>\]/DEC\[\<low\>,\<high\>\] deg'
# - Telescope with versions not included yet (not sure how to solve that).

RA = 265.97845833  # http://odahub.io/ontology#PointOfInterestRA
DEC = -29.74516667  # http://odahub.io/ontology#PointOfInterestDEC
start_time = "2017-03-06T13:26:48.0"  # http://odahub.io/ontology#StartTime
end_time = "2017-03-06T15:32:27.0"  # http://odahub.io/ontology#EndTime

survey = "gleam 72e6-231e6 Hz, RA[0.0,360.0]/DEC[-90.0,30.0] deg"  # http://odahub.io/ontology#String; oda:allowed_value 'gleam 72e6-231e6 Hz, RA[0.0,360.0]/DEC[-90.0,30.0] deg', 'mals-dr1v3 902e6-1644e6 Hz, RA[0.0,360.0]/DEC[-78.80,32.35] deg', 'mightee-l1 1304e6-1375e6 Hz, RA[149.40,150.84]/DEC[1.49,2.93] deg'; http://purl.org/dc/elements/1.1/description 'Make sure that the defined survey covers the chosen location!'
start_freq = 72e6  # http://odahub.io/ontology#FrequencyHz; oda:lower_limit 30000; oda:upper_limit 30000000000; http://purl.org/dc/elements/1.1/description 'Consider `survey` frequencies!'
num_of_channels = 1  # http://odahub.io/ontology#Integer; oda:lower_limit 1
freq_inc_hz = (
    8e6  # http://odahub.io/ontology#FrequencyHz; oda:lower_limit 1000
)
num_of_time_steps = 1  # http://odahub.io/ontology#Integer; oda:lower_limit 1
telescope = "ASKAP"  # http://odahub.io/ontology#String; oda:allowed_value 'ASKAP','MeerKAT','LOFAR','MKATPlus','PDBI','SKA1LOW','SKA1MID','VLBA','WSRT'; http://purl.org/dc/elements/1.1/description 'Make sure that the telescope is compatible with the frequencies you observe!'
field_of_view = 1.0  # http://odahub.io/ontology#AngleDegrees; oda:lower_limit 0.1; oda:upper_limit 10.0; http://purl.org/dc/elements/1.1/description 'Consider small enough field-of-view to avoid time & memory issues! Image is 2048x2048.'

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

# some sanity-checks which couldn't be cached using ontology-restrictions

imaging_npixel = 2048  # fixed
start_date_time = pd.to_datetime(
    start_time
).to_pydatetime()  # robust datetime-parser
end_date_time = pd.to_datetime(
    end_time
).to_pydatetime()  # robust datetime-parser
if end_date_time <= start_date_time:
    raise RuntimeError(f"{end_time=} <= {start_time=} is not allowed!")
observation_length = end_date_time - start_date_time
observation_length_hours = observation_length.total_seconds() / 3600
obs_length_hours_min, obs_length_hours_max = 0.5, 8.0
if (
    observation_length_hours < obs_length_hours_min
    or observation_length_hours > obs_length_hours_max
):
    raise RuntimeError(
        "Choose observation-length `end_time` - `start_time` between "
        + f"{obs_length_hours_min} and {obs_length_hours_max} hours, "
        + f"and not {observation_length_hours=}!"
    )

def get_sky_from_survey(
    survey: str,
    ra: float,
    dec: float,
    start_freq: float,
    num_of_channels: int,
    freq_inc_hz: float,
    field_of_view: float,
) -> SkyModel:
    def extract_survey_freqs(survey: str) -> tuple[float, float]:
        extract = lambda idx: float(
            survey.split(sep=" ")[1].split(sep="-")[idx]
        )
        return extract(0), extract(1)

    end_freq = start_freq + freq_inc_hz * num_of_channels
    survey_start_freq, survey_end_freq = extract_survey_freqs(survey=survey)
    if start_freq < survey_start_freq:
        raise RuntimeError(
            f"Chosen {start_freq=} is not in the {survey=} freq range!"
        )
    if end_freq > survey_end_freq:
        raise RuntimeError(
            f"Chosen {start_freq=}, {freq_inc_hz=} & {num_of_channels=} results in "
            + f"{end_freq=} which is not in the {survey=} freq range!"
        )
    if survey == "gleam 72e6-231e6 Hz, RA[0.0,360.0]/DEC[-90.0,30.0] deg":
        survey_sky = SkyModel.get_GLEAM_Sky(
            min_freq=start_freq, max_freq=end_freq
        )
    elif (
        survey
        == "mightee-l1 1304e6-1375e6 Hz, RA[149.40,150.84]/DEC[1.49,2.93] deg"
    ):
        survey_sky = SkyModel.get_MIGHTEE_Sky(
            min_freq=start_freq, max_freq=end_freq
        )
    elif (
        survey
        == "mals-dr1v3 902e6-1644e6 Hz, RA[0.0,360.0]/DEC[-78.80,32.35] deg"
    ):
        survey_sky = SkyModel.get_MALS_DR1V3_Sky(
            min_freq=start_freq, max_freq=end_freq
        )
    else:
        # unreachable if `allowed_value` works with according Karabo version
        raise NotImplementedError(f"{survey=} is not available!")

    outer_filter_radius = (
        field_of_view * np.sqrt(2) / 2
    ) + 2  # calc radius + 2 of outer edge of dirty-image in deg
    survey_sky = survey_sky.filter_by_radius(
        inner_radius_deg=0,
        outer_radius_deg=outer_filter_radius,
        ra0_deg=ra,
        dec0_deg=dec,
    )
    if survey_sky.num_sources <= 0:
        raise RuntimeError(
            f"No sourcs available at {RA=}, {DEC=} for chosen {field_of_view=} of {survey=}. "
            + f"It is also possible, that {start_freq=} & {end_freq=} filtering results "
            + "in having no sources. Please check the survey's frequencies and area-coverage."
        )
    return survey_sky

sky = get_sky_from_survey(
    survey=survey,
    ra=RA,
    dec=DEC,
    start_freq=start_freq,
    num_of_channels=num_of_channels,
    freq_inc_hz=freq_inc_hz,
    field_of_view=field_of_view,
)

if telescope not in get_args(OSKARTelescopesWithoutVersionType):
    # theoretically unreachable if `allowed_value` works as intended with according Karabo version
    raise NotImplementedError(f"Provided {telescope=} is not valid.")
tel = Telescope.constructor(telescope)

observation_settings = Observation(
    start_frequency_hz=start_freq,
    start_date_and_time=start_date_time,
    length=observation_length,
    phase_centre_ra_deg=RA,
    phase_centre_dec_deg=DEC,
    number_of_channels=num_of_channels,
    number_of_time_steps=num_of_time_steps,
)

vis = Visibility()
vis_path = vis.vis_path

interferometer_sim = InterferometerSimulation(
    vis_path=vis_path,
    channel_bandwidth_hz=freq_inc_hz,  # `freq_inc_hz` should be fine
)
visibility = interferometer_sim.run_simulation(
    tel,
    sky,
    observation_settings,
)

# ### Dirty Images
#
# We can create dirty images of visibilites and display them as shown below

imaging_cellsize = np.deg2rad(field_of_view) / imaging_npixel

imager = Imager(
    visibility,
    imaging_npixel=imaging_npixel,
    imaging_cellsize=imaging_cellsize,
)
dirty = imager.get_dirty_image()
if np.sum(np.abs(dirty.data)) == 0:
    err_msg = (
        "No sources in dirty-image. This is probably because the beam is below the horizon. "
        + "You may try again at a different time of day."
    )
    raise RuntimeError(err_msg)
fname_dirty = os.path.join(os.path.split(dirty.path)[0], "dirty.png")
plt.rcParams["savefig.dpi"] = 500
survey_name = survey.split(sep=" ")[0]
dirty.plot(  # saves dirty-image as matplotlib plot (has less pixels than 2048x2048)
    title=f"DIRTY Image {survey_name}, FOV:{field_of_view}",
    filename=fname_dirty,
)
bin_image = PictureProduct.from_file(fname_dirty)

# The next cell is tagged as `outputs`

dirty_png = bin_image  # http://odahub.io/ontology#ODAPictureProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    ("out_karabo_workflow_dirty_png", "dirty_png_galaxy.output", dirty_png)
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
