#!/usr/bin/env python
# coding: utf-8

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

# flake8: noqa

import json
import os

from pipeline_ra_dec import rule_based_ra_dec_detector

text = "The candidate Be/X-ray Binary system RX J0032.9-7348 (ATEL #16880) has been detected in 2 observations conducted by the Swift SMC Survey (S-CUBED) over the past month. This system was first detected by S-CUBED on 22 October 2024 with a 0.3-10 keV band flux of 1.46 x 10^-11 erg/cm^2/s. The most recent detection, obtained by Swift on 11 November 2024, shows that the source has since brightened to a 0.3-10 keV flux of 2.77 x 10^-11 erg/cm^2/s. The time-averaged spectrum of both observations is best fit to an absorbed power law with a photon index of 0.7 (+0.5, -0.4) and a column density of NH = 4.2 (+29.9, - 4.2) x 10^20 cm^-2.Using S-CUBED data, we can refine the localization of this source to a new position of: R.A. = 00:32:59.30 (8.22459d), Dec. = -73:48:28.1 (-73.8078d) with a 90% confidence error region of 4.3 arc-seconds. This new position places RX J0032.9-7348 within 5â of two early-type stars. One of these stars is the B0.5V star GSC 09141-01338 (Evans et al., 2004). The other is marked as Object 1 in Figure 1a of Stevens el al. (1999), where it is observed to have an emission feature that corresponds to the Balmer Series Hydrogen alpha line. However, those authors were searching a large ROSAT error region and warned that there were many such H-alpha emitters in the region. At this time, the exact optical counterpart cannot be determined. Further observations are required to distinguish between the two stars. "
text_id = "16900"

rule_based_ra_dec_detector(text_id, text)

from astropy import units as u
from astropy.coordinates import SkyCoord

SkyCoord(ra="00:32:59.30", dec="-73:48:28.1", unit=(u.hourangle, u.deg))

# output gathering
_galaxy_meta_data = {}

with open(os.path.join(_galaxy_wd, "galaxy.json"), "w") as fd:
    json.dump(_galaxy_meta_data, fd)
print("*** Job finished successfully ***")
