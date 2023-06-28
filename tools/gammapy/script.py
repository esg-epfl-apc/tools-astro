import os
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.coordinates import SkyCoord, Angle
from regions import CircleSkyRegion
from gammapy.data import DataStore, Observation
from gammapy.datasets import MapDataset, MapDatasetEventSampler
from gammapy.estimators import LightCurveEstimator
from gammapy.maps import MapAxis, WcsGeom, Map, RegionGeom
from gammapy.irf import load_cta_irfs
from gammapy.makers import MapDatasetMaker, SafeMaskMaker
from gammapy.modeling import Fit
from gammapy.modeling.models import (
    Model,
    Models,
    SkyModel,
    PowerLawSpectralModel,
    PowerLawNormSpectralModel,
    TemplateSpectralModel,
    PointSpatialModel,
    GaussianSpatialModel,
    TemplateSpatialModel,
    ConstantTemporalModel,
    ExpDecayTemporalModel,
    LightCurveTemplateTemporalModel,
    FoVBackgroundModel,
)
from regions import CircleSkyRegion, PointSkyRegion

def save_events(events, dataset, fn):
    """Save simulated event list.

    Input

    events: ~gammapy.EventList
    dataset: ~gammapy.Dataset
    number: identifier
    """
    primary_hdu = fits.PrimaryHDU()
    hdu_evt = fits.BinTableHDU(events.table)
    hdu_gti = fits.BinTableHDU(dataset.gti.table, name="GTI")
    hdu_all = fits.HDUList([primary_hdu, hdu_evt, hdu_gti])
    hdu_all.writeto(fn, overwrite=True)

filename = os.path.join(os.getenv("GAMMAPY_DATA"), "cta-caldb/Prod5-North-20deg-AverageAz-4LSTs09MSTs.180000s-v0.1.fits.gz")
IRFS = load_cta_irfs(filename)

# Define sky model to used simulate the data.
# Here we use a Gaussian spatial model and a Power Law spectral model.
spatial_model = GaussianSpatialModel(
    lon_0="0.2 deg", lat_0="0.1 deg", sigma="0.3 deg", frame="galactic"
)
spectral_model = PowerLawSpectralModel(
    index=3, amplitude="1e-10 cm-2 s-1 TeV-1", reference="1 TeV"
)
model_simu = SkyModel(
    spatial_model=spatial_model,
    spectral_model=spectral_model,
    name="model-simu",
)

bkg_model = FoVBackgroundModel(dataset_name="dataset-simu")

models = Models([model_simu, bkg_model])
print(models)

   
def synth_for_pointing(i, pointing, livetime):
    print(f"Make the observation for pointing {i}...")
    obs = Observation.create(
                      obs_id="{:06d}".format(i), pointing=pointing, 
                      livetime=livetime, 
                      irfs=IRFS
                      )

    print(f"Create the dataset for {pointing}")
    energy_axis = MapAxis.from_energy_bounds(
        "0.012 TeV", "100 TeV", nbin=10, per_decade=True
        )
    energy_axis_true = MapAxis.from_energy_bounds(
        "0.001 TeV", "300 TeV", nbin=20, per_decade=True, name="energy_true"
        )
    migra_axis = MapAxis.from_bounds(
        0.5, 2, nbin=50, node_type="edges", name="migra"
        )

    geom = WcsGeom.create(
        skydir=pointing,
        width=(12, 12),
        binsz=0.02,
        frame="icrs",
        axes=[energy_axis],
    )

    empty = MapDataset.create(
            geom,
            energy_axis_true=energy_axis_true,
            migra_axis=migra_axis,
            name="my-dataset",
                )

    # Make the MapDataset
    maker = MapDatasetMaker(selection=["exposure", "background", "psf", "edisp"])

    maker_safe_mask = SafeMaskMaker(methods=["offset-max"], offset_max=4.0 * u.deg)

    dataset = maker.run(empty, obs)
    dataset = maker_safe_mask.run(dataset, obs)
    print(dataset)

    print("Simulate...")

    # Add the model on the dataset and Poission fluctuate
    dataset.models = models
    dataset.fake()
    # Do a print on the dataset - there is now a counts maps
    print(dataset)

    # To plot, eg, counts:
    # dataset.counts.smooth(0.05 * u.deg).plot(
    #     add_cbar=True, stretch="linear"
    # )

    sampler = MapDatasetEventSampler(random_state=i)
    events = sampler.run(dataset, obs)

    print(f"Save events {i}...")
    save_events(events, dataset, "$events")


livetime = 1.0 * u.hr
pointing = SkyCoord(0, 0, unit="deg", frame="galactic")

synth_for_pointing(0, pointing, livetime)