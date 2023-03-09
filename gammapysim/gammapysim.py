import matplotlib.pyplot as plt
from pathlib import Path
import click

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
from gammapy.makers import MapDatasetMaker
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

def save_events(events, dataset, number):
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
    hdu_all.writeto(f"./events_{number}.fits", overwrite=True)



from gammapy.modeling.models import PowerLawSpectralModel
from gammapy.modeling.models import GaussianSpatialModel

def prepare():
    filename = "data/Prod5-North-20deg-AverageAz-4LSTs09MSTs.180000s-v0.1.fits.gz"
    IRFS = load_cta_irfs(filename)


    pwl = PowerLawSpectralModel(amplitude="2.7e-12 TeV-1 cm-2 s-1", index=2.2)
    gauss = GaussianSpatialModel(
        lon_0="0 deg", lat_0="0 deg", sigma="0.2 deg", frame="galactic"
    )

    modelsky = [SkyModel(spectral_model=pwl, spatial_model=gauss, name="my-source")]
    print(modelsky)

    bkg_model = FoVBackgroundModel(dataset_name="my-dataset")
    modelsky.append(bkg_model)

    return modelsky, IRFS
    
   
def synth_for_pointing(i, pointing):
    modelsky, IRFS = prepare()

    print(f"Make the observation for pointing {i}...")
    observation = Observation.create(
                      obs_id="{:06d}".format(i), pointing=pointing, 
                      livetime="0.01 hr", 
                      irfs=IRFS
                      )

    print(f"Create the dataset for {pointing}")
    energy_axis = MapAxis.from_energy_bounds(
        "0.012 TeV", "100 TeV", nbin=2, per_decade=True
        )
    energy_axis_true = MapAxis.from_energy_bounds(
        "0.001 TeV", "300 TeV", nbin=5, per_decade=True, name="energy_true"
        )
    migra_axis = MapAxis.from_bounds(
        0.5, 2, nbin=15, node_type="edges", name="migra"
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
    maker = MapDatasetMaker(selection=["exposure", "background", "psf", "edisp"])
    dataset = maker.run(empty, observation)

    rad_size = 1
    region_sky = CircleSkyRegion(center=pointing, radius=rad_size * u.deg)
    mask_map = dataset.geoms["geom"].region_mask(region_sky)
    # mod = modelsky.select_mask(mask_map)
    mod = modelsky    

    # bkg_idx = np.where(np.array(modelsky.names) == 'my-dataset-bkg')
    # mod.append(modelsky[int(bkg_idx[0][0])])

    dataset.models = mod

    for m in dataset.models[:-1]:
        sep = m.spatial_model.position.separation(pointing).deg
        print(f"This is the spatial separation of {m.name} from the pointing direction: {sep}")

    print("Simulate...")
    sampler = MapDatasetEventSampler(random_state=i)
    events = sampler.run(dataset, observation)

    print(f"Save events {i}...")
    save_events(events, dataset, "{:06d}".format(i))

@click.command()
@click.option("--ra", default=0, type=float)
@click.option("--dec", default=0, type=float)
def main(ra, dec):
    click.echo(f"generating pointing for RA, Dec = {ra}, {dec}")
    pointing = SkyCoord(ra, dec, unit="deg", frame="galactic")

    synth_for_pointing(0, pointing)


if __name__ == "__main__":
    main()