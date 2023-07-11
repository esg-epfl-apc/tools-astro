from pathlib import Path
import click

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from astropy.io import fits

# %matplotlib inline
from IPython.display import display
from gammapy.data import Observation, observatory_locations
from gammapy.datasets import MapDataset, MapDatasetEventSampler
from gammapy.irf import load_irf_dict_from_file
from gammapy.makers import MapDatasetMaker, SafeMaskMaker
from gammapy.maps import MapAxis, WcsGeom, Map
from gammapy.modeling import Fit
from gammapy.modeling.models import (
    FoVBackgroundModel,
    GaussianSpatialModel,
    Models,
    PowerLawSpectralModel,
    SkyModel,
)

   
def synth_for_pointing(livetime, pointing, output_events, output_peek_png, output_peek_fits):
    # Loading IRFs
    irfs = load_irf_dict_from_file(
        "$GAMMAPY_DATA/cta-1dc/caldb/data/cta/1dc/bcf/South_z20_50h/irf_file.fits"
    )

    # Define map geometry for binned simulation
    energy_reco = MapAxis.from_edges(
        np.logspace(-1.0, 1.0, 10), unit="TeV", name="energy", interp="log"
    )
    geom = WcsGeom.create(
        skydir=(0, 0),
        binsz=0.02,
        width=(6, 6),
        frame="galactic",
        axes=[energy_reco],
    )
    # It is usually useful to have a separate binning for the true energy axis
    energy_true = MapAxis.from_edges(
        np.logspace(-1.5, 1.5, 30), unit="TeV", name="energy_true", interp="log"
    )

    migra_axis = MapAxis.from_bounds(0.5, 2, nbin=150, node_type="edges", name="migra")

    empty = MapDataset.create(geom, name="dataset-simu", energy_axis_true=energy_true, migra_axis=migra_axis)

    # Define sky model to used simulate the data.
    # Here we use a Gaussian spatial model and a Power Law spectral model.
    spatial_model = GaussianSpatialModel(
        lon_0="0.2 deg", lat_0="0.1 deg", sigma="0.3 deg", frame="galactic"
    )
    spectral_model = PowerLawSpectralModel(
        index=3, amplitude="1e-11 cm-2 s-1 TeV-1", reference="1 TeV"
    )
    model_simu = SkyModel(
        spatial_model=spatial_model,
        spectral_model=spectral_model,
        name="model-simu",
    )

    bkg_model = FoVBackgroundModel(dataset_name="dataset-simu")

    models = Models([model_simu, bkg_model])
    print(models)

    file_model = "point-pwl.yaml"
    models.write(file_model, overwrite=True)

    # Create an in-memory observation
    location = observatory_locations["cta_south"]
    obs = Observation.create(
        pointing=pointing, livetime=livetime, irfs=irfs, location=location
    )
    print(obs)

    # Make the MapDataset
    maker = MapDatasetMaker(selection=["exposure", "background", "psf", "edisp"])

    maker_safe_mask = SafeMaskMaker(methods=["offset-max"], offset_max=4.0 * u.deg)

    dataset = maker.run(empty, obs)
    dataset = maker_safe_mask.run(dataset, obs)
    print(dataset)

    dataset.write(output_events, overwrite=True)

    # Add the model on the dataset and Poisson fluctuate
    dataset.models = models
    dataset.fake()
    # Do a print on the dataset - there is now a counts maps
    print(dataset)

    sampler = MapDatasetEventSampler(random_state=0)
    events = sampler.run(dataset, obs)

    print(f"Source events: {(events.table['MC_ID'] == 1).sum()}")
    print(f"Background events: {(events.table['MC_ID'] == 0).sum()}")

    events.peek()
    plt.savefig(output_peek_png)
    events.select_offset([0, 1] * u.deg).peek()
    # plt.savefig("peek-focus.png")

    counts = Map.from_geom(geom)
    counts.fill_events(events)
    counts.sum_over_axes().plot(add_cbar=True)

    counts.write(output_peek_fits)

    # To plot, eg, counts:
    # dataset.counts.smooth(0.05 * u.deg).plot_interactive(add_cbar=True, stretch="linear")
    # plt.savefig("counts.png")


@click.command()
@click.option("--obs_id", default=1, help="Observation ID")
@click.option("--livetime-hr", default=1, type=float, help="Livetime (hours)")
@click.option("--pointing-coord", default="0 0", help="Pointing coordinate (SkyCoord)")
@click.option("--output-events", default="events.fits", help="Output events file")
@click.option("--output-peek-png", default="peek.png", help="Output peek file")
@click.option("--output-peek-fits", default="peek.fits")
def main(obs_id, livetime_hr, pointing_coord, output_events, output_peek_png, output_peek_fits):
    livetime = livetime_hr * u.hr
    pointing = SkyCoord(pointing_coord, unit="deg", frame="galactic")

    synth_for_pointing(livetime, pointing, output_events, output_peek_png, output_peek_fits)

if __name__ == "__main__":
    main()