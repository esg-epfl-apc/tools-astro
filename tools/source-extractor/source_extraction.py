#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python

# This script is generated with nb2galaxy

# flake8: noqa

import json
import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
import sep
import tifffile
from astropy.io import fits
from astropy.table import Table
from matplotlib import rcParams
from matplotlib.patches import Ellipse
from oda_api.json import CustomJSONEncoder

get_ipython().run_line_magic("matplotlib", "inline")   # noqa: F821

rcParams["figure.figsize"] = [10.0, 8.0]

input_file = "./input.fits"  # oda:POSIXPath; oda:label "Input file"

### These params are for both functions
mask_file = None  # oda:POSIXPath, oda:optional; oda:label "Mask file"

### These params are for sep.extract()
thresh = 1.5  # oda:Float
err_option = "float_globalrms"  # oda:String; oda:allowed_value 'float_globalrms','array_rms', 'none'
# gain = None
maskthresh = 0.0  # oda:Float
minarea = 5  # oda:Integer
filter_case = "default"  # oda:String; oda:label "Filter Case"; oda:allowed_value 'none', 'default', 'file'
filter_file = None  # oda:POSIXPath, oda:optional; oda:label "Filter file"
filter_type = "matched"  # oda:String; oda:allowed_value 'matched','conv'
deblend_nthresh = 32  # oda:Integer
deblend_cont = 0.005  # oda:Float
clean = True  # oda:Boolean
clean_param = 1.0  # oda:Float

### These params are for sep.Background()
# maskthresh = 0.0
bw = 64  # oda:Integer
bh = 64  # oda:Integer
fw = 3  # oda:Integer
fh = 3  # oda:Integer
fthresh = 0.0  # oda:Float

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "C_data_product_" in inp_dic.keys():
    inp_pdic = inp_dic["C_data_product_"]
else:
    inp_pdic = inp_dic
input_file = str(inp_pdic["input_file"])

mask_file = (
    str(inp_pdic["mask_file"])
    if inp_pdic.get("mask_file", None) is not None
    else None
)

thresh = float(inp_pdic["thresh"])
err_option = str(inp_pdic["err_option"])
maskthresh = float(inp_pdic["maskthresh"])
minarea = int(inp_pdic["minarea"])
filter_case = str(inp_pdic["filter_case"])

filter_file = (
    str(inp_pdic["filter_file"])
    if inp_pdic.get("filter_file", None) is not None
    else None
)

filter_type = str(inp_pdic["filter_type"])
deblend_nthresh = int(inp_pdic["deblend_nthresh"])
deblend_cont = float(inp_pdic["deblend_cont"])
clean = bool(inp_pdic["clean"])
clean_param = float(inp_pdic["clean_param"])
bw = int(inp_pdic["bw"])
bh = int(inp_pdic["bh"])
fw = int(inp_pdic["fw"])
fh = int(inp_pdic["fh"])
fthresh = float(inp_pdic["fthresh"])

try:
    hdul = fits.open(input_file)
    data = hdul[0].data
    data = data.astype(data.dtype.newbyteorder("=")).astype(float)
except:
    try:
        data = tifffile.imread(input_file).astype(float)
    except:
        raise RuntimeError(
            "The input file should have the FITS or TIFF format."
        )

print("INFO: Data shape:", data.shape)

if mask_file is not None:
    try:
        hdul = fits.open(mask_file)
        mask = hdul[0].data
        mask = mask.astype(mask.dtype.newbyteorder("="))
    except:
        try:
            mask = tifffile.imread(mask_file)
        except:
            raise RuntimeError(
                "The mask file should have the FITS or TIFF format."
            )
else:
    mask = None

print("INFO: Mask type:", type(mask))

filter_kernel = None
if filter_case == "none":
    filter_kernel = None
elif filter_case == "default":
    filter_kernel = np.array([[1, 2, 1], [2, 4, 2], [1, 2, 1]])
elif filter_case == "file":
    try:
        filter_kernel = np.loadtxt(filter_file)
    except:
        raise RuntimeError(
            "The filter file should be a text file that is loaded with numpy.loadtxt"
        )

print("INFO: Filter kernel:", filter_kernel)

# measure a spatially varying background on the image
bkg = sep.Background(
    data,
    mask=mask,
    maskthresh=maskthresh,
    bw=bw,
    bh=bh,
    fw=fw,
    fh=fh,
    fthresh=fthresh,
)

# evaluate background as 2-d array, same size as original image
bkg_array = bkg.back()
hdu_bkg = fits.PrimaryHDU(bkg_array)
hdu_bkg.writeto("bkg_array.fits")

# evaluate the background noise as 2-d array, same size as original image
bkg_rms = bkg.rms()
hdu_rms = fits.PrimaryHDU(bkg_rms)
hdu_rms.writeto("bkg_rms.fits")

# subtract the background
data_sub = data - bkg

if err_option == "float_globalrms":
    err = bkg.globalrms
elif err_option == "array_rms":
    err = bkg_rms
else:
    err = None

# extract sources:
objects, segmap = sep.extract(
    data_sub,
    thresh,
    err=err,
    gain=None,
    mask=mask,
    maskthresh=maskthresh,
    minarea=minarea,
    filter_kernel=filter_kernel,
    filter_type=filter_type,
    deblend_nthresh=deblend_nthresh,
    deblend_cont=deblend_cont,
    clean=clean,
    clean_param=clean_param,
    segmentation_map=True,
)

# show the background
fig, ax = plt.subplots()
im = ax.imshow(bkg_array, interpolation="nearest", cmap="gray", origin="lower")
fig.colorbar(im)
fig.savefig("./bkg_image.png", format="png", bbox_inches="tight")

# show the background noise
fig, ax = plt.subplots()
im = ax.imshow(bkg_rms, interpolation="nearest", cmap="gray", origin="lower")
ax.set_title(f"This is array_rms. While float_globalrms={bkg.globalrms}")
fig.colorbar(im)
fig.savefig("./bkg_rms.png", format="png", bbox_inches="tight")

# plot image
fig, ax = plt.subplots()
m, s = np.mean(data), np.std(data)
im = ax.imshow(
    data,
    interpolation="nearest",
    cmap="gray",
    vmin=m - s,
    vmax=m + s,
    origin="lower",
)
fig.colorbar(im)
fig.savefig("./fits2image.png", format="png", bbox_inches="tight")

# show the segmentation map
fig, ax = plt.subplots()
im = ax.imshow(
    segmap,
    interpolation="nearest",
    cmap="gray",
    origin="lower",
    vmin=0,
    vmax=1,
)
fig.colorbar(im)
fig.savefig("./segmap.png", format="png", bbox_inches="tight")

# plot background-subtracted image
fig, ax = plt.subplots()
m, s = np.mean(data_sub), np.std(data_sub)
im = ax.imshow(
    data_sub,
    interpolation="nearest",
    cmap="gray",
    vmin=m - s,
    vmax=m + s,
    origin="lower",
)

# plot an ellipse for each object
for i in range(len(objects)):
    e = Ellipse(
        xy=(objects["x"][i], objects["y"][i]),
        width=6 * objects["a"][i],
        height=6 * objects["b"][i],
        angle=objects["theta"][i] * 180.0 / np.pi,
    )
    e.set_facecolor("none")
    e.set_edgecolor("red")
    ax.add_artist(e)

fig.savefig("./sources.png", format="png", bbox_inches="tight")

plt.show()

from oda_api.data_products import ODAAstropyTable

cat = ODAAstropyTable(Table(data=objects))
hdu_rms = fits.PrimaryHDU(segmap.astype("uint32"))
hdu_rms.writeto("segmentation_map.fits")

bkg_picture = "./bkg_image.png"  # oda:POSIXPath
rms_picture = "./bkg_rms.png"  # oda:POSIXPath
data_picture = "./fits2image.png"  # oda:POSIXPath
sources_picture = "./sources.png"  # oda:POSIXPath
segmentation_map_picture = "./segmap.png"  # oda:POSIXPath
segmentation_map = "./segmentation_map.fits"  # oda:POSIXPath
bkg_array = "./bkg_array.fits"  # oda:POSIXPath
rms_array = "./bkg_rms.fits"  # oda:POSIXPath
catalog_table = cat  # oda:ODAAstropyTable

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_source_extraction_catalog_table",
        "catalog_table_galaxy.output",
        catalog_table,
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
_simple_outs = []
_simple_outs.append(
    (
        "out_source_extraction_bkg_picture",
        "bkg_picture_galaxy.output",
        bkg_picture,
    )
)
_simple_outs.append(
    (
        "out_source_extraction_rms_picture",
        "rms_picture_galaxy.output",
        rms_picture,
    )
)
_simple_outs.append(
    (
        "out_source_extraction_data_picture",
        "data_picture_galaxy.output",
        data_picture,
    )
)
_simple_outs.append(
    (
        "out_source_extraction_sources_picture",
        "sources_picture_galaxy.output",
        sources_picture,
    )
)
_simple_outs.append(
    (
        "out_source_extraction_segmentation_map_picture",
        "segmentation_map_picture_galaxy.output",
        segmentation_map_picture,
    )
)
_simple_outs.append(
    (
        "out_source_extraction_segmentation_map",
        "segmentation_map_galaxy.output",
        segmentation_map,
    )
)
_simple_outs.append(
    ("out_source_extraction_bkg_array", "bkg_array_galaxy.output", bkg_array)
)
_simple_outs.append(
    ("out_source_extraction_rms_array", "rms_array_galaxy.output", rms_array)
)
_numpy_available = True

for _outn, _outfn, _outv in _simple_outs:
    _galaxy_outfile_name = os.path.join(_galaxy_wd, _outfn)
    if isinstance(_outv, str) and os.path.isfile(_outv):
        shutil.move(_outv, _galaxy_outfile_name)
        _galaxy_meta_data[_outn] = {"ext": "_sniff_"}
    elif _numpy_available and isinstance(_outv, np.ndarray):
        with open(_galaxy_outfile_name, "wb") as fd:
            np.savez(fd, _outv)
        _galaxy_meta_data[_outn] = {"ext": "npz"}
    else:
        with open(_galaxy_outfile_name, "w") as fd:
            json.dump(_outv, fd)
        _galaxy_meta_data[_outn] = {"ext": "expression.json"}

with open(os.path.join(_galaxy_wd, "galaxy.json"), "w") as fd:
    json.dump(_galaxy_meta_data, fd)
print("*** Job finished successfully ***")
