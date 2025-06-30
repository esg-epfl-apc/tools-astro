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
from astropy.io import fits
from astropy.table import Table
from matplotlib import rcParams
from matplotlib.patches import Ellipse
from oda_api.json import CustomJSONEncoder

get_ipython().run_line_magic("matplotlib", "inline")   # noqa: F821

rcParams["figure.figsize"] = [10.0, 8.0]

input_file = "./input.fits"  # oda:POSIXPath

### sep.extract()
thresh = 1.5  # oda:Float
err_option = "float_globalrms"  # oda:String; oda:allowed_value 'float_globalrms','array_rms'
# gain = None
# mask = None
minarea = 5  # oda:Integer
# filter_kernel=default_kernel
filter_type = "matched"  # oda:String; oda:allowed_value 'matched','conv'
deblend_nthresh = 32  # oda:Integer
deblend_cont = 0.005  # oda:Float
clean = True  # oda:Boolean ; oda:allowed_value True,False
clean_param = 1.0  # oda:Float
# segmentation_map=False

### sep.Background()
# mask = None
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
thresh = float(inp_pdic["thresh"])
err_option = str(inp_pdic["err_option"])
minarea = int(inp_pdic["minarea"])
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

hdul = fits.open(input_file)
data = hdul[0].data
data = data.astype(data.dtype.newbyteorder("="))

# measure a spatially varying background on the image
bkg = sep.Background(
    data,
    mask=None,
    maskthresh=0.0,
    bw=bw,
    bh=bh,
    fw=fw,
    fh=fh,
    fthresh=fthresh,
)

# evaluate background as 2-d array, same size as original image
bkg_image = bkg.back()

# evaluate the background noise as 2-d array, same size as original image
bkg_rms = bkg.rms()

# subtract the background
data_sub = data - bkg

if err_option == "float_globalrms":
    err = bkg.globalrms
else:
    err = bkg.rms()

# extract sources:
objects = sep.extract(
    data_sub,
    thresh,
    err=err,
    gain=None,
    mask=None,
    minarea=minarea,
    filter_type=filter_type,
    deblend_nthresh=deblend_nthresh,
    deblend_cont=deblend_cont,
    clean=clean,
    clean_param=clean_param,
)

# show the background
fig, ax = plt.subplots()
im = ax.imshow(bkg_image, interpolation="nearest", cmap="gray", origin="lower")
fig.colorbar(im)
fig.savefig("bkg_image.png", format="png", bbox_inches="tight")

# show the background noise
fig, ax = plt.subplots()
im = ax.imshow(bkg_rms, interpolation="nearest", cmap="gray", origin="lower")
fig.colorbar(im)
fig.savefig("bkg_rms.png", format="png", bbox_inches="tight")

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
fig.savefig("fits2image.png", format="png", bbox_inches="tight")

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

fig.savefig("sources.png", format="png", bbox_inches="tight")

plt.show()

from oda_api.data_products import ODAAstropyTable, PictureProduct

bin_bkg_image = PictureProduct.from_file("bkg_image.png")
bin_rms_image = PictureProduct.from_file("bkg_rms.png")
bin_dat_image = PictureProduct.from_file("fits2image.png")
bin_sor_image = PictureProduct.from_file("sources.png")

cat = ODAAstropyTable(Table(data=objects))

picture1 = bin_bkg_image  # oda:ODAPictureProduct
picture2 = bin_rms_image  # oda:ODAPictureProduct
picture3 = bin_dat_image  # oda:ODAPictureProduct
picture4 = bin_sor_image  # oda:ODAPictureProduct
catalog_table = cat  # oda:ODAAstropyTable

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    ("out_source_extraction_picture1", "picture1_galaxy.output", picture1)
)
_oda_outs.append(
    ("out_source_extraction_picture2", "picture2_galaxy.output", picture2)
)
_oda_outs.append(
    ("out_source_extraction_picture3", "picture3_galaxy.output", picture3)
)
_oda_outs.append(
    ("out_source_extraction_picture4", "picture4_galaxy.output", picture4)
)
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

with open(os.path.join(_galaxy_wd, "galaxy.json"), "w") as fd:
    json.dump(_galaxy_meta_data, fd)
print("*** Job finished successfully ***")
