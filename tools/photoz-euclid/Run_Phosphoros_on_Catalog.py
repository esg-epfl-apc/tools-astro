#!/usr/bin/env python
# coding: utf-8

# flake8: noqa

import json
import os
import shutil
from os.path import join

import astropy.table
import numpy as np
from astropy.coordinates import FK5, SkyCoord
from astropy.io import ascii, fits
from oda_api.json import CustomJSONEncoder

# # Run Phosphoros
#
# Use pre-computed Model Grids (MGs) and Galactic Correction Coefficient Grids (GCCGs)
#
# - read input catalog from url or upload file
#
# - read pre-computed MGs and GCCGs for the instruments selected
#
# - build new MGs and GCCGs
#
# - map filters to catalog columns
#
# - run Phosphoros

get_ipython().run_line_magic("matplotlib", "inline")   # noqa: F821
import subprocess

import matplotlib.pyplot as plt
from oda_api.data_products import (
    BinaryProduct,
    ODAAstropyTable,
    PictureProduct,
)

"""
List of filters/fluxes available in pre-computed ModelGrids

CFHT|MegaCam i,g,gri,r,u,z
DECam|DECam Y,g,i,r,u,z
Euclid|NISP Y,J,H  |VIS vis
GAIA|GAIA3 G,Gbp,Grp
GALEX|GALEX FUV,NUV
LSST|LSST u,g,r,i,z,y
PANSTARRS|PS1 g,r,i,z,y,w
VISTA|VISTA Z,Y,J,H,Ks
SDSS|SDSS u,g,r,i,z
HSC|HSC g,r,i,z,Y
UKIRT|WFCAM Z,Y,J,H,K
WISE|WISE W1,W2,W3,W4

Names to be corrected (from name in Frontend to name in the code)
CTIO --> DECam
PAN-STARRS --> PS1
VIRCAM --> VISTA
SLOAN --> SDSS
Subaru --> HSC

New Instruments + new filters

JWST|MIRI.* NIRCam.*
HST|ACS_WFC.* WFC3_IR.*
Spitzer|IRAC.I1,IRAC.I2
DESI|bass.g,bass.r,MzLS.z
HSCinter|HSC.*_filter
SuprimeCam

"""

print("Start")

# inputs & parameters

# input catalog
catalog_URL = "https://www.astro.unige.ch/~tucci/Phosphoros/MultiBands_Catalog_1k.fits"  # http://odahub.io/ontology#FileReference ; oda:group "Catalog filter" ; oda:label ""
# for tests
# catalog_URL = 'data/Catalogs/Catalog_Galaxy_n1K.fits' #http://odahub.io/ontology#POSIXPath ; oda:group "Catalog filter" ; oda:label ""
# catalog_URL = 'data/Catalogs/Star_Cat_COSMOS.fits' #http://odahub.io/ontology#FileURL ; oda:group "Catalog filter" ; oda:label ""
# catalog_URL = 'MultiBands_Catalog_1k.fits' #http://odahub.io/ontology#FileReference ; oda:group "Catalog filter" ; oda:label ""
# catalog_URL = 'data/Catalogs/SpecZ_PhotoZ_South.fits'

# oda:oda_token_access oda:InOdaContext .

# filters - flux/error_flux column names
filters_table = {
    "filter": [
        "Euclid|VIS.vis",
        "Euclid|NISP.Y",
        "Euclid|NISP.J",
        "Euclid|NISP.H",
        "LSST|LSST.g",
        "LSST|LSST.r",
        "LSST|LSST.i",
        "LSST|LSST.z",
    ],
    "flux": [
        "FLUX_DETECTION_TOTAL",
        "FLUX_Y_TOTAL",
        "FLUX_J_TOTAL",
        "FLUX_H_TOTAL",
        "FLUX_G_EXT_LSST_TOTAL",
        "FLUX_R_EXT_LSST_TOTAL",
        "FLUX_I_EXT_LSST_TOTAL",
        "FLUX_Z_EXT_LSST_TOTAL",
    ],
    "flux_error": [
        "FLUXERR_DETECTION_TOTAL",
        "FLUXERR_Y_TOTAL",
        "FLUXERR_J_TOTAL",
        "FLUXERR_H_TOTAL",
        "FLUXERR_G_EXT_LSST_TOTAL",
        "FLUXERR_R_EXT_LSST_TOTAL",
        "FLUXERR_I_EXT_LSST_TOTAL",
        "FLUXERR_Z_EXT_LSST_TOTAL",
    ],
}  # http://odahub.io/ontology#PhosphorosFiltersTable ; oda:group "Catalog filter"
# filters_table = {'filter': ['Euclid|VIS.vis','Euclid|NISP.Y','Euclid|NISP.J','Euclid|NISP.H','LSST|LSST.g','LSST|LSST.r','LSST|LSST.i','LSST|LSST.z'], 'flux': ['VIS_vis','NISP_Y','NISP_J','NISP_H','LSST_g','LSST_r','LSST_i','LSST_z'], 'flux_error': ['eVIS_vis','eNISP_Y','eNISP_J','eNISP_H','eLSST_g','eLSST_r','eLSST_i','eLSST_z']} # http://odahub.io/ontology#PhosphorosFiltersTable ; oda:group "Catalog filter"
# for tests
# filters_table = {'filter': ['SLOAN|SDSS.g','Subaru|HSC.r','CTIO|DECam.i','PAN-STARRS|PS1.z','VIRCAM|VISTA.J','Euclid|VIS.vis','Euclid|NISP.Y',], 'flux': ['SDSS_g','HSC_r','DECam_i','PS1_z','VISTA_J','VIS_vis','NISP_Y'], 'flux_error': ['eSDSS_g','eHSC_r','eDECam_i','ePS1_z','eVISTA_J','eVIS_vis','eNISP_Y']} # http://odahub.io/ontology#PhosphorosFiltersTable ; oda:group "Catalog filter"
# filters_table = {'filter': ['SDSS|SDSS.g','HSC|HSC.r','DECam|DECam.i','PANSTARRS|PS1.z','VISTA|VISTA.J','Euclid|VIS.vis','Euclid|NISP.Y',], 'flux': ['SDSS_g','HSC_r','DECam_i','PS1_z','VISTA_J','VIS_vis','NISP_Y'], 'flux_error': ['eSDSS_g','eHSC_r','eDECam_i','ePS1_z','eVISTA_J','eVIS_vis','eNISP_Y']} # http://odahub.io/ontology#PhosphorosFiltersTable ; oda:group "Catalog filter"
# star catalog example
# filters_table = {'filter': ['CTIO|DECam.g','CTIO|DECam.r','CTIO|DECam.i','CTIO|DECam.z','VISTA|VISTA.Y','VISTA|VISTA.J','VISTA|VISTA.H'], 'flux': ['AsinhMag_g','AsinhMag_r','AsinhMag_i','AsinhMag_z','AsinhMag_Y','AsinhMag_J','AsinhMag_H'], 'flux_error': ['AsinhMagErr_g','AsinhMagErr_r','AsinhMagErr_i','AsinhMagErr_z','AsinhMagErr_Y','AsinhMagErr_J','AsinhMagErr_H']} # http://odahub.io/ontology#PhosphorosFiltersTable ; oda:group "Catalog filter"
# filters_table = {'filter': ['CTIO|DECam.g','CTIO|DECam.r','CTIO|DECam.i','CTIO|DECam.z','VISTA|VISTA.Y','VISTA|VISTA.J','VISTA|VISTA.H'],
#                 'flux': ['bdf_flux_dered_calib_g','bdf_flux_dered_calib_r','bdf_flux_dered_calib_i','bdf_flux_dered_calib_z','bdf_flux_dered_calib_Y','bdf_flux_dered_calib_J','bdf_flux_dered_calib_H'],
#                 'flux_error': ['bdf_flux_err_dered_calib_g','bdf_flux_err_dered_calib_r','bdf_flux_err_dered_calib_i','bdf_flux_err_dered_calib_z','bdf_flux_err_dered_calib_Y','bdf_flux_err_dered_calib_J','bdf_flux_err_dered_calib_H']} # http://odahub.io/ontology#PhosphorosFiltersTable ; oda:group "Catalog filter"
#
#
# SPECZ
# filters_table = {'filter': ['DECam|DECam.g','DECam|DECam.r','DECam|DECam.i','DECam|DECam.z','Euclid|VIS.vis','Euclid|NISP.Y','Euclid|NISP.J','Euclid|NISP.H',],
#                 'flux': ['flux_g_ext_decam_unif', 'flux_r_ext_decam_unif', 'flux_i_ext_decam_unif', 'flux_z_ext_decam_unif', 'flux_vis_unif', 'flux_y_unif', 'flux_j_unif', 'flux_h_unif'],
#                 'flux_error': ['fluxerr_g_ext_decam_unif', 'fluxerr_r_ext_decam_unif', 'fluxerr_i_ext_decam_unif', 'fluxerr_z_ext_decam_unif', 'fluxerr_vis_unif', 'fluxerr_y_unif', 'fluxerr_j_unif', 'fluxerr_h_unif']} # http://odahub.io/ontology#PhosphorosFiltersTable ; oda:group "Catalog filter"

# TESTS NEW FILTERS
# filters_table = {'filter': ['JWST|MIRI.F1800W','JWST|NIRCam.F277W','HST|ACS_WFC.F555W','HST|ACS_WFC.G800L','HST|WFC3_IR.F140W','Euclid|VIS.vis'],
#                 'flux': ['FLUX_J_TOTAL','FLUX_H_TOTAL','FLUX_G_EXT_LSST_TOTAL','FLUX_R_EXT_LSST_TOTAL','FLUX_Z_EXT_LSST_TOTAL','FLUX_DETECTION_TOTAL'],
#                 'flux_error': ['FLUXERR_J_TOTAL','FLUXERR_H_TOTAL','FLUXERR_G_EXT_LSST_TOTAL','FLUXERR_R_EXT_LSST_TOTAL','FLUXERR_Z_EXT_LSST_TOTAL','FLUXERR_DETECTION_TOTAL']}
# filters_table = {'filter': ['Spitzer|IRAC.I1','Spitzer|IRAC.I2','DESI|bass.g','DESI|bass.r','DESI|MzLS.z','Euclid|VIS.vis'],
#                 'flux': ['FLUX_J_TOTAL','FLUX_H_TOTAL','FLUX_G_EXT_LSST_TOTAL','FLUX_R_EXT_LSST_TOTAL','FLUX_Z_EXT_LSST_TOTAL','FLUX_DETECTION_TOTAL'],
#                 'flux_error': ['FLUXERR_J_TOTAL','FLUXERR_H_TOTAL','FLUXERR_G_EXT_LSST_TOTAL','FLUXERR_R_EXT_LSST_TOTAL','FLUXERR_Z_EXT_LSST_TOTAL','FLUXERR_DETECTION_TOTAL']}
# filters_table = {'filter': ['VISTA|VISTA.NB118','SuprimeCam|Suprime.B','SuprimeCam|Suprime.r','SuprimeCam|Suprime.IB709','HSCinter|HSC.NB387_filter','HSCinter|HSC.NB973_filter'],
#                 'flux': ['FLUX_J_TOTAL','FLUX_H_TOTAL','FLUX_G_EXT_LSST_TOTAL','FLUX_R_EXT_LSST_TOTAL','FLUX_Z_EXT_LSST_TOTAL','FLUX_DETECTION_TOTAL'],
#                 'flux_error': ['FLUXERR_J_TOTAL','FLUXERR_H_TOTAL','FLUXERR_G_EXT_LSST_TOTAL','FLUXERR_R_EXT_LSST_TOTAL','FLUXERR_Z_EXT_LSST_TOTAL','FLUXERR_DETECTION_TOTAL']}

# flux or magnitude?
ab_magnitude = "Flux"  # http://odahub.io/ontology#String ; oda:allowed_value "Flux", "Magnitude" ; oda:label "Flux [microJy] or AB Magnitude"

# Object Type (Galaxy/Star/QSO)
object_type = "Galaxy"  # http://odahub.io/ontology#String ; oda:allowed_value "Galaxy", "HighZ", "Star", "QSO" ; oda:label "Object Type" ; oda:description "Object type"

# first and last line
skip_sources = 0  # http://odahub.io/ontology#Integer ; oda:label "Skip the first N sources"
proc_sources = 0  # http://odahub.io/ontology#Integer ; oda:label "N sources to process (all if = 0)"

# more input parameters (OPTIONAL)
# MW extinction correction (e.g., 'MW_EBV', 'RA', 'DEC')
column_name_MW_EBV = ""  # http://odahub.io/ontology#String ; oda:label "E(B-V)" ; oda:group "Milky Way Extinction" ; oda:description "E(B-V)"
column_name_RA = ""  # http://odahub.io/ontology#String ; oda:label "Right Ascension" ; oda:group "Milky Way Extinction" ; oda:description "Right Ascension"
column_name_DEC = ""  # http://odahub.io/ontology#String ; oda:label "Declination" ; oda:group "Milky Way Extinction" ; oda:description "Declination"

# include Priors (only for object_type=Galaxy)
priors = "None"  # http://odahub.io/ontology#String ; oda:allowed_value "None", "Volume", "Redshift (Benitez 2000)", "Top-Hat LF" ; oda:label "Priors (only for Galaxies)"
# must be equal to one element of filters_table['flux'], e.g. 'LSST_i'
column_name_Nz_prior_I = (
    ""  # http://odahub.io/ontology#String ; oda:label "Redshift Prior: I Band"
)

# Column with 'spectroscopic" redshifts (for comparison; e.g. 'Ztrue')
column_name_Ztrue = ""  # http://odahub.io/ontology#String ; oda:label "Spectroscopic (Reference) Redshift" ; oda:group "Redshift" ; oda:description "Spectroscopic (Reference) Redshift"
ZP_correction = False  # http://odahub.io/ontology#Boolean ; oda:allowed_value "False", "True" ; oda:label "Zero Point Correction" ; oda:group "Redshift" ; oda:description "Zero Point Correction"

# Add in output redshift PDF
PDZ_output = False  # http://odahub.io/ontology#Boolean ; oda:allowed_value "False", "True" ; oda:label "Output: Redshift PDF" ; oda:description "Output: Redshift PDF"

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "_data_product" in inp_dic.keys():
    inp_pdic = inp_dic["_data_product"]
else:
    inp_pdic = inp_dic

for _vn in [
    "catalog_URL",
    "ab_magnitude",
    "object_type",
    "skip_sources",
    "proc_sources",
    "column_name_MW_EBV",
    "column_name_RA",
    "column_name_DEC",
    "priors",
    "column_name_Nz_prior_I",
    "column_name_Ztrue",
    "ZP_correction",
    "PDZ_output",
]:
    globals()[_vn] = type(globals()[_vn])(inp_pdic[_vn])

for _vn in ["filters_table"]:
    with open(inp_pdic[_vn], "r") as _this_fd:
        _vv = json.load(_this_fd)
        globals()[_vn] = _vv

# check filters_table inputs

# remove WISE3 & WISE4 filters (if present): problems with ModelGrids
flt = [x.split("|")[1] for x in filters_table["filter"]]
if "W3" in flt or "W4" in flt:
    f1 = []
    f2 = []
    f3 = []
    for i, x in enumerate(flt):
        if x != "W3" and x != "W4":
            f1.append(filters_table["filter"][i])
            f2.append(filters_table["flux"][i])
            f3.append(filters_table["flux_error"][i])
    filters_table = {"filter": f1, "flux": f2, "flux_error": f3}
    print("New Filter Table:", filters_table)

# at least 2 filters; same length in 'filter', 'flux', 'fluxerr'
if len(filters_table["filter"]) < 2:
    raise RuntimeError("0 or 1 filter provided; at least 2 required!")
else:
    print("Number of Filters provided =", len(filters_table["filter"]))
if len(filters_table["filter"]) != len(filters_table["flux"]):
    raise RuntimeError("Input filters/fluxes NOT correct")
if len(filters_table["filter"]) != len(filters_table["flux_error"]):
    raise RuntimeError("Input filters/flux_erros NOT correct")

# Fix some Phosphoros parameters for STARS
if object_type == "Star":
    # no galactic correction is applied
    column_name_MW_EBV = ""
    column_name_RA = ""
    column_name_DEC = ""
    # no priors, ZPC, PDZ, spec-z
    priors = "None"
    column_name_Ztrue = ""
    ZP_correction = False
    PDZ_output = False
# priors only for Galaxy & QSO
if object_type == "HighZ":
    priors = "None"

# root path for Phosphoros
workdir = os.getcwd()  #'/home/jovyan/work/photoz-euclid/'
repo_basedir = os.environ.get("BASEDIR", os.getcwd())

# path to Model Grid (MG)
if object_type == "Star":
    path_data_MG = repo_basedir + "/data/StarModelGrids/"
elif object_type == "QSO":
    path_data_MG = repo_basedir + "/data/QSOModelGrids/"
elif object_type == "HighZ":
    path_data_MG = workdir + "/tmp/ModelGrids/"
else:
    path_data_MG = repo_basedir + "/data/ModelGrids/"
print("MG for " + object_type + ":", path_data_MG)
# path to Temporary files
path_tmp = workdir + "/tmp/"
# path to Phosphoros Auxiliary Data
path_auxiliary = repo_basedir + "/Phosphoros/AuxiliaryData/"
# path to Phosphoros Outputs
path_out = path_tmp + "Input_Catalog/"
# path to Results
path_results = workdir + "/Results/"

os.makedirs(path_out, exist_ok=True)
os.makedirs(path_results, exist_ok=True)

# get the input catalog and save it into tmp/ directory as Input_Catalog.fits
input_file = path_tmp + "Input_Catalog.fits"

read_from_url = False
try:
    output = subprocess.check_output(
        "cp " + catalog_URL + " " + input_file, shell=True
    ).decode()
except:
    print("NOT a file")
    read_from_url = True

if read_from_url:
    try:
        # output = subprocess.check_output('wget -nv -O ' + input_file + catalog_URL, shell=True).decode()
        output = subprocess.check_output(
            "wget -O " + input_file + ' "' + catalog_URL + '"', shell=True
        ).decode()
    except:
        raise RuntimeError("File NOT found")

# use SFD dust map
def ComputeMWext_SFD_SkyCoord(ra, dec):
    from dustmaps.config import config
    from dustmaps.sfd import SFDQuery

    config["data_dir"] = path_tmp
    import dustmaps.sfd

    dustmaps.sfd.fetch()

    print("Use SFD map to compute EBV")
    coords = SkyCoord(ra, dec, unit="deg", frame=FK5)
    sfd = SFDQuery()
    return sfd(coords)

# use Planck dust map
def ComputeMWext_Planck_SkyCoord(ra, dec):
    from dustmaps.config import config
    from dustmaps.planck import PlanckGNILCQuery

    config["data_dir"] = path_tmp
    import dustmaps.planck

    # dustmaps.planck.fetch() # HFI dust map = 1.2Gb
    dustmaps.planck.fetch(which="GNILC")  # 350 Mb

    print("Use Planck dust map to compute EBV")
    coords = SkyCoord(ra, dec, unit="deg", frame=FK5)
    planck = PlanckGNILCQuery()
    return planck(coords)

# Phosphoros parameters to fix & options
flag_MissingPhotometry = -99
flag_UpperLimit = -99
column_name_SourceID = "OBJECT_ID"

# read input file
filename = "Input_Catalog.fits"
try:
    catalog = ascii.read(join(path_tmp, filename))
except:
    hdul = fits.open(join(path_tmp, filename))
    catalog = astropy.table.Table(hdul[1].data)

# check not empty & not larger than allowed max number of sources
if len(catalog) == 0:
    raise RuntimeError("Empty File")

n_max_sources = 500000
if object_type == "HighZ":
    n_max_sources = 200000
print("Number of Rows in the Catalog =", len(catalog))
if len(catalog) > n_max_sources:
    if (proc_sources <= 0) | (proc_sources > n_max_sources):
        raise RuntimeError(
            f"Catalog DIM: Number of objects in the catalog is {len(catalog)}, more than the allowed value ({n_max_sources})"
        )

# define Source ID column
catalog[column_name_SourceID] = np.arange(len(catalog))

# check for: flux column names; missing fluxes; negative errors
for x in filters_table["flux"]:
    if x not in catalog.columns:
        raise RuntimeError(f"FLUX: {x} column NOT in the catalog")
    catalog[x][~np.isfinite(catalog[x])] = -99
for x in filters_table["flux_error"]:
    if x not in catalog.columns:
        raise RuntimeError(f"FLUXERR: {x} column NOT in the catalog")
    catalog[x][catalog[x] < 0] = -99
    catalog[x][~np.isfinite(catalog[x])] = -99

# Galactic Extinction Correction
column_name_DustColumnDensity = "NONE"
if column_name_MW_EBV == "":
    if column_name_RA != "" and column_name_DEC != "":
        if column_name_RA not in catalog.columns:
            raise RuntimeError(
                f"RA: {column_name_RA} column NOT in the catalog"
            )
        if column_name_DEC not in catalog.columns:
            raise RuntimeError(
                f"DEC: {column_name_DEC} column NOT in the catalog"
            )
        # catalog['GAL_EBV'] = ComputeMWext_Planck_SkyCoord(catalog[column_name_RA], catalog[column_name_DEC])
        catalog["GAL_EBV"] = ComputeMWext_SFD_SkyCoord(
            catalog[column_name_RA], catalog[column_name_DEC]
        )
        column_name_DustColumnDensity = "GAL_EBV"
        print("Compute MW extinction correction from dust map")
    else:
        print("No correction for MW extinction")
else:
    if column_name_MW_EBV not in catalog.columns:
        raise RuntimeError(
            f"MW_EBV: {column_name_MW_EBV} column NOT in the catalog"
        )
    catalog[column_name_MW_EBV][~np.isfinite(catalog[column_name_MW_EBV])] = 0
    column_name_DustColumnDensity = column_name_MW_EBV
    print("Compute MW extinction correction using provided EBV")

# save the catalog
catalog.write(
    join(path_tmp, "Input_Catalog.fits"), format="fits", overwrite=True
)

# filter name to be corrected (from front-end)
name2change = {
    "CTIO": "DECam",
    "PAN-STARRS": "PANSTARRS",
    "VIRCAM": "VISTA",
    "SLOAN": "SDSS",
    "Subaru": "HSC",
}

def correct_name(filt):
    for x in name2change:
        if x == filt:
            return name2change[x]
    return filt

print("Input Filters", filters_table["filter"])
newname = []
for x in filters_table["filter"]:
    xsplit = x.split("|")
    newname.append(correct_name(xsplit[0]) + "|" + xsplit[1])

filters_table["filter"] = newname

# list of filter groups & names
filter_groups = list(
    np.unique([x.split("|")[0] for x in filters_table["filter"]])
)
filter_list = [x.replace("|", "/") for x in filters_table["filter"]]
filter_list_name = [x.split("|")[1] for x in filters_table["filter"]]
mg_groups = [x + "_MG.fits" for x in filter_groups]
gccg_groups = [x + "_GCCG.fits" for x in filter_groups]

print("Instruments", filter_groups)
print("Path + Filters", filter_list)
print("Filters", filter_list_name)
print("ModelGrids", mg_groups)
print("GCCGrids", gccg_groups)

# For Phosphoros CR Action
Phospho_para = {"input-catalog-file": "Input_Catalog.fits"}

# other Phosphoros CR parameters
Phospho_para["missing-photometry-flag"] = flag_MissingPhotometry
Phospho_para["upper-limit-use-threshold-flag"] = flag_UpperLimit

# specify columns in input catalog
Phospho_para["dust-column-density-column-name"] = column_name_DustColumnDensity
if column_name_DustColumnDensity != "NONE":
    Phospho_para["galactic-correction-coefficient-grid-file"] = (
        path_tmp + "ModelGrids/GCCGrid.txt"
    )
Phospho_para["source-id-column-name"] = column_name_SourceID

# MG to use in Phosphoros
Phospho_para["model-grid-file"] = path_tmp + "ModelGrids/ModelGrid.txt"

# for ZPC------
Phospho_ZPC_para = Phospho_para.copy()
# Phospho_ZPC_para.pop('create-output-pdf', None)
# Phospho_ZPC_para.pop('copy-columns', None)

# Parameters not in ZPC-----

# skip lines
if skip_sources > 0:
    print("skip the first {} lines".format(skip_sources))
    Phospho_para["input-skip-head"] = skip_sources
if proc_sources < len(catalog) and proc_sources > 0:
    print("process {} lines".format(proc_sources))
    Phospho_para["input-process-max"] = proc_sources

# add PDZs in output
if PDZ_output:
    Phospho_para["create-output-pdf"] = "Z"

# add Spectroscopic Redshift
if column_name_Ztrue != "":
    if column_name_Ztrue not in catalog.columns:
        raise RuntimeError(
            f"Reference Redshift: {column_name_Ztrue} column NOT in the catalog"
        )
    Phospho_para["copy-columns"] = column_name_Ztrue

# Priors
if priors == "Volume":
    print("Volume prior")
    Phospho_para["volume-prior"] = "YES"
    Phospho_para["volume-prior-effectiveness"] = 0.3
elif priors == "Redshift":
    print("Redshift prior")
    # get the I filter from I flux
    for i, x in enumerate(filters_table["flux"]):
        if x == column_name_Nz_prior_I:
            column_name_Nz_prior_I = filters_table["filter"][i].split("|")[1]
            break
    print("I filter", column_name_Nz_prior_I)

    if column_name_Nz_prior_I not in filter_list_name:
        raise RuntimeError(
            f"Redshift prior: {column_name_Nz_prior_I} filter NOT used in the catalog"
        )
    Phospho_para["Nz-prior"] = "YES"
    Phospho_para["Nz-prior_B_Filter"] = "CFHT/MegaCam.g"
    i_filter = [
        i
        for i, x in enumerate(filter_list_name)
        if x == column_name_Nz_prior_I
    ][0]
    Phospho_para["Nz-prior_I_Filter"] = filter_list[i_filter]
elif priors == "Top-Hat LF":
    Phospho_para["volume-prior"] = "YES"
    Phospho_para["volume-prior-effectiveness"] = 0.3
    print("TopHat LF prior")
    Phospho_para["luminosity-prior"] = "YES"
    Phospho_para["luminosity-function-expressed-in-magnitude"] = "YES"
    Phospho_para["luminosity-sed-group-Cosmos"] = (
        "CosmosEll,CosmosSB,CosmosSp,EuclidQSO"
    )
    Phospho_para["luminosity-function-sed-group-1"] = "Cosmos"
    Phospho_para["luminosity-function-min-z-1"] = 0
    Phospho_para["luminosity-function-max-z-1"] = 6
    if object_type == "QSO":
        print("Range [-30,-20]")
        Phospho_para["luminosity-function-curve-1"] = "TopHat_-30_-20"
    else:
        print("Range [-24,0]")
        Phospho_para["luminosity-function-curve-1"] = (
            "TopHat_-24_0"  #'TopHat_-28.000000_-8.000000'
        )
    # Phospho_para['luminosity-prior-per-mpc3']='YES'
else:
    print("No Priors")

"""
Create proper Model Grid (and GCC Grid):
1) read model grids of instruments (in .fits format) and join model grids
2) select filters and save in one .fits file
3) read the new fits model grid and save it into Phosphoros format
"""

columns_to_keep = [
    "ID",
    "Model_SED",
    "Model_RedCurve",
    "Model_EBV",
    "Model_Z",
] + filter_list_name
print("IMP: Columns to write into MG", columns_to_keep)

def read_fits(f):
    hdul = fits.open(f)
    # get data
    mg = astropy.table.Table(hdul[1].data)
    # get comments in header
    s = str(hdul[1].header["COMMENT"])
    s = s.replace("\n", "")
    s = s.replace("_new_el/", "/")
    return mg, s

def download_mg_highz(url):
    tmp_path = path_data_MG + "highz.tar.gz "
    output = subprocess.check_output(
        "wget -nv -O " + tmp_path + url, shell=True
    ).decode()
    output = subprocess.check_output(
        "tar -xvf " + tmp_path + "-C " + path_data_MG, shell=True
    ).decode()
    return

# MG_web = ['HSTACS','HSTWFC3','JWSTMIRI','JWSTNIRCam','HSCinter','SuprimeCam']
MG_web = ["HST", "JWST", "SuprimeCam"]

def download_mg_gal(group):
    url = "https://www.astro.unige.ch/~tucci/Phosphoros/" + group + ".tar.gz"
    tmp_path = path_data_MG + group + ".tar.gz "
    output = subprocess.check_output(
        "wget -nv -O " + tmp_path + url, shell=True
    ).decode()
    output = subprocess.check_output(
        "tar -xvf " + tmp_path + "-C " + path_data_MG, shell=True
    ).decode()
    return

if object_type == "HighZ":
    download_mg_highz(
        "https://www.astro.unige.ch/~tucci/Phosphoros/highz.tar.gz"
    )

if object_type == "Galaxy":
    for xg in filter_groups:
        if xg in MG_web:
            print("Download MGs from WEB:", xg)
            download_mg_gal(xg)

# step 1
for i, x in enumerate(mg_groups):
    print(x)
    mg, comment = read_fits(join(path_data_MG, filter_groups[i] + "_MG.fits"))
    if i == 0:
        mg_full = mg.copy()
    else:
        # check MG are from the same ParameterSpace
        if len(mg_full) != len(mg):
            raise RuntimeError(
                f"Error: ModelGrid {x} not from the same Parameter Space"
            )
        # join the two tables (matched on columns with same name)
        mg_full = astropy.table.join(mg_full, mg)

# only if MW correction is applied
if column_name_DustColumnDensity != "NONE":
    for i, x in enumerate(gccg_groups):
        gccg, comment2 = read_fits(
            join(path_data_MG, filter_groups[i] + "_GCCG.fits")
        )
        if i == 0:
            gccg_full = gccg.copy()
        else:
            # check GCCGs are from the same ParameterSpace
            if len(gccg_full) != len(gccg):
                raise RuntimeError(
                    f"Error: GCCGrid {x} not from the same Parameter Space"
                )
            # join the two tables (matched on columns with same name)
            gccg_full = astropy.table.join(gccg_full, gccg)

# step 2
def write_fits(d, comment, f):
    table_hdu = fits.BinTableHDU(data=d)
    table_hdu.header["COMMENT"] = comment
    # print(table_hdu.header)
    table_hdu.writeto(join(path_tmp, f), overwrite=True)
    return

# select columns to use
mg_full = mg_full[columns_to_keep]
# change the name of filters to include the path
for x, y in zip(filter_list_name, filter_list):
    mg_full.rename_column(x, y)

print(mg_full.info)
# create a new .fits file
write_fits(mg_full, comment, "ModelGrid.fits")

# only if MW correction is applied
if column_name_DustColumnDensity != "NONE":
    # select columns to use
    gccg_full = gccg_full[columns_to_keep]
    # change the name of filters to include the path
    for x, y in zip(filter_list_name, filter_list):
        gccg_full.rename_column(x, y)

    print(gccg_full.info)
    # create a new .fits file
    write_fits(gccg_full, comment2, "GCCGrid.fits")

# step 3
if object_type == "Star":
    conf_file = "/Phosphoros/config/RenkuTest.FMG.STAR.conf"
elif object_type == "QSO":
    conf_file = "/Phosphoros/config/RenkuTest.FMG.QSO.conf"
elif object_type == "HighZ":
    conf_file = "/Phosphoros/config/RenkuTest.FMG.HighZ.conf"
else:
    conf_file = "/Phosphoros/config/RenkuTest.FMG.conf"
# create MG in Phosphoros format
run_FMG = "Phosphoros FMG"
run_FMG += " --config-file " + repo_basedir + conf_file
run_FMG += (
    " --phosphoros-root "
    + repo_basedir
    + " --intermediate-products-dir "
    + workdir
    + " --catalog-type tmp"
)
run_FMG += (
    " --aux-data-dir "
    + path_auxiliary
    + " --input-catalog-file "
    + path_tmp
    + "ModelGrid.fits"
)
for x in mg_full.columns:
    if (
        x in filter_list
    ):  # must be in the same order that found in the .fits file
        run_FMG += " --filter-name " + x

print("\n", run_FMG)
output = subprocess.check_output(run_FMG, shell=True).decode()

# only if MW correction is applied
if column_name_DustColumnDensity != "NONE":
    # create GCCG in Phosphoros format
    run_FMG = "Phosphoros FMG"
    run_FMG += " --config-file " + repo_basedir + conf_file
    run_FMG += (
        " --phosphoros-root "
        + repo_basedir
        + " --intermediate-products-dir "
        + workdir
        + " --catalog-type tmp"
    )
    run_FMG += (
        " --aux-data-dir "
        + path_auxiliary
        + " --input-catalog-file "
        + path_tmp
        + "GCCGrid.fits"
    )
    run_FMG += " --output-model-grid " + path_tmp + "ModelGrids/GCCGrid.txt"
    for x in gccg_full.columns:
        if (
            x in filter_list
        ):  # must be in the same order that found in the .fits file
            run_FMG += " --filter-name " + x

    print("\n", run_FMG)
    output = subprocess.check_output(run_FMG, shell=True).decode()

# Mapping Filters to Catalog Columns (write file: 'filter_mapping.txt')

# map columns to filters
map_filters = astropy.table.Table(
    np.chararray((len(filter_list), 6)),
    names=(
        "#Filter",
        "Flux",
        "Error",
        "Upper Limit/error ratio",
        "Convert from MAG",
        "Filter Shift Column",
    ),
)

flux = []
fluxerr = []
filt = []
for g in filter_groups:
    filt += ",".join([x for x in filter_list if x.startswith(g)]).split(",")
    flux += ",".join(
        [
            filters_table["flux"][i]
            for i, x in enumerate(filters_table["filter"])
            if x.startswith(g)
        ]
    ).split(",")
    fluxerr += ",".join(
        [
            filters_table["flux_error"][i]
            for i, x in enumerate(filters_table["filter"])
            if x.startswith(g)
        ]
    ).split(",")

map_filters["#Filter"] = filt
map_filters["Flux"] = flux
map_filters["Error"] = fluxerr

map_filters["Upper Limit/error ratio"] = 3

if ab_magnitude == "Magnitude":
    map_filters["Convert from MAG"] = 1
else:
    map_filters["Convert from MAG"] = 0

map_filters["Filter Shift Column"] = ["NONE"] * len(filter_list)

out_map_filters = join(path_tmp, "filter_mapping.txt")
map_filters.write(out_map_filters, format="ascii", overwrite=True)

# read file in order to save it as text (string)
# file_txt = open(join(path_tmp,'filter_mapping.txt'), 'r')
# filter_mapping_list = file_txt.readlines()
# filter_mapping_list

map_filters

zpc_table = astropy.table.Table()

# Run Zero-Point Correction, if required
if (ZP_correction) & (column_name_Ztrue != ""):
    print("Compute Zero-Point Correction")
    filename_zpc = "ZeroPointCorrection.txt"

    run_ZPC = (
        "Phosphoros CPC --config-file "
        + repo_basedir
        + "/Phosphoros/config/RenkuTest.CPC.conf"
    )
    run_ZPC += (
        " --phosphoros-root "
        + repo_basedir
        + " --catalogs-dir "
        + workdir
        + " --intermediate-products-dir "
        + workdir
    )
    run_ZPC += " --aux-data-dir " + path_auxiliary
    # CPC parameters
    run_ZPC += " --output-phot-corr-file " + filename_zpc
    run_ZPC += (
        " --phot-corr-iter-no "
        + "3"
        + " --phot-corr-selection-method "
        + "WEIGHTED_MEDIAN"
    )
    run_ZPC += " --spec-z-column-name " + column_name_Ztrue
    print(run_ZPC, "\n")

    CR_options = ""
    for x in Phospho_ZPC_para:
        CR_options += " --" + x + " " + str(Phospho_ZPC_para[x])
    print("\n", CR_options)

    output = subprocess.check_output(run_ZPC + CR_options, shell=True).decode()

    # add ZPC parameters to CR
    Phospho_para["enable-photometric-correction"] = "YES"
    Phospho_para["photometric-correction-file"] = filename_zpc

    # ZPC table to save
    zpc_table = ascii.read(join(path_tmp, filename_zpc))
    zpc_table.rename_column("col1", "Filter")
    zpc_table.rename_column("col2", "Correction")
    print(zpc_table)

# import time
# start_time = time.time()

run_Phosphoros = (
    "Phosphoros CR --config-file "
    + repo_basedir
    + "/Phosphoros/config/RenkuTest.CR.conf"
)
run_Phosphoros += (
    " --phosphoros-root "
    + repo_basedir
    + " --catalogs-dir "
    + workdir
    + " --intermediate-products-dir "
    + workdir
)
run_Phosphoros += (
    " --aux-data-dir " + path_auxiliary + " --results-dir " + workdir
)
print(run_Phosphoros, "\n")

CR_options = ""
for x in Phospho_para:
    # print(x)
    CR_options += " --" + x + " " + str(Phospho_para[x])
print("\n", CR_options)

output = subprocess.check_output(
    run_Phosphoros + CR_options, shell=True
).decode()

# print("--- %s seconds ---" % (time.time() - start_time))

out_fn = join(path_out, "phz_cat.fits")
cat = astropy.table.Table(fits.open(out_fn)[1].data)
print(len(cat))

if PDZ_output:
    print("zPDF & stats")
    # compute quantiles of zPDFs
    pdf = np.array(cat["Z-1D-PDF"])
    zed = np.arange(0, 6.02, 0.02)

    area = np.sum(pdf, axis=1)
    pdf = (pdf.T / area).T
    cdf = np.cumsum(pdf, axis=1)

    quant = np.zeros((pdf.shape[0], 5))
    for i in range(pdf.shape[0]):
        quant[i, :] = zed[
            np.searchsorted(cdf[i, :], [0.05, 0.25, 0.5, 0.75, 0.95])
        ]
    cat["Z-Quantiles"] = quant

    # histogram of medians
    # med = cat['Z-Quantiles']
    # fig, ax = plt.subplots()
    # ax.hist(med[:,2], bins=100, range=(0,6))

    # save the catalog
    cat.write(join(path_out, "phz_cat.fits"), format="fits", overwrite=True)

# histogram
zmin = np.min(cat["Z"])
if zmin > 0:
    zmin = 0
zmax = np.max(cat["Z"])
nsteps = int((zmax - zmin) / 0.1)

# name histogram output
out_photoZ = join(path_out, "plot_photoZ.png")

def z_dist(ax, z):
    ax.hist(z, bins=nsteps, range=(zmin, zmax))
    ax.set_xlabel("Photometric Redshift")
    ax.set_ylabel("Number of Objects")
    return

def chi2_dist(ax, chi2):
    ax.hist(chi2, bins=50)  # , range=(zmin,zmax))
    ax.set_xlabel("Phosphoros Chi2 (STAR ModelGrid)")
    ax.set_ylabel("Number of Objects")
    return

if column_name_Ztrue == "":
    # STAR: chi2 distribution
    if object_type == "Star":
        fig, ax = plt.subplots(1)
        chi2_dist(ax, np.log10(-2 * cat["Posterior-Log"]))
    else:
        # GALAXY,QSO,HighZ: Redshift distribution
        fig, ax = plt.subplots(1)
        z_dist(ax, cat["Z"])

    plt.savefig(out_photoZ, format="png")
else:
    if column_name_Ztrue not in catalog.columns:
        print(
            f"WARNING: {column_name_Ztrue} column NOT in the catalog; no PHZ vs SPECZ plot"
        )
        raise RuntimeError(
            f"Reference Z: {column_name_Ztrue} column NOT in the catalog"
        )
    else:
        fig, axs = plt.subplots(3, figsize=(6, 12), height_ratios=[1, 2, 1])

        # photo-Z distribution
        ax = axs[0]
        z_dist(ax, cat["Z"])

        # photo-Z vs Spec-Z
        ax = axs[1]
        zmax2 = max(max(cat[column_name_Ztrue]), zmax) + 0.1
        ax.set_xlim(0, zmax2)
        ax.set_ylim(0, zmax2)
        ax.plot([0, zmax2], [0, zmax2], "-", c="black", alpha=0.5)
        ax.plot(
            [0, zmax2],
            [0.15, zmax2 + 0.15 * (1 + zmax2)],
            ":",
            c="black",
            alpha=0.5,
        )
        ax.plot(
            [0, zmax2],
            [-0.15, zmax2 - 0.15 * (1 + zmax2)],
            ":",
            c="black",
            alpha=0.5,
        )
        ax.scatter(cat["Z"], cat[column_name_Ztrue], s=5, alpha=0.3)
        ax.set_xlabel("Photometric Redshift")
        ax.set_ylabel("Spectroscopic Redshift")

        # histogram (photoZ - specZ)/(1+speZ)
        dz = (cat["Z"] - cat[column_name_Ztrue]) / (1 + cat[column_name_Ztrue])
        mean = np.nanmean(dz)
        median = np.nanmedian(dz)
        mad = np.nanmedian(np.abs(dz - median))
        nmad = 1.4826 * mad
        sigma = np.nanstd(dz)
        outliers = len(dz[np.abs(dz) > 0.15]) / len(dz)
        # print(mean,median,mad,sigma,outliers)

        ax = axs[2]
        ndz = ax.hist(dz, bins=100, range=(-2, 2))
        ax.plot(
            [-0.15, -0.15], [0, 0.9 * max(ndz[0])], ":", c="red", alpha=0.5
        )
        ax.plot([0.15, 0.15], [0, 0.9 * max(ndz[0])], ":", c="red", alpha=0.5)

        ax.text(max(0.5 * ndz[1]), 0.8 * max(ndz[0]), f"Mean: {mean:.3f}")
        ax.text(max(0.5 * ndz[1]), 0.7 * max(ndz[0]), f"Median: {median:.3f}")
        ax.text(max(0.5 * ndz[1]), 0.6 * max(ndz[0]), f"NMAD: {nmad:.3f}")
        ax.text(max(0.5 * ndz[1]), 0.5 * max(ndz[0]), f"Sigma: {sigma:.3f}")
        ax.text(
            max(0.5 * ndz[1]),
            0.4 * max(ndz[0]),
            f"Outliers: {100*outliers:.1f}%",
        )

        ax.set_xlabel("(PhotoZ - SpecZ)/(1 + SpecZ)")
        ax.set_ylabel("Number of Objects")

        plt.savefig(out_photoZ, format="png")

output_catalog = BinaryProduct.from_file(
    out_fn, name="catalog.fits"
)  # http://odahub.io/ontology#ODABinaryProduct
z_hist = PictureProduct.from_file(
    out_photoZ
)  # http://odahub.io/ontology#ODAPictureProduct
# filter_mapping_output=filter_mapping_str # http://odahub.io/ontology#ODATextProduct
output_zpc = ODAAstropyTable(
    zpc_table
)  # http://odahub.io/ontology#ODAAstropyTable

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    (
        "out_Run_Phosphoros_on_Catalog_output_catalog",
        "output_catalog_galaxy.output",
        output_catalog,
    )
)
_oda_outs.append(
    ("out_Run_Phosphoros_on_Catalog_z_hist", "z_hist_galaxy.output", z_hist)
)
_oda_outs.append(
    (
        "out_Run_Phosphoros_on_Catalog_output_zpc",
        "output_zpc_galaxy.output",
        output_zpc,
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
