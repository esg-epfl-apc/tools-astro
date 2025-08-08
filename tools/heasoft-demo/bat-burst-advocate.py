#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python

# This script is generated with nb2galaxy

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

indir = ""  # oda:String; oda:description "Input observation directory name"
outdir = ""  # oda:String; oda:description "Output products directory name"
ra = "BAT"  # oda:String,oda:optional; oda:description "RA of GRB (or BAT,XRT,BLIND) [deg]"
dec = "BAT"  # oda:String,oda:optional; oda:description "Dec of GRB (or BAT,XRT,BLIND) [deg]"
trigtime = 0.0  # oda:Float,oda:optional; oda:description "Start time of trigger interval (MET) [s]"
trigstop = 0.0  # oda:Float,oda:optional; oda:description "Stop time of trigger interval (MET) [s]"
backstrt = 0.0  # oda:Float,oda:optional; oda:description "Start time of background interval (MET) [s]"
backstop = 0.0  # oda:Float,oda:optional; oda:description "Stop time of background interval (MET) [s]"
shortfix = "scaledmap,expand"  # oda:String,oda:optional; oda:description "How to deal with short trigger time errors?"
tbkgsub = True  # oda:Boolean,oda:optional; oda:description "Perform background subtraction for T50/T90 duration estimates?"
tnear = 2000.0  # oda:Float,oda:optional ; oda:lower_limit 0; oda:description "Number of seconds around trigger time to include?"
tbinmax = 1000.0  # oda:Float,oda:optional ; oda:lower_limit 0; oda:description "Maximum time bin size for light curves when searching for T50/T90?"
pcodethresh = 0.0  # oda:Float,oda:optional ; oda:lower_limit 0.0 ; oda:upper_limit 1.0; oda:description "Minimum allowed partial coding fraction"
imgpcodethresh = 0.05  # oda:Float,oda:optional ; oda:lower_limit 0.0 ; oda:upper_limit 1.0; oda:description "Minimum allowed partial coding fraction for images"
aperture = "CALDB:FLUX"  # oda:String,oda:optional; oda:description "User-requested aperture file (or INDEF for default)"
date_obs = "INDEF"  # oda:String,oda:optional; oda:description "User override of DATE-OBS keyword (or INDEF)"
extractor = "fextract-events"  # oda:String,oda:optional; oda:description "Event extractor to use (extractor or fextract-events)"
chatter = 1  # oda:Integer,oda:optional ; oda:lower_limit 0 ; oda:upper_limit 5; oda:description "Verbosity level"
clobber = True  # oda:Boolean,oda:optional; oda:description "Clobber existing output file?"
history = (
    True  # oda:Boolean,oda:optional; oda:description "Write history block?"
)
mode = "ql"  # oda:String,oda:optional; oda:description ""

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "C_data_product_" in inp_dic.keys():
    inp_pdic = inp_dic["C_data_product_"]
else:
    inp_pdic = inp_dic
indir = str(inp_pdic["indir"])
outdir = str(inp_pdic["outdir"])

ra = str(inp_pdic["ra"]) if inp_pdic.get("ra", None) is not None else None

dec = str(inp_pdic["dec"]) if inp_pdic.get("dec", None) is not None else None

trigtime = (
    float(inp_pdic["trigtime"])
    if inp_pdic.get("trigtime", None) is not None
    else None
)

trigstop = (
    float(inp_pdic["trigstop"])
    if inp_pdic.get("trigstop", None) is not None
    else None
)

backstrt = (
    float(inp_pdic["backstrt"])
    if inp_pdic.get("backstrt", None) is not None
    else None
)

backstop = (
    float(inp_pdic["backstop"])
    if inp_pdic.get("backstop", None) is not None
    else None
)

shortfix = (
    str(inp_pdic["shortfix"])
    if inp_pdic.get("shortfix", None) is not None
    else None
)

tbkgsub = (
    bool(inp_pdic["tbkgsub"])
    if inp_pdic.get("tbkgsub", None) is not None
    else None
)

tnear = (
    float(inp_pdic["tnear"])
    if inp_pdic.get("tnear", None) is not None
    else None
)

tbinmax = (
    float(inp_pdic["tbinmax"])
    if inp_pdic.get("tbinmax", None) is not None
    else None
)

pcodethresh = (
    float(inp_pdic["pcodethresh"])
    if inp_pdic.get("pcodethresh", None) is not None
    else None
)

imgpcodethresh = (
    float(inp_pdic["imgpcodethresh"])
    if inp_pdic.get("imgpcodethresh", None) is not None
    else None
)

aperture = (
    str(inp_pdic["aperture"])
    if inp_pdic.get("aperture", None) is not None
    else None
)

date_obs = (
    str(inp_pdic["date_obs"])
    if inp_pdic.get("date_obs", None) is not None
    else None
)

extractor = (
    str(inp_pdic["extractor"])
    if inp_pdic.get("extractor", None) is not None
    else None
)

chatter = (
    int(inp_pdic["chatter"])
    if inp_pdic.get("chatter", None) is not None
    else None
)

clobber = (
    bool(inp_pdic["clobber"])
    if inp_pdic.get("clobber", None) is not None
    else None
)

history = (
    bool(inp_pdic["history"])
    if inp_pdic.get("history", None) is not None
    else None
)

mode = (
    str(inp_pdic["mode"]) if inp_pdic.get("mode", None) is not None else None
)

# temporary fix to support python 3.13
get_ipython().system("pip install git+https://github.com/oda-hub/oda_api.git")   # noqa: F821

command = "bat-burst-advocate"
params_names = "infile,indir,outdir,ra,dec,trigtime,trigstop,backstrt,backstop,shortfix,tbkgsub,tnear,tbinmax,pcodethresh,imgpcodethresh,aperture,date_obs,extractor,chatter,clobber,history,mode"
optional_file_params = ""

outfile = "output"

def format(p):
    if type(p) == bool:
        return "YES" if p else "NO"
    # if type(p) == str and ' ' in p:
    #    return f'"{p}"'
    return p

optional_file_params = [
    p for p in optional_file_params.split(",") if len(p) > 0
]

l = locals()
params = {p: format(l[p]) for p in params_names.split(",")}

for param_name in optional_file_params:
    if str(params[param_name]).lower() == "none":
        print(param_name, "is not set")
        del params[param_name]

cline_pars = [f"{p}={v}" for p, v in params.items()]

cline_pars

import subprocess

command_line = [command] + cline_pars
process = subprocess.Popen(
    command_line, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
)
stdout, stderr = process.communicate(input="\n\n\n\n")
if process.returncode != 0:
    raise Exception(stderr)
else:
    print(stdout)

import os

from oda_api.data_products import BinaryProduct

result = None

if os.path.isfile(outfile):
    result = BinaryProduct.from_file(outfile, name="result")

result = result  # http://odahub.io/ontology#ODABinaryProduct

# output gathering
_galaxy_meta_data = {}
_oda_outs = []
_oda_outs.append(
    ("out_bat_burst_advocate_result", "result_galaxy.output", result)
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
