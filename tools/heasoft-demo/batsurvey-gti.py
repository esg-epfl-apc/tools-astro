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
dphfiles = "NONE"  # oda:POSIXPath; oda:description "Input survey DPH files names (or @filename)"
outdir = ""  # oda:String; oda:description "Output results directory"
gtifile = "NONE"  # oda:POSIXPath,oda:optional; oda:description "Name of user-requested good time interval file (or NONE)"
totexpo = (
    -1.0
)  # oda:Float,oda:optional; oda:description "Total exposure of input DPH data [sec] (set upon exit)"
goodexpo = (
    -1.0
)  # oda:Float,oda:optional; oda:description "Total good time [sec] (set upon exit)"
ngti = (
    -1
)  # oda:Integer,oda:optional; oda:description "Total number of good times (set upon exit)"
filters = "all"  # oda:String,oda:optional; oda:description "Filters to apply (or 'all')"
elimits = "14-195"  # oda:String,oda:optional; oda:description "Energy limits when computing DPH light curve? [kev]"
rateminthresh = 3000.0  # oda:Float,oda:optional ; oda:lower_limit 0; oda:description "Lower survey rate threshold? [ct/s]"
ratemaxthresh = 12000.0  # oda:Float,oda:optional ; oda:lower_limit 0; oda:description "Upper survey rate threshold? [ct/s]"
detthresh = 24000  # oda:Integer,oda:optional ; oda:lower_limit 0; oda:description "Minimum number of detectors?"
ra = "NONE"  # oda:String,oda:optional; oda:description "R.A. of source [deg] for occultation check (or NONE)"
dec = "NONE"  # oda:String,oda:optional; oda:description "Dec. of source [deg] for occultation check (or NONE)"
pcodethresh = 0.1  # oda:Float,oda:optional ; oda:lower_limit 0; oda:description "Partial coding threshold for occultation check [fraction]"
batoccultgti_opts = ""  # oda:String,oda:optional; oda:description "Additional options to batoccultgti (empty string allowed)"
stlossfcnthresh = 1e-09  # oda:Float,oda:optional ; oda:lower_limit 0; oda:description "Star tracker loss function upper threshold"
filtexpr = "NONE"  # oda:String,oda:optional; oda:description "Expression to apply to the filter file"
saofiltexpr = "ELV > 30.0"  # oda:String,oda:optional; oda:description "Expression to apply to the prefilter file"
dphfiltexpr = "DATA_FLAGS == 0"  # oda:String,oda:optional; oda:description "Expression to apply to the DPH files"
surveyslop = 5.0  # oda:Float,oda:optional; oda:description "Duration to add to begin and end of each DPH [sec]"
sepdph = True  # oda:Boolean,oda:optional; oda:description "Assign each survey DPH its own GTI"
maxgtigap = 5.0  # oda:Float,oda:optional ; oda:lower_limit 0; oda:description "Maximum gap that can be joined [sec]"
mingti = 30.0  # oda:Float,oda:optional ; oda:lower_limit 0; oda:description "Minimum good time interval size allowed [sec]"
mjdrefi = 51910  # oda:Integer,oda:optional ; oda:lower_limit 0; oda:description "Integer part of MJDREF"
mjdreff = 0.00074287037  # oda:Float,oda:optional ; oda:lower_limit 0; oda:description "Fractional part of MJDREF"
clobber = True  # oda:Boolean,oda:optional; oda:description "Overwrite existing output file?"
chatter = 2  # oda:Integer,oda:optional ; oda:lower_limit 0 ; oda:upper_limit 5; oda:description "Verbosity level"
history = True  # oda:Boolean,oda:optional; oda:description "Write HISTORY keywords in copied HDU?"
mode = "ql"  # oda:String,oda:optional; oda:description "Mode"

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "C_data_product_" in inp_dic.keys():
    inp_pdic = inp_dic["C_data_product_"]
else:
    inp_pdic = inp_dic
indir = str(inp_pdic["indir"])
dphfiles = str(inp_pdic["dphfiles"])
outdir = str(inp_pdic["outdir"])

gtifile = (
    str(inp_pdic["gtifile"])
    if inp_pdic.get("gtifile", None) is not None
    else None
)

totexpo = (
    float(inp_pdic["totexpo"])
    if inp_pdic.get("totexpo", None) is not None
    else None
)

goodexpo = (
    float(inp_pdic["goodexpo"])
    if inp_pdic.get("goodexpo", None) is not None
    else None
)

ngti = (
    int(inp_pdic["ngti"]) if inp_pdic.get("ngti", None) is not None else None
)

filters = (
    str(inp_pdic["filters"])
    if inp_pdic.get("filters", None) is not None
    else None
)

elimits = (
    str(inp_pdic["elimits"])
    if inp_pdic.get("elimits", None) is not None
    else None
)

rateminthresh = (
    float(inp_pdic["rateminthresh"])
    if inp_pdic.get("rateminthresh", None) is not None
    else None
)

ratemaxthresh = (
    float(inp_pdic["ratemaxthresh"])
    if inp_pdic.get("ratemaxthresh", None) is not None
    else None
)

detthresh = (
    int(inp_pdic["detthresh"])
    if inp_pdic.get("detthresh", None) is not None
    else None
)

ra = str(inp_pdic["ra"]) if inp_pdic.get("ra", None) is not None else None

dec = str(inp_pdic["dec"]) if inp_pdic.get("dec", None) is not None else None

pcodethresh = (
    float(inp_pdic["pcodethresh"])
    if inp_pdic.get("pcodethresh", None) is not None
    else None
)

batoccultgti_opts = (
    str(inp_pdic["batoccultgti_opts"])
    if inp_pdic.get("batoccultgti_opts", None) is not None
    else None
)

stlossfcnthresh = (
    float(inp_pdic["stlossfcnthresh"])
    if inp_pdic.get("stlossfcnthresh", None) is not None
    else None
)

filtexpr = (
    str(inp_pdic["filtexpr"])
    if inp_pdic.get("filtexpr", None) is not None
    else None
)

saofiltexpr = (
    str(inp_pdic["saofiltexpr"])
    if inp_pdic.get("saofiltexpr", None) is not None
    else None
)

dphfiltexpr = (
    str(inp_pdic["dphfiltexpr"])
    if inp_pdic.get("dphfiltexpr", None) is not None
    else None
)

surveyslop = (
    float(inp_pdic["surveyslop"])
    if inp_pdic.get("surveyslop", None) is not None
    else None
)

sepdph = (
    bool(inp_pdic["sepdph"])
    if inp_pdic.get("sepdph", None) is not None
    else None
)

maxgtigap = (
    float(inp_pdic["maxgtigap"])
    if inp_pdic.get("maxgtigap", None) is not None
    else None
)

mingti = (
    float(inp_pdic["mingti"])
    if inp_pdic.get("mingti", None) is not None
    else None
)

mjdrefi = (
    int(inp_pdic["mjdrefi"])
    if inp_pdic.get("mjdrefi", None) is not None
    else None
)

mjdreff = (
    float(inp_pdic["mjdreff"])
    if inp_pdic.get("mjdreff", None) is not None
    else None
)

clobber = (
    bool(inp_pdic["clobber"])
    if inp_pdic.get("clobber", None) is not None
    else None
)

chatter = (
    int(inp_pdic["chatter"])
    if inp_pdic.get("chatter", None) is not None
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

command = "batsurvey-gti"
params_names = "infile,indir,dphfiles,outdir,gtifile,totexpo,goodexpo,ngti,filters,elimits,rateminthresh,ratemaxthresh,detthresh,ra,dec,pcodethresh,batoccultgti_opts,stlossfcnthresh,filtexpr,saofiltexpr,dphfiltexpr,surveyslop,sepdph,maxgtigap,mingti,mjdrefi,mjdreff,clobber,chatter,history,mode"
optional_file_params = "gtifile"

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
_oda_outs.append(("out_batsurvey_gti_result", "result_galaxy.output", result))

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
