#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python

# This script is generated with nb2galaxy

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

infile = "NONE"  # oda:POSIXPath; oda:description "Input event file name(s), or @file.list"
outtype = ""  # oda:String ; oda:allowed_value "LC","PHA","DPH","DPI","DPITAB"; oda:description "Make light curve or spectrum (LC, PHA, DPH, DPI or DPITAB)"
timedel = 1.0  # oda:Float ; oda:lower_limit 0; oda:description "Histogram time bin size [sec]"
timebinalg = "uniform"  # oda:String ; oda:allowed_value "uniform","snr","highsnr","gti","infile","matchlc"; oda:description "Time binning algorithm (uniform,snr,highsnr,gti,infile,matchlc)"
energybins = "-"  # oda:String; oda:description "Energy bin list (comma separated [keV], file name or INFILE)"
gtifile = "NONE"  # oda:POSIXPath,oda:optional; oda:description "Good time interval file (or NONE)"
ecol = "ENERGY"  # oda:String,oda:optional ; oda:allowed_value "ENERGY","PI","PHA"; oda:description "Energy column to accumulate (ENERGY, PI or PHA)"
weighted = "INDEF"  # oda:String,oda:optional ; oda:allowed_value "yes","no","INDEF"; oda:description "Apply mask weighting or not (yes,no,INDEF)"
outunits = "INDEF"  # oda:String,oda:optional ; oda:allowed_value "RATE","COUNTS"; oda:description "Output units (RATE or COUNTS)"
timepixr = (
    -1.0
)  # oda:Float,oda:optional ; oda:lower_limit -1 ; oda:upper_limit 1; oda:description "Light curve time bin reference point (0.0=begin; 0.5=center; -1=default)"
maskwt = "NONE"  # oda:POSIXPath,oda:optional; oda:description "Mask weight data, if not in infile (or NONE)"
tstart = "INDEF"  # oda:String,oda:optional; oda:description "Histogram start time (MET, or INDEF)"
tstop = "INDEF"  # oda:String,oda:optional; oda:description "Histogram stop time (MET, or INDEF)"
snrthresh = 6.0  # oda:Float,oda:optional ; oda:lower_limit 0.01; oda:description "Signal to noise ratio"
detmask = "NONE"  # oda:POSIXPath,oda:optional; oda:description "Detector quality mask file name (or NONE)"
tcol = "TIME"  # oda:String,oda:optional; oda:description "Name of TIME column"
countscol = "DPH_COUNTS"  # oda:String,oda:optional; oda:description "Name of DPH counts column"
xcol = "DETX"  # oda:String,oda:optional; oda:description "Name of X column for image / DPH creation"
ycol = "DETY"  # oda:String,oda:optional; oda:description "Name of Y column for image / DPH creation"
maskwtcol = "MASK_WEIGHT"  # oda:String,oda:optional; oda:description "Name of mask weight column"
ebinquant = 0.1  # oda:Float,oda:optional ; oda:lower_limit 1e-6 ; oda:upper_limit 1e6; oda:description "Default energy bin quantization [keV]"
delzeroes = True  # oda:Boolean,oda:optional; oda:description "Delete time bins with zero flux?"
minfracexp = 0.1  # oda:Float,oda:optional ; oda:lower_limit 0 ; oda:upper_limit 1.0; oda:description "Minimum fractional exposure per time bin"
min_dph_frac_overlap = 0.999  # oda:Float,oda:optional ; oda:lower_limit 0 ; oda:upper_limit 1.0; oda:description "Minimum per-DPH overlap [fraction]"
min_dph_time_overlap = 0.0  # oda:Float,oda:optional ; oda:lower_limit 0; oda:description "Minimum per-DPH time overlap [sec]"
max_dph_time_nonoverlap = 0.5  # oda:Float,oda:optional ; oda:lower_limit 0; oda:description "Maximum per-DPH non-overlap time [sec]"
buffersize = 16384  # oda:Integer,oda:optional ; oda:lower_limit 8 ; oda:upper_limit 33554432; oda:description "Input read buffer size"
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
infile = str(inp_pdic["infile"])
outtype = str(inp_pdic["outtype"])
timedel = float(inp_pdic["timedel"])
timebinalg = str(inp_pdic["timebinalg"])
energybins = str(inp_pdic["energybins"])

gtifile = (
    str(inp_pdic["gtifile"])
    if inp_pdic.get("gtifile", None) is not None
    else None
)

ecol = (
    str(inp_pdic["ecol"]) if inp_pdic.get("ecol", None) is not None else None
)

weighted = (
    str(inp_pdic["weighted"])
    if inp_pdic.get("weighted", None) is not None
    else None
)

outunits = (
    str(inp_pdic["outunits"])
    if inp_pdic.get("outunits", None) is not None
    else None
)

timepixr = (
    float(inp_pdic["timepixr"])
    if inp_pdic.get("timepixr", None) is not None
    else None
)

maskwt = (
    str(inp_pdic["maskwt"])
    if inp_pdic.get("maskwt", None) is not None
    else None
)

tstart = (
    str(inp_pdic["tstart"])
    if inp_pdic.get("tstart", None) is not None
    else None
)

tstop = (
    str(inp_pdic["tstop"]) if inp_pdic.get("tstop", None) is not None else None
)

snrthresh = (
    float(inp_pdic["snrthresh"])
    if inp_pdic.get("snrthresh", None) is not None
    else None
)

detmask = (
    str(inp_pdic["detmask"])
    if inp_pdic.get("detmask", None) is not None
    else None
)

tcol = (
    str(inp_pdic["tcol"]) if inp_pdic.get("tcol", None) is not None else None
)

countscol = (
    str(inp_pdic["countscol"])
    if inp_pdic.get("countscol", None) is not None
    else None
)

xcol = (
    str(inp_pdic["xcol"]) if inp_pdic.get("xcol", None) is not None else None
)

ycol = (
    str(inp_pdic["ycol"]) if inp_pdic.get("ycol", None) is not None else None
)

maskwtcol = (
    str(inp_pdic["maskwtcol"])
    if inp_pdic.get("maskwtcol", None) is not None
    else None
)

ebinquant = (
    float(inp_pdic["ebinquant"])
    if inp_pdic.get("ebinquant", None) is not None
    else None
)

delzeroes = (
    bool(inp_pdic["delzeroes"])
    if inp_pdic.get("delzeroes", None) is not None
    else None
)

minfracexp = (
    float(inp_pdic["minfracexp"])
    if inp_pdic.get("minfracexp", None) is not None
    else None
)

min_dph_frac_overlap = (
    float(inp_pdic["min_dph_frac_overlap"])
    if inp_pdic.get("min_dph_frac_overlap", None) is not None
    else None
)

min_dph_time_overlap = (
    float(inp_pdic["min_dph_time_overlap"])
    if inp_pdic.get("min_dph_time_overlap", None) is not None
    else None
)

max_dph_time_nonoverlap = (
    float(inp_pdic["max_dph_time_nonoverlap"])
    if inp_pdic.get("max_dph_time_nonoverlap", None) is not None
    else None
)

buffersize = (
    int(inp_pdic["buffersize"])
    if inp_pdic.get("buffersize", None) is not None
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

command = "batbinevt"
params_names = "infile,infile,outfile,outtype,timedel,timebinalg,energybins,gtifile,ecol,weighted,outunits,timepixr,maskwt,tstart,tstop,snrthresh,detmask,tcol,countscol,xcol,ycol,maskwtcol,ebinquant,delzeroes,minfracexp,min_dph_frac_overlap,min_dph_time_overlap,max_dph_time_nonoverlap,buffersize,clobber,chatter,history,mode"
optional_file_params = "gtifile,maskwt,detmask"

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
_oda_outs.append(("out_batbinevt_result", "result_galaxy.output", result))

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
