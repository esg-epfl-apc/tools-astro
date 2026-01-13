#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python

# This script is generated with nb2galaxy

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

gtifile = "NONE"  # oda:POSIXPath; oda:description "Input GTI filename for one processing interval"
attfile = "NONE"  # oda:POSIXPath; oda:description "Attitude file covering processing interval"
outattfile = "NONE"  # oda:POSIXPath; oda:description "Output filtered attitude filename"
outgtifile = "NONE"  # oda:POSIXPath; oda:description "Output pointing-filtered GTI filename (or NONE)"
med_ra = (
    -999.0
)  # oda:Float,oda:optional; oda:description "Median R.A. [deg] (set upon exit)"
med_dec = (
    -999.0
)  # oda:Float,oda:optional; oda:description "Median Dec. [deg] (set upon exit)"
med_roll = (
    -999.0
)  # oda:Float,oda:optional; oda:description "Median Roll [deg] (set upon exit)"
expotot = (
    -1.0
)  # oda:Float,oda:optional; oda:description "Total image exposure [sec] (set upon exit)"
expobad = (
    -1.0
)  # oda:Float,oda:optional; oda:description "Bad-attitude duration [sec] (set upon exit)"
point_toler = 0.025  # oda:Float,oda:optional ; oda:lower_limit 0; oda:description "Pointing error tolerance [deg]"
roll_toler = 0.083  # oda:Float,oda:optional ; oda:lower_limit 0; oda:description "Roll error tolerance [deg]"
alignfile = "CALDB"  # oda:POSIXPath,oda:optional; oda:description "Spacecraft-level alignment file (or CALDB)"
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
gtifile = str(inp_pdic["gtifile"])
attfile = str(inp_pdic["attfile"])
outattfile = str(inp_pdic["outattfile"])
outgtifile = str(inp_pdic["outgtifile"])

med_ra = (
    float(inp_pdic["med_ra"])
    if inp_pdic.get("med_ra", None) is not None
    else None
)

med_dec = (
    float(inp_pdic["med_dec"])
    if inp_pdic.get("med_dec", None) is not None
    else None
)

med_roll = (
    float(inp_pdic["med_roll"])
    if inp_pdic.get("med_roll", None) is not None
    else None
)

expotot = (
    float(inp_pdic["expotot"])
    if inp_pdic.get("expotot", None) is not None
    else None
)

expobad = (
    float(inp_pdic["expobad"])
    if inp_pdic.get("expobad", None) is not None
    else None
)

point_toler = (
    float(inp_pdic["point_toler"])
    if inp_pdic.get("point_toler", None) is not None
    else None
)

roll_toler = (
    float(inp_pdic["roll_toler"])
    if inp_pdic.get("roll_toler", None) is not None
    else None
)

alignfile = (
    str(inp_pdic["alignfile"])
    if inp_pdic.get("alignfile", None) is not None
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

command = "batsurvey-aspect"
params_names = "infile,gtifile,attfile,outattfile,outgtifile,med_ra,med_dec,med_roll,expotot,expobad,point_toler,roll_toler,alignfile,clobber,chatter,history,mode"
optional_file_params = "alignfile"

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
    ("out_batsurvey_aspect_result", "result_galaxy.output", result)
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
