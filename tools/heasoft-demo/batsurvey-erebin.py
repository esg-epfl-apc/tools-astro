#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python

# This script is generated with nb2galaxy

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

infile = "NONE"  # oda:POSIXPath; oda:description "Input survey file name"
calfile = "NONE"  # oda:POSIXPath; oda:description "Name gain/offset file or directory"
residfile = "CALDB"  # oda:POSIXPath,oda:optional; oda:description "BAT energy residual file (or CALDB)"
pulserfile = "CALDB"  # oda:POSIXPath,oda:optional; oda:description "Pulser DAC to energy calibration file name (or CALDB)"
fltpulserfile = "CALDB"  # oda:POSIXPath,oda:optional; oda:description "As-flown pulser DAC to energy calibration file name (or CALDB)"
outmap = "NONE"  # oda:String,oda:optional; oda:description "Name of output detector quality map (or NONE)"
baterebin_opts = (
    ""  # oda:String,oda:optional; oda:description "Extra options to baterebin"
)
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
infile = str(inp_pdic["infile"])
calfile = str(inp_pdic["calfile"])

residfile = (
    str(inp_pdic["residfile"])
    if inp_pdic.get("residfile", None) is not None
    else None
)

pulserfile = (
    str(inp_pdic["pulserfile"])
    if inp_pdic.get("pulserfile", None) is not None
    else None
)

fltpulserfile = (
    str(inp_pdic["fltpulserfile"])
    if inp_pdic.get("fltpulserfile", None) is not None
    else None
)

outmap = (
    str(inp_pdic["outmap"])
    if inp_pdic.get("outmap", None) is not None
    else None
)

baterebin_opts = (
    str(inp_pdic["baterebin_opts"])
    if inp_pdic.get("baterebin_opts", None) is not None
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

command = "batsurvey-erebin"
params_names = "infile,infile,outfile,calfile,residfile,pulserfile,fltpulserfile,outmap,baterebin_opts,chatter,clobber,history,mode"
optional_file_params = "residfile,pulserfile,fltpulserfile"

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
    ("out_batsurvey_erebin_result", "result_galaxy.output", result)
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
