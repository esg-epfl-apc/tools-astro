#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python

# This script is generated with nb2galaxy

# flake8: noqa

import json
import os
import shutil

from oda_api.json import CustomJSONEncoder

infile = "NONE"  # oda:POSIXPath; oda:description "Input detector plane image file name"
detmask = "NONE"  # oda:POSIXPath,oda:optional; oda:description "Input quality mask file name (or NONE)"
row = (
    -1
)  # oda:Integer,oda:optional ; oda:lower_limit -1; oda:description "Input row/image number to process"
keepfract = 0.98  # oda:Float,oda:optional ; oda:lower_limit 0.01 ; oda:upper_limit 1.0; oda:description "Fraction of non-zero pixels to keep"
guardfract = 0.25  # oda:Float,oda:optional ; oda:lower_limit 0.0; oda:description "Fractional size of guard band (w.r.t. median)"
guardval = 5.0  # oda:Float,oda:optional ; oda:lower_limit 0.0; oda:description "Further additive guard band beyond guardfract"
applygaps = True  # oda:Boolean,oda:optional; oda:description "Ignore positions of BAT detector gaps?"
mergemask = True  # oda:Boolean,oda:optional; oda:description "Logically AND input detmask with result?"
zerothresh = 20.0  # oda:Float,oda:optional ; oda:lower_limit 0; oda:description "Mean counts threshold above which zero counts is considered COLD"
goodval = (
    0  # oda:Integer,oda:optional; oda:description "Mask value for good pixel"
)
badval = (
    1  # oda:Integer,oda:optional; oda:description "Mask value for bad pixel"
)
hotval = (
    1  # oda:Integer,oda:optional; oda:description "Mask value for hot pixel"
)
coldval = (
    1  # oda:Integer,oda:optional; oda:description "Mask value for cold pixel"
)
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

detmask = (
    str(inp_pdic["detmask"])
    if inp_pdic.get("detmask", None) is not None
    else None
)

row = int(inp_pdic["row"]) if inp_pdic.get("row", None) is not None else None

keepfract = (
    float(inp_pdic["keepfract"])
    if inp_pdic.get("keepfract", None) is not None
    else None
)

guardfract = (
    float(inp_pdic["guardfract"])
    if inp_pdic.get("guardfract", None) is not None
    else None
)

guardval = (
    float(inp_pdic["guardval"])
    if inp_pdic.get("guardval", None) is not None
    else None
)

applygaps = (
    bool(inp_pdic["applygaps"])
    if inp_pdic.get("applygaps", None) is not None
    else None
)

mergemask = (
    bool(inp_pdic["mergemask"])
    if inp_pdic.get("mergemask", None) is not None
    else None
)

zerothresh = (
    float(inp_pdic["zerothresh"])
    if inp_pdic.get("zerothresh", None) is not None
    else None
)

goodval = (
    int(inp_pdic["goodval"])
    if inp_pdic.get("goodval", None) is not None
    else None
)

badval = (
    int(inp_pdic["badval"])
    if inp_pdic.get("badval", None) is not None
    else None
)

hotval = (
    int(inp_pdic["hotval"])
    if inp_pdic.get("hotval", None) is not None
    else None
)

coldval = (
    int(inp_pdic["coldval"])
    if inp_pdic.get("coldval", None) is not None
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

command = "bathotpix"
params_names = "infile,infile,outfile,detmask,row,keepfract,guardfract,guardval,applygaps,mergemask,zerothresh,goodval,badval,hotval,coldval,clobber,chatter,history,mode"
optional_file_params = "detmask"

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
_oda_outs.append(("out_bathotpix_result", "result_galaxy.output", result))

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
