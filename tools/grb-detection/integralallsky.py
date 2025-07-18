#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python

# This script is generated with nb2galaxy

# flake8: noqa

import json
import os
import shutil
import time

import integralclient as ic
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from matplotlib import pylab
from scipy import stats

# %matplotlib notebook

RA = 293.732
Dec = 21.8967222
tstart_rel_mseconds = 300.0
tstop_rel_seconds = 300.0
t0_utc = "2023-06-27T01:53:37.000"
# t0_utc=Time(Time("2022-10-14T19:21:47").mjd - 8.632259375000002/24/3600, format='mjd').isot.replace(" ", "T") # hard x-ray
# t0_utc="2022-10-14T19:21:47"
# rt=1
# nrt=1
# arc=0
required_completeness = 0.3
# mode="rt" # scw|rt|arc
mode = "scw"  # scw|rt|arc
global_snr_threshold = 3.0
negative_excesses = 0

_galaxy_wd = os.getcwd()

with open("inputs.json", "r") as fd:
    inp_dic = json.load(fd)
if "C_data_product_" in inp_dic.keys():
    inp_pdic = inp_dic["C_data_product_"]
else:
    inp_pdic = inp_dic
RA = float(inp_pdic["RA"])
Dec = float(inp_pdic["Dec"])
tstart_rel_mseconds = float(inp_pdic["tstart_rel_mseconds"])
tstop_rel_seconds = float(inp_pdic["tstop_rel_seconds"])
t0_utc = str(inp_pdic["t0_utc"])
required_completeness = float(inp_pdic["required_completeness"])
mode = str(inp_pdic["mode"])
global_snr_threshold = float(inp_pdic["global_snr_threshold"])
negative_excesses = int(inp_pdic["negative_excesses"])

t0_utc

if mode == "scw":
    rt = 0
    nrt = 1
    arc = 0
elif mode == "rt":
    rt = 1
    nrt = 0
    arc = 0
elif mode == "arc":
    rt = 0
    nrt = 0
    arc = 1
elif mode == "flags":
    print("mode set by flags")
else:
    raise Exception("unknown mode: {}, allowed: scw, rt".format(mode))

source_coord = SkyCoord(RA, Dec, unit="deg")

import importlib

import integralenv

# /home/savchenk/work/transients/workflows/integral-all-sky

importlib.reload(integralenv)

arc_root_prefix = integralenv.get_arc_root_prefix()

now_ijd = float(
    ic.converttime("UTC", time.strftime("%Y-%m-%dT%H:%M:%S"), "IJD")
)
t0_ijd = float(ic.converttime("UTC", t0_utc, "IJD"))

tstart_ijd = t0_ijd - tstart_rel_mseconds / 24.0 / 3600
tstop_ijd = t0_ijd + tstop_rel_seconds / 24.0 / 3600

now_ijd, t0_ijd, tstart_ijd, tstop_ijd

lcs = {}

def get_rtlc(t0_utc, tstop_rel_seconds):
    url = f"https://www.astro.unige.ch/mmoda/dispatch-data/gw/integralhk/api/v1.0/rtlc/{t0_utc}/{tstop_rel_seconds}?json&prophecy"
    print("url:", url)
    r = requests.get(url)
    return np.array(r.json()["lc"]["data"])

if nrt == 1:
    import isdcclient

    IC = isdcclient.ISDCClient()

    lcs["ACS"] = IC.genlc(
        "ACS",
        t0_utc,
        "%.10lg" % max(tstart_rel_mseconds, tstop_rel_seconds),
        format="numpy",
    )

    print("got lc from service of", lcs["ACS"])

    lcs["ACS"][:, 1] = 0.05
    print("got ACS", lcs["ACS"])

    # lcs['IBIS/Veto'] = IC.genlc("IBIS_VETO", t0_utc, "%.10lg"%max(tstart_rel_mseconds,tstop_rel_seconds),format='numpy')
    # lcs['IBIS/Veto'][:,1] = 8.

import requests

if rt == 1:
    lc = get_rtlc(t0_utc, tstop_rel_seconds)
    lcs["ACS"] = np.vstack([lc[:, 1], lc[:, 1] * 0 + 0.05, lc[:, 0]]).T

import isdcclient

IC = isdcclient.ISDCClient()

from astropy.time import Time

try:
    t_ref = Time(Time(t0_utc).mjd - 6 / 24, format="mjd").isot

    scwlc_ref = IC.genlc("ACS", t_ref, 10, format="numpy")
    rtlc_ref = get_rtlc(t_ref, 10)

    from matplotlib import pyplot as plt

    plt.figure(figsize=(10, 5))
    plt.plot(rtlc_ref[:, 1], rtlc_ref[:, 0], label="ACS RT")
    plt.plot(scwlc_ref[:, 0], scwlc_ref[:, 2], label="ACS")
    # TODO: for extra timing accuracy, extract here
except Exception as e:
    print("failed to get ref lc", e)

from matplotlib import pylab as plt

plt.figure()

plt.step(
    (lcs["ACS"][:, 0] - t0_ijd) * 24 * 3600,
    lcs["ACS"][:, 2],
)

plt.xlim(-5, 5)

def rebin(lc, n, av=False):
    if n == 0:
        return lc

    N = int(lc.shape[0] / n) * n
    if av:
        return lc[:N].reshape((int(lc.shape[0] / n), n)).mean(1)
    else:
        return lc[:N].reshape((int(lc.shape[0] / n), n)).sum(1)

import time

# if rt == 1:

#     got_data = False

#     while not got_data:
#         current_rev=float(ic.converttime("UTC",t0_utc,"REVNUM"))

#         print("current rev", current_rev)

#         rtdata_roots=[
#             '/unsaved/astro/savchenk/dockers/realtimeacs/docker-ibas/spiacs-lcdump',
#             '/rtdata',
#             '/mnt/sshfs/isdc-in01//unsaved/astro/savchenk/dockers/realtimeacs/docker-ibas/spiacs-lcdump',
#         ]

#         for realtime_dump_root in rtdata_roots + [ None ]:
#             #print("probing",realtime_dump_root,"with",glob.glob(realtime_dump_root+"/lcdump-revol-*.csv"))
#             if realtime_dump_root and len(glob.glob(realtime_dump_root+"/lcdump-revol-*.csv"))>0:
#                 print("this",realtime_dump_root)
#                 break

#         if not realtime_dump_root:
#             raise Exception("no realtime archvie found")

#         for rt_fn in reversed(sorted([l for l in glob.glob(realtime_dump_root+"/lcdump-revol-*.csv") if
#                        float(re.search("lcdump-revol-(\d{4}).*.csv",l).groups()[0])<=current_rev+1])):

#             print(rt_fn)

#             rt_lc = np.genfromtxt(rt_fn)

#             lcs['ACS']=rt_lc[:,(3,0,2,0)]
#             lcs['ACS'][:,1] = 0.05

#             first_data = lcs['ACS'][:,0][0]
#             last_data = lcs['ACS'][:,0][-1]

#             print("now", now_ijd,
#                   "first data in file", first_data,
#                   "last data", last_data,
#                   "requested", t0_ijd,
#                   "have margin", (last_data-t0_ijd)*24*3600,"s",
#                   "data delay", (now_ijd-last_data)*24*3600,"s")

#             if t0_ijd<first_data:
#                 print("data in the previous file")
#                 continue

#             print("margin",(last_data-now_ijd)*24*3600-tstop_rel_seconds*1.5 + 100)
#             if  (last_data-t0_ijd)*24*3600>tstop_rel_seconds*1.5 + 100:
#                 print("this margin is sufficient")
#                 got_data=True
#                 break
#             else:
#                 print("this margin is NOT sufficient, waiting")
#             #    if (now_ijd-last_data)*24*3600>1000:
#             #        raise RuntimeError('margin insufficent, data too old: no more hope')

#                 time.sleep(30)
#                 break

# lcs['ACS']

lcs["ACS"][:, 0].min()

summary = dict()

for n, lc in lcs.items():

    try:
        rel_s = (lc[:, 0] - t0_ijd) * 24 * 3600
    except:
        continue

    m = rel_s > -tstart_rel_mseconds
    m &= rel_s < tstop_rel_seconds

    print("total lc", lc.shape)
    print("min", lc[:, 0].min() - t0_ijd)
    print("max", lc[:, 0].max() - t0_ijd)

    lc = lc[m]

    b_tb = np.mean(lc[:, 1])

    rel_s = (lc[:, 0] - t0_ijd) * 24 * 3600

    expected_telapse = tstop_rel_seconds + tstop_rel_seconds

    if len(rel_s) == 0:
        telapse = 0
        ontime = 0
    else:
        telapse = rel_s.max() - rel_s.min()
        ontime = np.sum(lc[:, 1])

    print(
        "expected telapse",
        expected_telapse,
        "telapse",
        telapse,
        "ontime",
        ontime,
    )

    if float(ontime) / expected_telapse < required_completeness:
        raise Exception(
            "data not available: exected %.5lg elapsed %.5lg ontime %.5lg completeness %s requireed %s"
            % (
                expected_telapse,
                telapse,
                ontime,
                ontime / expected_telapse,
                required_completeness,
            )
        )

    lc_summary = dict()
    summary[n.replace("/", "_")] = lc_summary

    print("size", lc.shape, rel_s.shape)

    if np.sum(m) == 0:
        continue

    pylab.figure(figsize=(8, 6))

    for ascale in [0.05, 0.5, 1, 10]:
        summary_scale = dict()
        lc_summary[("s_%.5lg" % ascale).replace(".", "_")] = summary_scale

        print("requested scale", ascale)
        print("b_tb", b_tb)

        if b_tb > ascale:
            ascale = b_tb

        nscale = int(ascale / b_tb)
        scale = nscale * b_tb

        print("acceptable, will be", nscale, scale)

        rate = rebin(lc[:, 2], nscale, False) / scale
        rate_err = rebin(lc[:, 2], nscale, False) ** 0.5 / scale

        print("rebinned to", rate.shape)

        pylab.errorbar(
            rebin(rel_s, nscale, True), rate, rate_err, xerr=scale / 4.0
        )

        summary_scale["meanrate"] = np.mean(rate)
        summary_scale["maxrate"] = np.max(rate)
        summary_scale["stdvar"] = np.std(rate)
        summary_scale["meanerr"] = np.mean(rate_err**2) ** 0.5
        summary_scale["excvar"] = (
            summary_scale["stdvar"] / summary_scale["meanerr"]
        )

        summary_scale["maxsnr"] = np.max(
            (rate - np.mean(rate)) / rate_err / summary_scale["excvar"]
        )

        summary_scale["localfar"] = (
            stats.norm.sf(summary_scale["maxsnr"]) * rate.shape[0]
        )

        summary_scale["localfar_s"] = (
            stats.norm.isf(summary_scale["localfar"] / 2.0)
            if summary_scale["localfar"] < 1
            else 0
        )

        # add FAR spike here

        if (
            "best" not in lc_summary
            or summary_scale["localfar_s"] > lc_summary["best"]["localfar_s"]
        ):
            lc_summary["best"] = dict(
                localfar_s=summary_scale["localfar_s"],
                scale=ascale,
            )

        print(summary_scale)

    # tight_layout()
    pylab.grid()

    pylab.xlim(-tstart_rel_mseconds, tstop_rel_seconds)
    # pylab.axhspan(0,10,alpha=0.2,color="red")
    # pylab.axhspan(10,15,alpha=0.2,color="green")
    # pylab.axhspan(15,20,alpha=0.2,color="blue")
    pylab.ylabel(n + ", count s$^{-1}$")
    # ylim([0,50])
    pylab.xlabel("seconds since %s (IJD %.10lg)" % (t0_utc, t0_ijd))

    fn = n.replace("/", "_") + "_lc.png"
    pylab.savefig(fn)
    print("saving as", fn)
    break

# below S/N of 4 FAR is determined primarily by poisson, above - by spikes

def approx_FAR_spike_hz(snr, scale):
    lim_snr = 2

    spike_rate_snr6 = 60.0 / 3600.0 / 24.0
    if scale >= 0.1:
        spike_rate_snr6 *= (scale / 0.1) ** -1

    approx_FAR_hz = snr * 0 + spike_rate_snr6 * (lim_snr / 6.0) ** -2.7

    try:
        if snr > lim_snr:
            approx_FAR_hz = spike_rate_snr6 * (np.abs(snr) / 6.0) ** -2.7
    except:
        m = snr > lim_snr
        approx_FAR_hz[m] = (np.abs(snr[m]) / 6.0) ** -2.7 * spike_rate_snr6

    return approx_FAR_hz

def approx_FAR_norm_hz(snr, scale_s):
    return stats.norm.sf(snr) / scale_s

def approx_FAP(snr, t, scale_s):

    try:
        t_scaled = t[:]
        t_scaled[abs(t) < scale_s] = scale_s
    except:
        if abs(t) < scale_s:
            t_scaled = scale_s
        else:
            t_scaled = t

    approx_FAP = (
        2
        * (
            approx_FAR_norm_hz(snr, scale_s)
            + approx_FAR_spike_hz(snr, scale_s)
        )
        * abs(t_scaled)
        * (1 + np.log(30 / 0.1))
    )

    return approx_FAP

pylab.figure()

x = np.linspace(-5, 10, 100)

for scale_s in 0.05, 0.1, 1, 10:

    c = pylab.plot(x, approx_FAR_norm_hz(x, scale_s), ls="--")
    pylab.plot(x, approx_FAR_spike_hz(x, scale_s), c=c[0].get_color(), ls=":")
    pylab.plot(
        x,
        approx_FAR_spike_hz(x, scale_s) + approx_FAR_norm_hz(x, scale_s),
        c=c[0].get_color(),
    )

    pylab.semilogy()

pylab.ylim([1e-5, 30])

timescales = sorted(
    set(
        [
            0.05 * ns
            for ns in sorted(
                set(list(map(int, np.logspace(0, np.log10(20 * 30), 100))))
            )
        ]
        + list(np.linspace(1, 31, 30 * 2 + 1))
    )
)
timescales

summary = dict()
all_excesses = []

best_lc = None

for n, lc in lcs.items():

    # rel_s = lc[:,0]
    rel_s = (lc[:, 0] - t0_ijd) * 24 * 3600

    m = rel_s > -tstart_rel_mseconds
    m &= rel_s < tstop_rel_seconds

    print("total lc", lc.shape)
    print("min", lc[:, 0].min() - t0_ijd)
    print("max", lc[:, 0].max() - t0_ijd)

    lc = lc[m]
    # rel_s = lc[:,0]

    b_tb = np.mean(lc[:, 1])

    rel_s = (lc[:, 0] - t0_ijd) * 24 * 3600

    expected_telapse = tstop_rel_seconds + tstop_rel_seconds

    if len(rel_s) == 0:
        telapse = 0
        ontime = 0
    else:
        telapse = rel_s.max() - rel_s.min()
        ontime = np.sum(lc[:, 1])

    print(
        "expected telapse",
        expected_telapse,
        "telapse",
        telapse,
        "ontime",
        ontime,
    )

    if ontime / expected_telapse < required_completeness:
        raise Exception(
            "data not available: exected %.5lg elapsed %.5lg ontime %.5lg"
            % (expected_telapse, telapse, ontime)
        )

    lc_summary = dict()
    summary[n.replace("/", "_")] = lc_summary

    print("size", lc.shape, rel_s.shape)

    if np.sum(m) == 0:
        continue

    pylab.figure(figsize=(8, 6))

    best_lc_byscale = {}

    # for ascale in [0.05, 0.1, 0.2, 0.5, 1, 2, 10]:
    for ascale in timescales:
        # for ascale in [0.05, 0.1, 0.5]:
        # for ascale in [0.05*i for i in range(20)] + [0.5*i for i in range(20)] + [15, 20, 25, 30]:
        # for ascale in [0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 1, 2, 8, 10]:
        # for ascale in [1,]:
        s_scale_mo = {}
        lc_summary[("s_%.5lg" % ascale).replace(".", "_")] = s_scale_mo

        print("requested scale", ascale)
        #        print("b_tb",b_tb)

        if b_tb > ascale:
            ascale = b_tb

        nscale = int(round(ascale / b_tb))
        scale = nscale * b_tb

        print("true scale", scale)

        #        print("acceptable, will be", nscale, scale)

        c = None

        # for offset in range(0,nscale):
        # for offset in (, round(nscale/2)):

        if nscale < 20:
            offsets = range(0, round(nscale / 2) + 1)
        else:
            offsets = range(
                0, round(nscale / 2) + 1, max(round(round(nscale / 2) / 20), 1)
            )

        for offset in offsets:
            summary_scale = dict()
            s_scale_mo[offset] = summary_scale

            rel_s_scale = rebin(rel_s[offset:], nscale, True)
            rate = rebin(lc[offset:, 2], nscale, False) / scale
            rate_err = rebin(lc[offset:, 2], nscale, False) ** 0.5 / scale

            # print("rebinned to",rate.shape)
            print("offset", offset, "rebinned to", rate.shape)

            summary_scale["scale_s"] = scale
            summary_scale["meanrate"] = np.mean(rate)
            summary_scale["maxrate"] = np.max(rate)
            summary_scale["stdvar"] = np.std(rate)
            summary_scale["meanerr"] = np.mean(rate_err**2) ** 0.5
            summary_scale["excvar"] = (
                summary_scale["stdvar"] / summary_scale["meanerr"]
            )

            print("summary_scale['excvar']", summary_scale["excvar"])

            if negative_excesses == 1:
                snr = (
                    -(rate - np.mean(rate))
                    / rate_err
                    / summary_scale["excvar"]
                )
            else:
                snr = (
                    (rate - np.mean(rate)) / rate_err / summary_scale["excvar"]
                )

            i_max = np.argmax(snr)

            print(i_max, snr[i_max], rel_s_scale[i_max])

            summary_scale["maxsnr"] = snr[i_max]
            summary_scale["maxsnr_t"] = rel_s_scale[i_max]

            summary_scale["localfar"] = (
                stats.norm.sf(summary_scale["maxsnr"]) * rate.shape[0]
            )

            summary_scale["localfar_s"] = (
                stats.norm.isf(summary_scale["localfar"] / 2.0)
                if summary_scale["localfar"] < 1
                else 0
            )

            m_over_threshold = snr > global_snr_threshold

            excesses = dict(
                snr=snr[m_over_threshold],
                rel_s_scale=rel_s_scale[m_over_threshold],
                rate=rate[m_over_threshold],
                rate_err=rate_err[m_over_threshold],
                rate_overbkg=rate[m_over_threshold] - np.mean(rate),
            )

            summary_scale["excesses"] = [
                dict(zip(excesses.keys(), er))
                for er in zip(*excesses.values())
            ]

            for e in summary_scale["excesses"]:
                e["FAP"] = approx_FAP(e["snr"], e["rel_s_scale"], scale)

            all_excesses += [
                dict(scale=scale, offset=offset, excess=e)
                for e in summary_scale["excesses"]
            ]

            print(
                "scale",
                scale,
                "offset",
                offset,
                "found excesses",
                len(summary_scale["excesses"]),
            )

            # r=pylab.errorbar(
            #    rebin(rel_s[offset:],nscale,True),
            #    rate,
            #    rate_err,
            #    xerr=scale/4.,
            #    c=c,
            #    alpha=0.7
            # )

            #    print(rel_s_scale.shape, snr.shape)

            r = pylab.errorbar(
                rel_s_scale, snr, snr * 0 + 1, xerr=scale / 4.0, c=c, alpha=0.7
            )

            pylab.axvline(summary_scale["maxsnr_t"], c="k")

            c = r[0].get_color()

            # add FAR spike here

            if (
                "best" not in lc_summary
                or summary_scale["localfar_s"]
                > lc_summary["best"]["localfar_s"]
            ):
                lc_summary["best"] = dict(
                    localfar_s=summary_scale["localfar_s"],
                    scale=ascale,
                    summary_scale=summary_scale,
                )
                best_lc = rel_s_scale, rate, rate_err

            if (
                "best" not in s_scale_mo
                or summary_scale["localfar_s"]
                > s_scale_mo["best"]["localfar_s"]
            ):
                s_scale_mo["best"] = dict(
                    localfar_s=summary_scale["localfar_s"],
                    scale=ascale,
                    summary_scale=summary_scale,
                )
            #  best_lc=rel_s_scale,rate,rate_err

            if (
                ascale not in best_lc_byscale
                or summary_scale["localfar_s"]
                > best_lc_byscale[ascale]["localfar_s"]
            ):
                best_lc_byscale[ascale] = dict(
                    localfar_s=summary_scale["localfar_s"],
                    scale=ascale,
                    summary_scale=summary_scale,
                    best_lc=(rel_s_scale, rate, rate_err),
                )

            # print(summary_scale)
        s_scale_mo.update(s_scale_mo["best"]["summary_scale"])

    # tight_layout()
    pylab.grid()

    # pylab.xlim(-tstart_rel_mseconds, tstop_rel_seconds)
    # pylab.axhspan(0,10,alpha=0.2,color="red")
    # pylab.axhspan(10,15,alpha=0.2,color="green")
    # pylab.axhspan(15,20,alpha=0.2,color="blue")
    pylab.ylabel(n + ", S/N")
    # ylim([0,50])
    pylab.xlabel("seconds since %s (IJD %.10lg)" % (t0_utc, t0_ijd))
    pylab.xlim([-10, 10])

    detfn = n.replace("/", "_") + "_det_lc.png"
    pylab.savefig(detfn)
    print("saving as", detfn)

summary["ACS"]["best"]

plt.figure()

snr = (rate - np.mean(rate)) / rate_err / summary_scale["excvar"]

plt.step(rel_s_scale, rate)
plt.xlim(-200, 200)

for e in all_excesses:
    if np.abs(e["scale"] - 0.05) < 0.001:
        # if np.abs(e['excess']['rel_s_scale']) < 5:
        print(e["excess"]["rel_s_scale"], e)

grouped_excesses = []

for i in sorted(all_excesses, key=lambda x: x["excess"]["FAP"]):
    if i["excess"]["FAP"] < 1 or True:
        print(
            i["scale"],
            i["offset"],
            i["excess"]["snr"],
            i["excess"]["rel_s_scale"],
            i["excess"]["FAP"],
        )

        grouped = False
        for g in grouped_excesses:
            if abs(
                i["excess"]["rel_s_scale"] - g["excess"]["rel_s_scale"]
            ) < max(i["scale"], g["scale"]):
                print("to group", g["excess"]["rel_s_scale"])
                if i["excess"]["snr"] > g["excess"]["snr"]:
                    print("group takeover")
                    g.update(i)
                grouped = True

        if not grouped:
            print("new group")
            # i['group']=[i]
            grouped_excesses.append(i)

grouped_excesses = sorted(grouped_excesses, key=lambda x: x["excess"]["FAP"])

for i in grouped_excesses:
    print(
        f"timescale {i['scale']:4.2f}   S/N {i['excess']['snr']:5.2f}   T0+{i['excess']['rel_s_scale']:7.1f}   FAP {i['excess']['FAP']:7.5f}"
    )

import json

len(json.dumps(grouped_excesses))

summary["ACS"]["best"]

# T

summary["ACS"]["s_8"]

excvar_summary = dict()

for k, s in summary["ACS"].items():
    if "scale_s" in s:
        print("%.5lg" % s["scale_s"], "%5.4lg" % s["excvar"])

        if s["scale_s"] <= 0.200:
            kg = "hf_200ms"
        elif s["scale_s"] <= 2.00:
            kg = "mf_200ms_2s"
        elif s["scale_s"] <= 10.00:
            kg = "mf_2s_10s"
        else:
            kg = "lf_10s"

        if kg not in excvar_summary:
            excvar_summary[kg] = [s["excvar"]]
        else:
            excvar_summary[kg] += [s["excvar"]]

for k, v in excvar_summary.items():
    print(k, min(v), max(v))

fig_names = []

for limit_group in 0.02, 0.1, 1:
    figs = dict()

    for n, lc in lcs.items():
        rel_s = (lc[:, 0] - t0_ijd) * 24 * 3600

        m = rel_s > -tstart_rel_mseconds
        m &= rel_s < tstop_rel_seconds

        print("total lc", lc.shape)
        print("min", lc[:, 0].min() - t0_ijd)
        print("max", lc[:, 0].max() - t0_ijd)

        lc = lc[m]
        rel_s = rel_s[m]

        for excess in grouped_excesses:
            # if excess['excess']['FAP'] > 0.02: continue
            if excess["excess"]["FAP"] > limit_group:
                continue

            print(excess)

            offset = excess["offset"]
            nscale = int(excess["scale"] / b_tb)
            scale = excess["scale"]

            s_figs = sorted(figs.items(), key=lambda x: abs(x[0] - scale))

            if (
                len(s_figs) == 0
                or s_figs[0][0] < scale * 0.5
                or s_figs[0][0] > scale * 1.5
            ):
                fig = pylab.figure(figsize=(8, 6))
                figs[scale] = fig
                pylab.xlim([-2, 2])
                pylab.xlabel("seconds since " + t0_utc)
                pylab.ylabel("counts/s")
                pylab.title("FAP threshold %.5lg" % limit_group)
            else:
                print("good match", s_figs[0][0], scale)
                pylab.figure(s_figs[0][1].number)
                pylab.xlabel("seconds since " + t0_utc)
                pylab.ylabel("counts/s")

            rel_s_scale = rebin(rel_s[offset:], nscale, True)
            rate = rebin(lc[offset:, 2], nscale, False) / scale
            rate_err = rebin(lc[offset:, 2], nscale, False) ** 0.5 / scale

            bkg = np.mean(rate)

            m_on = (
                np.abs(rel_s_scale - excess["excess"]["rel_s_scale"])
                < excess["scale"] * 1.5
            )

            pylab.grid(False)

            pylab.axhline(0, alpha=0.2, ls=":", color="gray")

            cr = pylab.errorbar(
                rel_s_scale,
                (rate - bkg),
                (rate_err),
                alpha=0.5,
                ls="",
            )[0].get_color()

            pylab.step(
                rel_s_scale,
                (rate - bkg),
                #      (rate_err),
                alpha=0.5,
                where="mid",
                c=cr,
            )

            pylab.axhline(np.std(rate) * 3, alpha=0.2, ls="--", c=cr)
            pylab.axhline(np.std(rate) * 5, alpha=0.2, ls="--", lw=2, c=cr)

            pylab.errorbar(
                rel_s_scale[m_on],
                (rate - bkg)[m_on],
                (rate_err)[m_on],
                lw=2.0,
                alpha=1,
                label="S/N %.3lg FAP %.3lg scale %.3lg s"
                % (
                    excess["excess"]["snr"],
                    excess["excess"]["FAP"],
                    excess["scale"],
                ),
                c=cr,
            )

            newlim = [
                min(
                    [
                        excess["excess"]["rel_s_scale"] * 1.3
                        - excess["scale"] * 5,
                        -excess["scale"] * 5,
                    ]
                ),
                max(
                    [
                        excess["excess"]["rel_s_scale"] * 1.3
                        + excess["scale"] * 5,
                        excess["scale"] * 5,
                    ]
                ),
            ]

            oldlim = pylab.gca().get_xlim()

            print(oldlim)

            pylab.xlim(
                [
                    min([oldlim[0], newlim[0]]),
                    max([oldlim[1], newlim[1]]),
                ]
            )

    for f_i, (s, f) in enumerate(figs.items()):
        f.legend()
        f.gca().axvline(0, ls="--", c="r", lw=3)
        fn = "excess_%.5lg_%i.png" % (s, len(fig_names))
        f.savefig(fn)
        fig_names.append(fn)

cols = 1
rows = int(np.ceil(len(fig_names) / cols))

if rows > 0:
    f, axes = pylab.subplots(rows, cols, figsize=(12, 8 * rows))
    print("axes", axes, axes.__class__)

    if rows > 1:
        axes = axes.flatten()
    else:
        axes = [axes]

    for i, fn in enumerate(fig_names):
        # f.add_subplot(len(fig_names), 2, i+1)
        axes[i].axis("off")
        axes[i].imshow(pylab.imread(fn))  # , extent=(0,1,0,1))
        # pylab.imshow(pylab.imread(fn), extent=(0,1,(i-1)/len(fig_names),i/len(fig_names)))

    f.tight_layout()
else:
    f = pylab.figure()

f.savefig("excesses_mosaic.png")

if rt == 1:
    summary["ACS_rt"] = summary["ACS"]

summary["ACS"]["s_1"]["meanerr"]

import json

json.dump(
    dict(
        summary=summary,
        reportable_excesses=grouped_excesses,
        excvar_summary=excvar_summary,
    ),
    open("integral_all_sky.json", "w"),
    indent=4,
)

acs_lc_png = "ACS_lc.png"
acs_rt_lc_png = "ACS_lc.png"
acs_rt_det_lc_png = "ACS_det_lc.png"
ibis_veto_lc_png = "IBIS_Veto_lc.png"
excesses_mosaic_png = "excesses_mosaic.png"
summary = summary
reportable_excesses = grouped_excesses
excvar_summary = excvar_summary

# output gathering
_galaxy_meta_data = {}
_simple_outs = []
_simple_outs.append(
    ("out_integralallsky_acs_lc_png", "acs_lc_png_galaxy.output", acs_lc_png)
)
_simple_outs.append(
    (
        "out_integralallsky_acs_rt_lc_png",
        "acs_rt_lc_png_galaxy.output",
        acs_rt_lc_png,
    )
)
_simple_outs.append(
    (
        "out_integralallsky_acs_rt_det_lc_png",
        "acs_rt_det_lc_png_galaxy.output",
        acs_rt_det_lc_png,
    )
)
_simple_outs.append(
    (
        "out_integralallsky_ibis_veto_lc_png",
        "ibis_veto_lc_png_galaxy.output",
        ibis_veto_lc_png,
    )
)
_simple_outs.append(
    (
        "out_integralallsky_excesses_mosaic_png",
        "excesses_mosaic_png_galaxy.output",
        excesses_mosaic_png,
    )
)
_simple_outs.append(
    ("out_integralallsky_summary", "summary_galaxy.output", summary)
)
_simple_outs.append(
    (
        "out_integralallsky_reportable_excesses",
        "reportable_excesses_galaxy.output",
        reportable_excesses,
    )
)
_simple_outs.append(
    (
        "out_integralallsky_excvar_summary",
        "excvar_summary_galaxy.output",
        excvar_summary,
    )
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
