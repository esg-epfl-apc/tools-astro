import numpy as np
import sys
import os

import matplotlib.pyplot as plt

Mpc_in_h = 2.8590868063e10


def weighted_quantile(
    values, quantiles, sample_weight=None, values_sorted=False, old_style=False
):
    """Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of initial array
    :param old_style: if True, will correct output to be consistent with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(
        quantiles <= 1
    ), "quantiles should be in [0, 1]"

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with np.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)


def get_distance_Mpc(mc_file):
    # reading the source distance in Mpc from the data file
    with open(mc_file) as f_mc_lines:
        for line in f_mc_lines:
            if line.startswith("#"):
                cols = list(line.split())
                idx = cols.index("T/Mpc")
                if idx >= 0:
                    return float(cols[idx + 2])
    raise ValueError('Unexpected mc file format: "T/Mpc" not found')


default_params = {
    "weight_col": 2,
    "E_src_col": 8,
    "delay_col": 6,
    "comoving_distance": None,
    "n_points": 100,
    "logscale": False,
    "suffix": "",
    "Emin": 1e6,
    "Emax": 1e20,
    "psf": 1.0,
    "cut_0": False,
    "rounding_error": 0.0007,  # 1./600,
    "min_n_particles": 100,
    "out_ext": "png",
    "max_t": 0.05,  # use negative for percentile
    "show": False,
    "add_alpha": 0.0,  # multiply weights be E_src^{-add_alpha}
    # 'frac_time': 2000/3600,  # calculate fraction of flux coming within this time in hours
    "verbose": 2,
    "format": ".2g",  # float format
}


def filter_data(data, **kwargs):
    params = default_params.copy()
    params.update(kwargs)

    Emax = params["Emax"]
    Emin = params["Emin"]
    data = data[data[:, 0] <= Emax]
    data = data[data[:, 0] >= Emin]
    psf = params["psf"]
    data = data[data[:, 2] <= psf]
    cut0 = params["cut_0"]
    if cut0:
        col = params["delay_col"] - 1
        data = data[data[:, col] != 0.0]
    return data


def get_counts(rotated_mc_file, **kwargs):
    params = default_params.copy()
    params.update(kwargs)
    data = np.loadtxt(rotated_mc_file)
    data_filtered = filter_data(data, **kwargs)
    verbose = params["verbose"]
    if verbose > 1:
        print(len(data_filtered), "of", len(data), "has passed the filter")

    weight_col = params["weight_col"] - 1
    col = params["delay_col"] - 1
    comoving_distance = params["comoving_distance"]
    if not comoving_distance:
        comoving_distance = get_distance_Mpc(rotated_mc_file)

    x_scale = comoving_distance * Mpc_in_h  # convert to hours

    delay = data_filtered[:, col] * x_scale

    idxs = np.argsort(delay)
    delay = delay[idxs]

    if weight_col >= 0:
        assert weight_col < data_filtered.shape[1]
        weights = data_filtered[idxs, weight_col]
    else:
        weights = np.ones(len(idxs))

    add_alpha = params["add_alpha"]
    if add_alpha != 0:
        E_src = data_filtered[:, params["E_src_col"]]
        av_Esrc = np.exp(np.log(E_src).mean())
        weights *= np.power(E_src / av_Esrc, -add_alpha)

    return delay, weights


def light_curve(delay, weights, **kwargs):
    params = default_params.copy()
    params.update(kwargs)
    min_n_particles = params["min_n_particles"]
    min_bin_size = params["rounding_error"]
    max_t = params["max_t"]

    if max_t < 0:
        max_t = weighted_quantile(delay, [-0.01 * max_t], sample_weight=weights)[0]

    f = []
    t = []
    N = []

    bin_idx = 0
    if delay[0] < min_bin_size:
        bin_idx = np.where(delay < min_bin_size)[0][-1]
        if bin_idx + 1 < min_n_particles:
            bin_idx = min_n_particles
        wsum = np.sum(weights[:bin_idx])
        _t = np.sum(delay[:bin_idx] * weights[:bin_idx]) / wsum
        _t = max(_t, 0.5 * min_bin_size)
        t.append(_t)
        bin_size = max(delay[bin_idx] - delay[0], min_bin_size)
        f.append(wsum / bin_size)
        N.append(bin_idx)

    while True:
        xmin = (
            0.5 * (delay[bin_idx - 1] + delay[bin_idx])
            if bin_idx > 0
            else delay[bin_idx]
        )
        if xmin > max_t:
            break
        bin_idx2 = np.where(delay[bin_idx + min_n_particles :] > xmin + min_bin_size)[0]
        if len(bin_idx2) == 0:
            break
        bin_idx2 = bin_idx2[0] + bin_idx + min_n_particles
        _delay = delay[bin_idx:bin_idx2]
        _weights = weights[bin_idx:bin_idx2]
        wsum = _weights.sum()
        t.append(np.sum(_delay * _weights) / wsum)
        xmax = 0.5 * (delay[bin_idx2 - 1] + delay[bin_idx2])
        f.append(wsum / (xmax - xmin))
        N.append(bin_idx2 - bin_idx)
        bin_idx = bin_idx2

    return [np.array(x) for x in (t, f, N)]


def make_plot(path, label="verbose", **kwargs):
    params = default_params.copy()
    params.update(kwargs)
    path = os.path.expanduser(path)
    fmt = params["format"]
    col = params["delay_col"] - 1
    logscale = params["logscale"]
    verbose = params["verbose"]
    Emax = params["Emax"]
    Emin = params["Emin"]
    psf = params["psf"]
    cut0 = params["cut_0"]

    delay, weights = get_counts(path, **kwargs)
    t, f, _ = light_curve(delay, weights, **kwargs)

    suffix = params["suffix"]

    if cut0:
        suffix += "_cut0"
    if logscale:
        suffix += "_log"

    max_t = params["max_t"]
    out_name = (
        f"{path}_f{col + 1}_E{Emin:{fmt}}-{Emax:{fmt}}TeV_th{psf}_r{max_t}{suffix}."
    )

    x_label = "t, [h]"
    y_label = "dN/dt [a.u.]"

    if verbose > 1:
        for key, val in params.items():
            print(f"\t{key} = {val}")

    if max_t < 0:
        median, max_t = weighted_quantile(
            delay, [0.5, -0.01 * max_t], sample_weight=weights
        )
    else:
        median = weighted_quantile(delay, [0.5], sample_weight=weights)[0]

    # the histogram of the data
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    if label == "verbose":
        label = f"50% with delay<{median:{fmt}} h"
    plot_args = {}
    if label:
        plot_args["label"] = label

    plt.plot(t, f, "g-", **plot_args)

    if logscale:
        plt.xscale("log")
        plt.yscale("log")

    if label:
        plt.legend(loc="upper right")

    file_name = out_name + params["out_ext"].strip()
    plt.savefig(file_name)
    if verbose > 0:
        print("saved to", file_name)
    if params["show"]:
        plt.show()

    return file_name


if __name__ == "__main__":

    def usage_exit(reason=""):
        print(reason, file=sys.stderr)
        print(
            "usage: python",
            sys.argv[0],
            "data_file",
            *["{}={}".format(k, v) for k, v in default_params.items()],
            file=sys.stderr,
        )
        exit(1)

    if len(sys.argv) < 2:
        usage_exit()

    kwargs = {}
    for par in sys.argv[2:]:
        kv = par.split("=")
        if len(kv) != 2:
            usage_exit("invalid argument " + par)
        if kv[0] not in default_params:
            usage_exit("unknown parameter " + kv[0])
        constr = type(default_params[kv[0]])
        value = kv[1]
        if constr == bool:
            value = value.lower() == "false" or value == "0"
            value = not value
        else:
            value = constr(value)
        kwargs[kv[0]] = value

    make_plot(sys.argv[1], **kwargs)
