# TODO: There looks to be a bug in the median values because we are not rescaling the copy of the final data used for computing metrics
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

import rpy2.robjects as ro
import numpy as np
import pandas as pd
import statsmodels.robust

pandas2ri.activate()


def main(args):
    data = pd.read_csv(args.in_file, sep="\t")

    with pd.HDFStore(args.out_file, "w") as store:
        store["multipliers"] = pd.Series(args.multipliers)

        store["data"] = data

        if np.all(np.isnan(data[data["ideal"]]["cor_gc"])) or np.all(np.isnan(data[data["ideal"]]["copy"])):
            store["failed"] = pd.Series([])
            quit()

        for multiplier in args.multipliers:
            params = load_params(num_states=args.num_states)

            scaling_factor = compute_scaling_factor(data, multiplier, params)

            test = data.copy()

            test["multiplier"] = multiplier

            test["copy"] = scaling_factor * test["copy"]

            segmented = segment_data(test, params)

            modal_seg = segmented["segs"]

            modal_seg["multiplier"] = multiplier

            modal_seg["state"] = modal_seg["state"].astype(int) - 1

            test["state"] = segmented["state"]

            test["segment"] = (
                (test["chr"] != test["chr"].shift()) | (test["state"] != test["state"].shift())
            ).cumsum() - 1

            test["median"] = np.array(modal_seg["median"][x] for x in test["segment"])

            bins = test[
                [
                    "chr",
                    "start",
                    "end",
                    "reads",
                    "cor_gc",
                    "copy",
                    "ideal",
                    "median",
                    "state",
                    "segment",
                ]
            ]

            segs = modal_seg[["chr", "start", "end", "median", "state"]]

            metrics = get_metrics(bins, segs)

            metrics["log_likelihood"] = segmented["loglik"][-1]

            metrics["multiplier"] = multiplier

            metrics["true_multiplier"] = scaling_factor

            store["/{}/metrics".format(multiplier)] = pd.Series(metrics)

            store["/{}/params/init".format(multiplier)] = params

            for p in ["lambdas", "mus", "pi", "rho"]:
                store["/{m}/params/{p}".format(m=multiplier, p=p)] = pd.DataFrame(segmented[p]).T

            store["/{}/bins".format(multiplier)] = bins

            store["/{}/segs".format(multiplier)] = segs


def load_params(num_states=12):
    strength = 1000 * np.ones(num_states)
    e = 0.999999 * np.ones(num_states)
    mu = np.arange(0, num_states)
    l = 20 * np.ones(num_states)
    nu = 2.1 * np.ones(num_states)
    k = 25 * np.ones(num_states)
    k[0] = 100
    k[1] = 100
    k[2] = 700
    k[3] = 100
    m = np.arange(0, num_states)
    eta = 50000 * np.ones(num_states)
    g = 3 * np.ones(num_states)
    s = np.ones(num_states)
    df = np.column_stack([strength, e, mu, l, nu, k, m, eta, g, s])
    df = pd.DataFrame(
        df,
        columns=[
            "strength",
            "e",
            "mu",
            "lambda",
            "nu",
            "kappa",
            "m",
            "eta",
            "gamma",
            "S",
        ],
    )
    return df


def compute_scaling_factor(data, multiplier, params):
    test = data.copy()

    test["multiplier"] = multiplier

    test["copy"] = multiplier * test["cor_gc"]

    test.loc[~test["ideal"], "copy"] = np.nan

    segmented = segment_data(test, params)

    test["state"] = segmented["state"]

    ideal = test[test["ideal"]]

    meds = ideal.groupby("state")["copy"].apply(lambda x: {"median": np.nanmedian(x), "n": len(x)})

    meds = meds.unstack().reset_index()

    meds["fix"] = meds["state"] / meds["median"]

    # TODO: In some cases this could be NAN
    scaling_factor = multiplier * np.nanmean(meds.loc[meds["n"] > 200, "fix"])

    return scaling_factor


def segment_data(corrected_counts, params):
    hc = importr("HMMcopy")

    segmented = hc.HMMsegment(corrected_counts, params, maxiter=200)

    with (ro.default_converter + pandas2ri.converter).context():
        segmented = ro.conversion.get_conversion().rpy2py(segmented)

    segmented["state"] = segmented["state"].astype(int)

    segmented["state"] = segmented["state"] - 1

    return segmented


def get_metrics(bins, segs):
    ideal = bins[bins["ideal"]]

    state_metrics = ideal.groupby("state").apply(
        lambda x: pd.Series(
            {
                "state_mads": statsmodels.robust.mad(x["cor_gc"].dropna(), c=1),
                "state_vars": np.nanvar(x["cor_gc"]),
            }
        )
    )

    metrics = {
        "MSRSI_non_integerness": np.nanmedian(np.abs(segs["median"] - segs["state"])),
        "MBRSI_dispersion_non_integerness": np.nanmedian(np.abs(ideal["copy"] - ideal["state"])),
        "MBRSM_dispersion": np.nanmedian(np.abs(ideal["copy"] - ideal["median"])),
        "autocorrelation_hmmcopy": ideal["cor_gc"].autocorr(),
        "cv_hmmcopy": np.nanstd(ideal["cor_gc"]) / np.nanmean(ideal["cor_gc"]),
        "empty_bins_hmmcopy": np.sum(ideal["reads"] == 0),
        "mad_hmmcopy": statsmodels.robust.mad(ideal["cor_gc"].dropna(), c=1),
        "mean_hmmcopy_reads_per_bin": np.nanmean(ideal["reads"]),
        "median_hmmcopy_reads_per_bin": np.nanmedian(ideal["reads"]),
        "std_hmmcopy_reads_per_bin": np.nanstd(ideal["reads"]),
        "total_mapped_reads_hmmcopy": np.nansum(ideal["reads"]),
        "total_halfiness": np.nansum(compute_halfiness(ideal)),
        "scaled_halfiness": np.nansum(compute_scaled_halfiness(ideal)),
        "mean_state_mads": np.nanmean(state_metrics["state_mads"]),
        "mean_state_vars": np.nanmean(state_metrics["state_vars"]),
        "breakpoints": len(segs) - segs["chr"].nunique(),
        "mean_copy": np.nanmean(ideal["copy"]),
        "state_mode": ideal["state"].value_counts().index[0],
    }

    ones = ideal["state"] == 1

    # Penalize solutions which are nearly haploid
    if np.sum(ones) / len(ones) > 0.7:
        metrics["scaled_halfiness"] = np.inf

    try:
        metrics["mad_neutral_state"] = state_metrics.loc[2, "state_mads"]

    except KeyError:
        metrics["mad_neutral_state"] = np.nan

    return metrics


def compute_halfiness(bins):
    diff = np.abs(bins["median"] - bins["state"])

    diff = np.minimum(diff, 0.499)

    diff = np.abs(diff - 0.5)

    diff = diff.astype(np.float64)

    return -np.log2(diff) - 1


def compute_scaled_halfiness(bins):
    halfiness = compute_halfiness(bins)

    return halfiness / (bins["state"] + 1)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    parser.add_argument("-m", "--multipliers", default=[1, 2, 3, 4, 5, 6], nargs="+", type=int)

    parser.add_argument("-s", "--num-states", default=12, type=int)

    cli_args = parser.parse_args()

    main(cli_args)
