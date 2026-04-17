import numpy as np
import pandas as pd


def main(args):
    in_files = get_non_failed_files(args.in_files)

    write_metrics_file(args.cell_ids, in_files, args.metrics_file)

    write_reads_file(args.cell_ids, in_files, args.reads_file, min_mappability=args.low_mappability_threshold)

    if args.paths_file is not None:
        write_paths_file(args.cell_ids, args.sample_id, in_files, args.paths_file)

    if args.segs_file is not None:
        write_segs_file(args.cell_ids, in_files, args.segs_file)


def get_non_failed_files(in_files):
    good_files = []

    for file_name in in_files:
        with pd.HDFStore(file_name, "r") as store:
            if "/failed" not in store:
                good_files.append(file_name)

    return good_files


def save_df(df, file_name):
    df.to_csv(file_name, compression="gzip", index=False)


def write_paths_file(cell_ids, sample_id, files, out_file):
    df = []

    for cell_id, file_name in zip(cell_ids, files):
        df.append((sample_id, cell_id, file_name))

    df = pd.DataFrame(df, columns=["sample_id", "cell_id", "path"])

    df.to_csv(out_file, index=False, sep="\t")


def write_metrics_file(cell_ids, files, out_file):
    df = []

    for cell_id, file_name in zip(cell_ids, files):
        multiplier = get_best_multiplier(file_name)

        with pd.HDFStore(file_name, "r") as store:
            cell_df = pd.DataFrame([store["/{}/metrics".format(multiplier)]])

            cell_df.insert(0, "cell_id", cell_id)

            df.append(cell_df)

    df = pd.concat(df)

    save_df(df, out_file)


def write_reads_file(cell_ids, files, out_file, min_mappability=0.9):
    df = []

    for cell_id, file_name in zip(cell_ids, files):
        multiplier = get_best_multiplier(file_name)

        with pd.HDFStore(file_name, "r") as store:
            cell_df = pd.merge(
                store["data"],
                store["/{}/bins".format(multiplier)][["chr", "start", "end", "median", "state"]],
                on=["chr", "start", "end"],
            )

            cell_df.insert(0, "cell_id", cell_id)

            cell_df.insert(cell_df.shape[1] - 2, "multiplier", multiplier)

            df.append(cell_df)

    df = pd.concat(df)

    df["map"] = 1  # ADDED IN TO GET CODE RUNNING
    df["is_low_mappability"] = df["map"] <= min_mappability

    save_df(df, out_file)


def write_segs_file(cell_ids, files, out_file):
    df = []

    for cell_id, file_name in zip(cell_ids, files):
        multiplier = get_best_multiplier(file_name)

        with pd.HDFStore(file_name, "r") as store:
            cell_df = store["/{}/segs".format(multiplier)]

            cell_df.insert(0, "cell_id", cell_id)

            df.append(cell_df)

    df = pd.concat(df)

    save_df(df, out_file)


def get_best_multiplier(file_name):
    with pd.HDFStore(file_name, "r") as store:
        best_multiplier = np.argmin([store["/{}/metrics".format(x)]["scaled_halfiness"] for x in store["multipliers"]])

        best_multiplier += 1

    return best_multiplier


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-a", "--sample-id", required=True)

    parser.add_argument("-c", "--cell-ids", nargs="+", required=True)

    parser.add_argument("-i", "--in-files", nargs="+", required=True)

    parser.add_argument("-m", "--metrics-file", required=True)

    parser.add_argument("-r", "--reads-file", required=True)

    parser.add_argument("-p", "--paths-file", default=None)

    parser.add_argument("-s", "--segs-file", default=None)

    parser.add_argument("-l", "--low-mappability-threshold", default=0.9, type=float)

    cli_args = parser.parse_args()

    main(cli_args)
