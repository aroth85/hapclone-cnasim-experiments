import h5py
import numpy as np
import pandas as pd


def main(args):
    baf, bins = get_baf(args.bin_size, args.block_size, args.snps_file)

    rdr, reads = get_rdr(bins, args.cell_id, args.reads_file)

    bin_ids = bins.apply(lambda row: "{chrom}:{start}:{end}".format(**row.to_dict()), axis=1)

    with h5py.File(args.out_file, "w") as f:
        f.create_dataset("baf", data=baf, dtype=np.int32)

        f.create_dataset("rdr", data=rdr, dtype=np.float64)

        f.create_dataset("reads", data=reads, dtype=np.int32)

        f.create_dataset("bins", data=bin_ids, dtype=h5py.string_dtype(encoding="utf-8"))

        f.create_dataset("cells", data=[args.cell_id], dtype=h5py.string_dtype(encoding="utf-8"))


def get_baf(bin_size, block_size, snps_file):
    assert bin_size % block_size == 0

    num_blocks = bin_size // block_size

    df = pd.read_csv(snps_file, sep="\t")

    out_df = []

    bin_idx = 0

    bins = []

    for (chrom, start, end), bin_df in df.groupby(["chrom", "start", "end"]):
        bins.append({"chrom": chrom, "start": start, "end": end, "bin_idx": bin_idx})

        for i in range(num_blocks):
            block_start = start + i * block_size

            block_end = block_start + block_size

            block_df = bin_df[bin_df["pos"].between(block_start, block_end)]

            row = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "bin_idx": bin_idx,
                "block_idx": i,
                "a": block_df["a"].sum(),
                "b": block_df["b"].sum(),
            }

            out_df.append(row)

        bin_idx += 1

    bins = pd.DataFrame(bins)

    out_df = pd.DataFrame(out_df)

    baf = np.stack(
        [
            out_df.pivot(index="bin_idx", columns="block_idx", values="a").values,
            out_df.pivot(index="bin_idx", columns="block_idx", values="b").values,
        ],
        axis=-1,
    )

    return baf, bins


def get_rdr(bins, cell_id, reads_file):
    df = pd.read_csv(reads_file, sep="\t")

    df = df.rename(columns={"CELL": "cell_id"})

    df = df[df["cell_id"] == cell_id]

    df = pd.merge(df, bins, on=["chrom", "start", "end"])

    df = df.sort_values(by="bin_idx")

    df["reads"] = df["Acount"] + df["Bcount"]

    reads = df.pivot(index="cell_id", columns="bin_idx", values="reads").values

    reads = np.squeeze(reads)

    reads = np.expand_dims(reads, axis=1)

    df["rdr"] = df["reads"] / df["reads"].mean()

    rdr = df.pivot(index="cell_id", columns="bin_idx", values="rdr").values

    rdr = np.squeeze(rdr)

    rdr = np.expand_dims(rdr, axis=1)

    return rdr, reads


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--cell-id", required=True)

    parser.add_argument("-r", "--reads-file", required=True)

    parser.add_argument("-s", "--snps-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    parser.add_argument("-b", "--bin-size", default=int(5e5), type=int)

    parser.add_argument("-k", "--block-size", default=int(1e4), type=int)

    cli_args = parser.parse_args()

    main(cli_args)
