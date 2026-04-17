import h5py
import numpy as np
import pandas as pd
import scipy


def main(args):
    _, _, _, ref_bins = load_data_from_file(args.in_files[0])

    in_files = dict(zip(args.cell_ids, args.in_files))

    chroms = pd.Series([x.split(":")[0] for x in ref_bins])

    regions = (chroms != chroms.shift()).cumsum().values

    valid = np.ones(len(regions)).astype(bool)

    cells, baf, rdr, reads = load_data(in_files, ref_bins)

    with h5py.File(args.out_file, "w") as fh:
        fh.create_dataset("baf", data=baf, dtype=np.int32)

        fh.create_dataset("rdr", data=rdr, dtype=np.float64)

        fh.create_dataset("reads", data=reads, dtype=np.int32)

        fh.create_dataset("bins", data=ref_bins, dtype=h5py.string_dtype(encoding="utf-8"))

        fh.create_dataset("cells", data=cells, dtype=h5py.string_dtype(encoding="utf-8"))

        fh.create_dataset("region", data=regions, dtype=np.int32)

        fh.create_dataset("valid", data=valid, dtype=np.bool_)


def load_data(in_files, ref_bins):
    cells = []

    baf = []

    rdr = []

    reads = []

    for cell_id, file_name in in_files.items():
        cells.append(cell_id)

        cell_baf, cell_rdr, cell_reads, bins = load_data_from_file(file_name)

        assert bins == ref_bins

        baf.append(cell_baf)

        rdr.append(cell_rdr)

        reads.append(cell_reads)

    baf = np.stack(baf, axis=0)

    rdr = np.stack(rdr, axis=0)

    reads = np.stack(reads, axis=0)

    return cells, baf, rdr, reads


def load_data_from_file(file_name):
    with h5py.File(file_name, "r") as fh:
        baf = fh["baf"][()]

        rdr = fh["rdr"][()]

        reads = fh["reads"][()]

        bins = [x.decode() for x in fh["bins"][()]]

    return baf, rdr, reads, bins


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-files", nargs="+", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    parser.add_argument("-c", "--cell-ids", nargs="+", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
