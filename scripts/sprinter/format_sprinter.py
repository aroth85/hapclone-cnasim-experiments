import pandas as pd
import numpy as np


def main(args):

    rdr = pd.read_csv(args.in_file, delimiter="\t")

    rdr = rdr.rename(columns={"start": "START", "chrom": "CHR", "end": "END"})
    rdr["COUNT"] = rdr["Acount"] + rdr["Bcount"]

    rdr["NORM_COUNT"] = 0.1
    rdr["RAW_RDR"] = 0.1

    cells = rdr.CELL.unique()

    reads = []
    for i in range(len(cells)):
        # Get all entries with that cell
        cell = rdr[rdr["CELL"] == cells[i]]
        # Average reads across bin
        mean = cell["COUNT"].mean()
        # Create new column with reads/mean
        cell.loc[:, "NORM_COUNT"] = cell["COUNT"] / mean
        raw = (cell["NORM_COUNT"] - cell["NORM_COUNT"].mean()) + 1
        raw = np.where(raw > 1e-6, raw, 0)
        cell.loc[:, "RAW_RDR"] = raw
        reads.append(cell)

    rdr = pd.concat(reads, ignore_index=True)
    vals = rdr.START.unique()
    for i in range(len(vals)):
        norm = rdr[rdr["START"] == vals[i]]["COUNT"].sum()
        rdr.loc[rdr["START"] == vals[i], "NORM_COUNT"] = norm

    columns = ["CHR", "START", "END", "CELL", "NORM_COUNT", "COUNT", "RAW_RDR"]
    rdr = rdr.reindex(columns=columns)

    rdr["CHR"] = rdr["CHR"].astype(str)
    rdr["START"] = rdr["START"].astype(int)
    rdr["END"] = rdr["END"].astype(int)
    rdr["CELL"] = rdr["CELL"].astype(str)
    rdr["NORM_COUNT"] = rdr["NORM_COUNT"].astype(int)
    rdr["COUNT"] = rdr["COUNT"].astype(int)
    rdr["RAW_RDR"] = rdr["RAW_RDR"].astype(float)
    rdr["START"] = rdr["START"].astype(int)
    rdr["END"] = rdr["END"].astype(int)

    rdr.to_csv(args.out_file, sep="\t", index=None)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
