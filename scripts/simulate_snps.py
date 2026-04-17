import numpy as np
import pandas as pd


def main(args):
    reads = pd.read_csv(args.in_file, delimiter="\t")

    chrom_df = load_chrom_df(args.num_snps, reads)

    out_df = []

    phase = 0

    phase_set = 0

    for c, row in chrom_df.iterrows():
        bin_df = reads.loc[reads["chrom"] == c, ["start", "end"]].drop_duplicates()

        bin_df = bin_df.sort_values(by=["start", "end"])

        snp_pos = np.random.choice(row["length"], row["num_snps"], replace=False)

        snp_pos = pd.Series(sorted(snp_pos))

        for b, e in zip(bin_df["start"], bin_df["end"]):
            for p in snp_pos[snp_pos.between(b, e)]:
                if np.random.uniform() < args.phase_switch_prob:
                    phase = 1 - phase

                    phase_set += 1

                out_df.append(
                    {
                        "chrom": c,
                        "start": b,
                        "end": e,
                        "pos": p,
                        "phase": phase,
                        "phase_set": phase_set,
                    }
                )

    out_df = pd.DataFrame(out_df)

    out_df.to_csv(args.out_file, index=False, sep="\t")


def load_chrom_df(num_snps, reads):
    df = []

    for c in reads.chrom.unique():
        df.append(
            {
                "chrom": c,
                "length": reads.loc[reads["chrom"] == c, "end"].max(),
            }
        )

    df = pd.DataFrame(df).set_index("chrom")

    p = df["length"] / df["length"].sum()

    df["num_snps"] = np.random.multinomial(num_snps, p)

    return df


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    parser.add_argument("-n", "--num_snps", default=int(2e6), type=int)

    parser.add_argument("-p", "--phase-switch-prob", default=0.01, type=float)

    cli_args = parser.parse_args()

    main(cli_args)
