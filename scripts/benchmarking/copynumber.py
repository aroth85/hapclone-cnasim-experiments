import pandas as pd
import numpy as np
from pathlib import Path
from utils import *


def main(args):

    ploidy = []
    hamming = []

    methods, hapclone_files = arrange_hapclone(args.hapclone_files, args.num_replicates)
    methods.extend(["Chisel", "Signals", "HMMCopy"])

    for f in range(args.num_replicates):

        hapclone_file = hapclone_files[f]
        chisel_file = args.chisel_files[f]
        signals_file = args.signals_files[f] 
        cnasim_file = args.cnasim_files[f] 
        hmmcopy_file = args.hmmcopy_files[f] 

        cnasim = load_cnasim_profile(cnasim_file)

        cells = cnasim.CELL.unique()
        ploidy_results = []
        hamming_results = []

        # HapClone
        for file in hapclone_file:
            hapclone = load_hapclone_results(file)
            hapclone_results = benchmark_copynumber(
                hapclone, cnasim, "cell_id", "total", cells
            )
            ploidy_results.append(hapclone_results[0])
            hamming_results.append(hapclone_results[1])

        # Chisel
        if Path(chisel_file).exists():
            chisel = load_chisel_results(chisel_file)
            chisel_results = benchmark_copynumber(
                chisel, cnasim, "CELL", "total", cells
            )
            ploidy_results.append(chisel_results[0])
            hamming_results.append(chisel_results[1])
        else:
            ploidy_results.append( np.nan)
            hamming_results.append( np.nan)

        # Signals
        if Path(signals_file).exists():
            signals = load_signals_results(signals_file)
            signals_results = benchmark_copynumber(
                signals, cnasim, "cell_id", "state", cells
            )
            ploidy_results.append(signals_results[0])
            hamming_results.append(signals_results[1])
        else:
            ploidy_results.append( np.nan)
            hamming_results.append( np.nan)

        # HMMCopy
        if Path(hmmcopy_file).exists():
            hmmcopy = load_hmmcopy_results(hmmcopy_file)
            hmmcopy_results = benchmark_copynumber(
                hmmcopy, cnasim, "cell_id", "state", cells
            )
            ploidy_results.append(hmmcopy_results[0])
            hamming_results.append(hmmcopy_results[1])
        else:
            ploidy_results.append( np.nan)
            hamming_results.append( np.nan)

        hamming.append(hamming_results)
        ploidy.append(ploidy_results)

    save_results(hamming, methods, args.hamming_results)
    save_results(ploidy, methods, args.ploidy_results)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-ha", "--hapclone-files", nargs="+", required=True)

    parser.add_argument("-c", "--chisel-files", nargs="+", required=True)

    parser.add_argument("-s", "--signals-files", nargs="+", required=True)

    parser.add_argument("-hm", "--hmmcopy-files", nargs="+", required=True)

    parser.add_argument("-p", "--cnasim-files", nargs="+", required=True)

    parser.add_argument("-po", "--ploidy-results", required=True)

    parser.add_argument("-ho", "--hamming-results", required=True)

    parser.add_argument("-n", "--num-replicates", required=True, type=int)

    cli_args = parser.parse_args()

    main(cli_args)
