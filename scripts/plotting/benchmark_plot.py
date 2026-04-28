import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from utils import plot_results

def main(args):

    in_files = [[args.hamming_file, args.hamming_outfile, 'Hamming'], 
        [args.ploidy_file, args.ploidy_outfile, 'Ploidy'],
        [args.max_file, args.max_outfile, 'Max'], 
        [args.min_file, args.min_outfile, 'Min'],  
        [args.between_file,  args.between_outfile, 'Between'], 
        [args.within_file, args.within_outfile, 'Within']]

    for files in in_files:
        data = pd.read_csv(files[0], sep = '\t', index_col='Unnamed: 0')
        xlim = np.arange(len(data) + 1) - 0.5
        ticks = np.arange(len(data))
        plot = plot_results(data, (xlim[0], xlim[-1]) , 0, files[2], ticks)
        plot.savefig(files[1])


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-p", "--ploidy-file", required=True)

    parser.add_argument("-po", "--ploidy-outfile", required=True)

    parser.add_argument("-ha", "--hamming-file", required=True)

    parser.add_argument("-hao", "--hamming-outfile", required=True)

    parser.add_argument("-ma", "--max-file", required=True)

    parser.add_argument("-mao", "--max-outfile", required=True)

    parser.add_argument("-mi", "--min-file", required=True)

    parser.add_argument("-mio", "--min-outfile", required=True)

    parser.add_argument("-b", "--between-file", required=True)

    parser.add_argument("-bo", "--between-outfile", required=True)

    parser.add_argument("-w", "--within-file", required=True)

    parser.add_argument("-wo", "--within-outfile", required=True)

    cli_args = parser.parse_args()

    main(cli_args)