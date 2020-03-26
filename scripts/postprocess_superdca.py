#!/usr/bin/env python


import argparse
import pandas as pd
from Bio import SeqIO


def get_options():
    description = 'Prepare SuperDCA outputs for phylogenetic ranking'
    parser = argparse.ArgumentParser(description=description)

    #parser.add_argument('bed',
    #                    help='Alignment BED file (output from parsnp)')
    parser.add_argument('loci',
                        help='Loci file')
    parser.add_argument('couplings',
                        help='Couplings file')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    #m = pd.read_csv(options.bed, sep='\t')
    #m = m[m.columns[:2]]
    #m.columns = ['start', 'end']

    #d = {}
    #i = 1
    #for start, end in m.values:
    #    for x in range(start, end+1):
    #        d[i] = x
    #        i += 1

    d = {}
    for i, l in enumerate(open(options.loci)):
        d[int(l.rstrip())] = i+1

    for l in open(options.couplings):
        coup, start, end = l.rstrip().split()
        start = int(start)
        end = int(end)
        print('%s\t%d\t%d\t%d\t%d' % (coup,
                                      start, end,
                                      d[start], d[end]))

