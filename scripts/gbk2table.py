#!/usr/bin/env python


import argparse
from Bio import SeqIO


def get_options():
    description = 'Genbank to NT -> AA table'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('genbank',
                        help='Input GenBank file')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    print('molecule\tstrand\tname\tnt\taa')
    for s in SeqIO.parse(options.genbank, 'genbank'):
        for f in s.features:
            if f.type != 'CDS':
                continue
            if 'locus_tag' in f.qualifiers:
                gene = f.qualifiers['locus_tag'][0]
            else:
                gene = f.qualifiers['gene'][0]
            if f.strand > 0:
                pos = range(f.location.start+1, f.location.end+1)
                strand = '+'
            else:
                pos = range(f.location.end+1, f.location.start+1, -1)
                strand = '-'
            aa = 1
            for i, nt in enumerate(pos):
                print('\t'.join([s.id, strand, gene, str(nt), str(aa)]))
                if not (i+1) % 3:
                    aa += 1
