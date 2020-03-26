#!/usr/bin/env python


import argparse
from Bio import SeqIO


def get_options():
    description = 'Rename raw sequences'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('input',
                        help='Input FASTA file')
    parser.add_argument('output',
                        help='Output FASTA file')

    parser.add_argument('--lsplit',
                        default='hCoV-19/',
                        help='Left split for sequence ID '
                             '(default: %(default)s)')
    parser.add_argument('--rsplit',
                        default='|',
                        help='Right split for sequence ID '
                             '(default: %(default)s)')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    seen = set()
    seqs = []
    for s in SeqIO.parse(options.input, 'fasta'):
        tmp = s.description.split(options.lsplit)[1].split(options.rsplit)[0].replace(' ', '')
        s.description = ''
        s.id = tmp
        if s.id in seen:
            continue
        seqs.append(s)
        seen.add(s.id)
    SeqIO.write(seqs, options.output, 'fasta')
