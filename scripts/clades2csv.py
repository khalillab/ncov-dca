#!/usr/bin/env python


import json
import argparse


def get_options():
    description = 'From clades to a csv table'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('clades',
                        help='Clades json file')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    c = json.load(open(options.clades))
    i = 1
    print('sample,clade')
    for k in c['nodes']:
        if k.startswith('NODE_'):
            continue
        clade = c['nodes'][k]['clade_membership']
        if clade == 'unassigned':
            clade = 'unassigned_%d' % i
            i += 1
        print('%s,%s' % (k, clade))

