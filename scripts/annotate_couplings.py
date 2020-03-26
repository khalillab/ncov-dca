#!/usr/bin/env python


import argparse
import sys
import json
import numpy as np
import pandas as pd


def get_options():
    description = 'Annotate the ranked couplings table'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('ntaa',
                        help='nt -> aa conversion table')
    parser.add_argument('nt',
                        help='nt JSON file')
    parser.add_argument('aa',
                        help='aa JSON file')
    parser.add_argument('couplings',
                        help='Ranked couplings table')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    m = pd.read_csv(options.ntaa, sep='\t')

    nt_aa = {}
    for _, _, gene, nt_pos, aa_pos in m.values:
        nt_aa[nt_pos] = nt_aa.get(nt_pos, {})
        nt_aa[nt_pos][gene] = aa_pos

    nt = json.load(open(options.nt))

    nt_muts = {k: {int(x[1:-1]) for x in nt['nodes'][k]['muts']}
               for k in nt['nodes']
               if not k.startswith('NODE_')}

    aa = json.load(open(options.aa))

    aa_muts = {}
    for sample in aa['nodes']:
        if sample.startswith('NODE_'):
            continue
        aa_muts[sample] = {}
        for gene, muts in aa['nodes'][sample]['aa_muts'].items():
            if len(muts) > 0:
                aa_muts[sample][gene] = set()
            for mut in muts:
                aa_muts[sample][gene].add(int(mut[1:-1]))

    nt_aa_muts = {}
    for k, v in nt_muts.items():
        for m in v:
            for gene in nt_aa.get(m, ()):
                if gene not in aa_muts[k] or nt_aa[m][gene] not in aa_muts[k][gene]:
                    continue
                nt_aa_muts[m] = nt_aa_muts.get(m, set())
                nt_aa_muts[m].add((gene, nt_aa[m][gene]))

    r = pd.read_csv(options.couplings, sep='\t')

    pos = set(r[r.columns[1]]).union(r[r.columns[2]])

    aa_pos = pos.intersection(nt_aa_muts)

    r['both'] = False
    r.loc[r[(r[r.columns[1]].isin(aa_pos)) & (r[r.columns[2]].isin(aa_pos))].index, 'both'] = True

    genes = set()
    for gene in sorted({gene for x in aa_pos
                        for gene, pos in nt_aa_muts[x]}):
        r[gene + '_1'] = np.nan
        r[gene + '_2'] = np.nan
        genes.add(gene + '_1')
        genes.add(gene + '_2')

    for i, v in r.iterrows():
        for gene, pos in nt_aa_muts.get(int(v[1]), ()):
            r.loc[i, gene + '_1'] = pos
        for gene, pos in nt_aa_muts.get(int(v[2]), ()):
            r.loc[i, gene + '_2'] = pos

    r.to_csv(sys.stdout, sep='\t', index=False)
