#!/usr/bin/env python

import pandas as pd

if __name__ == '__main__':
	ack = pd.read_csv('data/metadata.tsv', sep='\t')

	for i in sorted({x.replace('  ', ' ').replace('/', '-') for x in ack[['submitting_lab']].values.flatten()}):
	    print('- ' + i)
