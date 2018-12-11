#!/usr/bin/env python

summary = '''
File: Run_GkNN.py
Project: QSAR
Date: May 2017
Author: Ilya Balabin <ilya.balabin@gmail.com>
Purpose: Calculate similarity matrix for two files
Dependencies: IOCommon, Fingerprints, GkNNModel
'''

import time, argparse

from IOCommon import *
from Fingerprints import *
from GkNNModel import *

# main parser
parser = argparse.ArgumentParser(
    description='''
    Calculate similarity matrix for two files with chemical sets
    Use Indigo fingerprints with Indigo similarity metrics
    and RDKit fingerprints with RDKit similarity metrics
    ''')

parser.add_argument('toolkit', type=str,
                    help='chemoinformatics toolkit (\'indigo\' '
                    + 'or \'rdkit\')')
parser.add_argument('fingerprint', type=str,
                    help='Chemical fingerprint: Indigo: full, sim, sub, sub-res, '
                    + 'sub-tau; RDKit (bivector): atompair, avalon, daylight, '
                    + 'maccs morgan, pharm2d, topological; RDKit: (non-bit): '
                    + 'sheridan, topotorsion; RDKit: (E-state): estate1, estate2')
parser.add_argument('similarity', type=str,
                    help='Similarity metrics: Indigo: euclid-sub, tanimoto, tversky;, '
                    + 'RDKit (bivector): allbit, asymmetric, asymmetric2, '
                    + 'braunblanquet, cosine, dice,  kulczynski, mcconnaughey, '
                    + 'russel, rogotgoldberg, sokal, tanimoto, tversky; RDKit: '
                    + '(non-bit): dince, tanimoto; RDKit: (E-state): pearson (SLOW!)')
parser.add_argument('ftrstr', type=str,
                    help='training set structure file')
parser.add_argument('trformat', type=str,
                    help='training set structure file '
                    + 'format (smiles or sdf)')
parser.add_argument('fevstr', type=str,
                    help='evaluation set structure file')
parser.add_argument('evformat', type=str,
                    help='evaluation set structure file '
                    + 'format (smiles or sdf)')
parser.add_argument('-o', type=str, default='running_sims.dat.gz',
                    help='output file (default \'running_sims.dat.gz\')')

if __name__ == '__main__':
    args = parser.parse_args()
    t0 = time.time()
    
    # obtaining molecular similarities
    trmols = ReadMolecules(args.ftrstr, args.toolkit, args.trformat)
    trfps = Fingerprints(trmols, args.fingerprint)
    evmols = ReadMolecules(args.fevstr, args.toolkit, args.evformat)
    evfps = Fingerprints(evmols, args.fingerprint)
    sims = Similarities(evfps, trfps, args.fingerprint, measure=args.similarity)
    TableToFile(args.o, sims)
    t1 = time.time()
    print >> sys.stderr, "Calculated and saved similarities in %f seconds" \
          % (t1-t0)

