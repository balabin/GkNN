#!/usr/bin/env python

summary = '''
File: Run_GkNN.py
Project: QSAR
Date: October 2016 - February 2017
Author: Ilya Balabin <ilya.balabin@gmail.com>
Purpose: Run GkNN model evaluation
Dependencies: IOCommon, Fingerprints, GkNNmodel
'''

import time, argparse

from IOCommon import *
from Fingerprints import *
from GkNNModel import *

# main parser
parser = argparse.ArgumentParser(description='Run the GkNN model')
subparsers = parser.add_subparsers()

# parser for structure-based similarities
struct_parser = subparsers.add_parser('structures')
struct_parser.add_argument('ftracts', type=str,
                           help='file with activities of training set '
                           + 'chemicals')
struct_parser.add_argument('toolkit', type=str,
                           help='chemoinformatics toolkit (\'indigo\' '
                           + 'or \'rdkit\')')
struct_parser.add_argument('fingerprint', type=str,
                           help='Chemical fingerprint: Indigo: full, sim, sub, sub-res, '
                           + 'sub-tau; RDKit (bivector): atompair, avalon, daylight, '
                           + 'maccs morgan, pharm2d, topological; RDKit: (non-bit): '
                           + 'sheridan, topotorsion; RDKit: (E-state): estate1, estate2')
struct_parser.add_argument('similarity', type=str,
                           help='Similarity metrics: Indigo: euclid-sub, tanimoto, tversky;, '
                           + 'RDKit (bivector): allbit, asymmetric, asymmetric2, '
                           + 'braunblanquet, cosine, dice,  kulczynski, mcconnaughey, '
                           + 'russel, rogotgoldberg, sokal, tanimoto, tversky; RDKit: '
                           + '(non-bit): dince, tanimoto; RDKit: (E-state): pearson (SLOW!)')
struct_parser.add_argument('ftrstr', type=str,
                           help='training set structure file')
struct_parser.add_argument('trformat', type=str,
                           help='training set structure file '
                           + 'format (smiles or sdf)')
struct_parser.add_argument('fevstr', type=str,
                           help='evaluation set structure file')
struct_parser.add_argument('evformat', type=str,
                           help='evaluation set structure file '
                           + 'format (smiles or sdf)')
struct_parser.add_argument('k', type=int,
                           help='number of nearest neighbors')
struct_parser.add_argument('x', type=float,
                           help='float parameter X')
struct_parser.add_argument('y', type=float,
                           help='float parameter Y')
struct_parser.add_argument('-q', type=str, default=None,
                           help='file with training set weights (optional)')
struct_parser.add_argument('-c', type=int, default=0,
                           help='Cross-validation: 0 (external eval), '
                           + '1 (leave-one-out), 2 or more (n-fold)')
struct_parser.add_argument('-o', type=str, default='running_sims.dat.gz',
                           help='output file (default \'running_sims.dat.gz\')')

# parser for pre-computed similarities
sims_parser = subparsers.add_parser('similarities')
sims_parser.add_argument('ftracts', type=str,
                         help='file with activities of training set '
                         + 'chemicals')
sims_parser.add_argument('fsims', type=str,
                         help='file with evaluation-to-training '
                         + 'set similarities')
sims_parser.add_argument('k', type=int,
                         help='number of nearest neighbors')
sims_parser.add_argument('x', type=float,
                         help='float parameter X')
sims_parser.add_argument('y', type=float,
                         help='float parameter Y')
sims_parser.add_argument('-q', type=str, default=None,
                         help='file with training set weights '
                         + '(optional)')
sims_parser.add_argument('-c', type=int, default=0,
                           help='Cross-validation: 0 (external eval), '
                           + '1 (leave-one-out), 2 or more (n-fold)')

if __name__ == '__main__':
    args = parser.parse_args()
    t0 = time.time()

    # obtaining molecular similarities
    if 'ftrstr' in dir(args):

        # obtaining molecular similarities
        trmols = ReadMolecules(args.ftrstr, args.toolkit, args.trformat)
        trfps = Fingerprints(trmols, args.fingerprint)
        evmols = ReadMolecules(args.fevstr, args.toolkit, args.evformat)
        evfps = Fingerprints(evmols, args.fingerprint)
        sims = Similarities(evfps, trfps, args.fingerprint, measure=args.similarity)

        # save SMILES for identifying chemicals if some fail to convert to molecule
        evsmiles = [Chem.MolToSmiles(mol) for mol in evmols]
        temp_file_smiles = 'running_smiles.dat'
        TableToFile(temp_file_smiles, evsmiles)
        print >> sys.stderr, "Saved evaluation SMILES for predicted chemicals to %s" \
              % temp_file_smiles

        # save similaities for possible future use
        temp_file_sims = 'running_sims.dat.gz'
        TableToFile(temp_file_sims, sims)
        print >> sys.stderr, "Saved similarities to %s for possible future use" \
              % temp_file_sims
        t1 = time.time()
        print >> sys.stderr, "Calculated and saved similarities in %f seconds" \
              % (t1-t0)
    elif 'fsims' in dir(args):
        print >> sys.stderr, "Using pre-computed similarities"
        sims = [map(float, row) for row in TableFromFile(args.fsims)]
        t1 = time.time()
        print >> sys.stderr, "Read similarities in %f seconds" % (t1-t0)
    else:
        print >> sys.stderr, "Unknown mode, stopping"; quit()
        
    # running GkNN model
    tracts = [map(float, row) for row in TableFromFile(args.ftracts)]
    tracts = np.array(tracts).T.tolist()[0]      ;# convert 2d table to 1d list

    # remap tracts to the range of 0-1
    tmin = np.min(tracts)
    tmax = np.max(tracts)
    tracts = (tracts-tmin)/(tmax-tmin)

    args.tracts = tracts
    args.sims = sims
    if args.q:
        q = [map(float, row) for row in TableFromFile(args.q)]
        q = np.array(q).T.tolist()[0]            ;# convert 2d table to 1d list
        args.q = q
    else:
        args.q = False

    gknn = GkNNmodel(args)
    results = gknn.GkNN()

    # remap results back
    results = np.array(results).T
    results[0] = tmin + results[0] * (tmax - tmin)
    results = results.T

    PrintTable(results)
    t2 = time.time()
    print >> sys.stderr, "Performed GkNN-based evaluation in %f seconds" % (t2-t1)
    print >> sys.stderr, "All done in %f seconds" % (t2-t0)
