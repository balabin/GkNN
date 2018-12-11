#!/usr/bin/env python

summary = '''
File: Run_kNN_exp.py
Project: QSAR
Date: July 2018
Author: Ilya Balabin <balabin@gmail.com>
Purpose: Run evaluation of kNN model with exponential averaging
Dependencies: IOCommon, Fingerprints, kNN_exp
'''

import time, argparse

from IOCommon import *
from Fingerprints import *
from kNN_exp import *

# main parser
parser = argparse.ArgumentParser(description='Run the kNN_exp model')
subparsers = parser.add_subparsers()

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
                         help='exponent value')
sims_parser.add_argument('-c', type=int, default=0,
                           help='Cross-validation: 0 (external eval), '
                           + '1 (leave-one-out), 2 or more (n-fold)')

if __name__ == '__main__':
    args = parser.parse_args()
    t0 = time.time()

    # reading molecular similarities
    if 'fsims' in dir(args):
        print >> sys.stderr, "Using pre-computed similarities"
        sims = [map(float, row) for row in TableFromFile(args.fsims)]
        t1 = time.time()
        print >> sys.stderr, "Read similarities in %f seconds" % (t1-t0)
    else:
        print >> sys.stderr, "Unknown mode, stopping"; quit()
        
    # running kNN model with exponential averaging
    tracts = [map(float, row) for row in TableFromFile(args.ftracts)]
    tracts = np.array(tracts).T.tolist()[0]      ;# convert 2d table to 1d list

    # remap tracts to the range of 0-1
    tmin = np.min(tracts)
    tmax = np.max(tracts)
    tracts = (tracts-tmin)/(tmax-tmin)

    args.tracts = tracts
    args.sims = sims

    knn = kNN_exp(args)
    results = knn.kNN_exp()

    # remap results back
    results = np.array(results).T
    results[0] = tmin + results[0] * (tmax - tmin)
    results = results.T

    PrintTable(results)
    t2 = time.time()
    print >> sys.stderr, "Performed kNN_exp-based evaluation in %f seconds" % (t2-t1)
    print >> sys.stderr, "All done in %f seconds" % (t2-t0)
