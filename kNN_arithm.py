#!/usr/bin/env python

summary = '''
File: kNN_arithm.py
Project: QSAR
Date: Juky 2017
Author: Ilya Balabin <ilya.balabin@gmail.com>
Purpose: Provide kNN_arithm model class
Dependencies: numpy
'''

import os, sys, re
import numpy as np

# placeholder for all QSAR models
class QSARModel:
    pass

# kNN_arithm model
class kNN_arithm(QSARModel):

    # initialize
    def __init__(self, params):
        self.tracts = np.array(params.tracts)
        self.sims = np.array(params.sims)
        self.k = params.k
        self.cross = params.c

    # run kNN_arithm model
    def kNN_arithm(self):

         # types of run
        if not self.cross:
            # external evaluation
            return kNN1_arithm(self.tracts, self.sims, self.k, 0)
        elif self.cross==1:
            # leabe-one-out cross-validation
            return kNN1_arithm(self.tracts, self.sims, self.k, 1)
        else:
            # N-fold cross-validation
            n = len(self.tracts)
            nfolds = self.cross
            dn = int(n/nfolds)
            preds = []
            for j in range(0, nfolds):
                n1 = j * dn
                n2 = n1 + dn
                if n2+dn>n: n2 = n
                tracts = np.delete(self.tracts, range(n1, n2), axis=0)
                sims = np.delete(np.delete(self.sims, range(n1, n2), axis=1),
                                 range(0, n1)+range(n2, n), axis=0)
                preds += kNN1_arithm(tracts, sims, self.k, 0)
            return preds

# One-fold kNN_arithm calculation
def kNN1_arithm(tracts, sims, k, cross=0):

    # numbers of evaluation and training chemicals, respectively
    (m, n) = sims.shape
    
    # sanity checks
    assert len(tracts)==n                             ;# correct number of training activities
    assert k>0 and k<=n                               ;# correct number of nearest neighbors
    assert min(tracts)>=0.0                           ;# non-negative training activities
    # allow negative similarities (experimental)
##    assert min(sims.reshape(1, m*n).tolist()[0])>=0   ;# non-negative similarities
    if cross>=n:
        cross = 1
        print >> sys.stderr, "\n*** WARNING: Cross-validation fold out of range, " \
              + "switching to leave-one-out cross-validation ***\n"
    if cross>0:
        assert m==n

    # kNN_arithm calculations
    preds = []
    for i in range(0, m):        
        
        # check for cross-validation
        if cross==1:
            myacts = np.delete(tracts, [i])
            mysims = np.delete(sims[i], [i])
        else:
            myacts = tracts
            mysims = sims[i]
         
        # make estimates
        refs = sorted(zip(myacts, mysims),
                      key=lambda t: t[1], reverse=True)[0:k]
        [myacts, mysims] = np.array(zip(*refs))
        maxsim = max(mysims)
        result = [sum(myacts)/float(k), maxsim]
        preds.append(result)
        
    # return predictions
    return preds

if __name__ == '__main__':
    print summary

    
    
    
