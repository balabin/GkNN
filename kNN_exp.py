#!/usr/bin/env python

summary = '''
File: kNN_exp.py
Project: QSAR
Date: Juky 2017
Author: Ilya Balabin <ilya.balabin@gmail.com>
Purpose: Provide kNN_exp model class
Dependencies: numpy
'''

import os, sys, re
import numpy as np

# placeholder for all QSAR models
class QSARModel:
    pass

# kNN_exp model
class kNN_exp(QSARModel):

    # initialize
    def __init__(self, params):
        self.tracts = np.array(params.tracts)
        self.sims = np.array(params.sims)
        self.k = params.k
        self.x = params.x
        self.cross = params.c

    # run kNN_exp model
    def kNN_exp(self):

         # types of run
        if not self.cross:
            # external evaluation
            return kNN1_exp(self.tracts, self.sims, self.k, self.x, 0)
        elif self.cross==1:
            # leabe-one-out cross-validation
            return kNN1_exp(self.tracts, self.sims, self.k, self.x, 1)
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
                preds += kNN1_exp(tracts, sims, self.k, self.x, 0)
            return preds

# One-fold kNN_exp calculation
def kNN1_exp(tracts, sims, k, x, cross=0):

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

    # kNN_exp calculations
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
        eps = 0.000000000000001     ;# avoiding division by zero
        mydists = np.array([max(0.0, 1.0/(sim+eps)-1.0) for sim in mysims])
        s1 = sum(myacts * np.exp(-x*mydists))
        s2 = sum(np.exp(-x*mydists))
        result = [s1/s2, maxsim]
        preds.append(result)
        
    # return predictions
    return preds

if __name__ == '__main__':
    print summary

    
    
    
