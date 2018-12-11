#!/usr/bin/env python

summary = '''
File: GkNN.py
Project: QSAR
Date: October 2016 - February 2017
Author: Ilya Balabin <ilya.balabin@gmail.com>
Purpose: Provide GkNN model class
Dependencies: numpy
'''

import os, sys, re
import numpy as np

# placeholder for all QSAR models
class QSARModel:
    pass

# GkNN model
class GkNNmodel(QSARModel):

    # initialize
    def __init__(self, params):
        self.tracts = np.array(params.tracts)
        self.sims = np.array(params.sims)
        self.k = params.k
        self.x = params.x
        self.y = params.y
        self.cross = params.c
        self.q = np.array(params.q)

    # run GkNN model
    def GkNN(self):

         # types of run
        if not self.cross:
            # external evaluation
            return GkNN1(self.tracts, self.sims, self.k, self.x, self.y, self.q, 0)
        elif self.cross==1:
            # leabe-one-out cross-validation
            return GkNN1(self.tracts, self.sims, self.k, self.x, self.y, self.q, 1)
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
                preds += GkNN1(tracts, sims, self.k, self.x, self.y, self.q, 0)
            return preds

# One-fold GkNN calculation
def GkNN1(tracts, sims, k, x, y, q, cross=0):

    # numbers of evaluation and training chemicals, respectively
    (m, n) = sims.shape
    
    # sanity checks
    assert len(tracts)==n                             ;# correct number of training activities
    assert k>0 and k<=n                               ;# correct number of nearest neighbors
    assert min(tracts)>=0.0                           ;# non-negative training activities
    # allow negative similarities (experimental)
##    assert min(sims.reshape(1, m*n).tolist()[0])>=0   ;# non-negative similarities
    assert x>=0.0                                      ;# positive parameter X
    assert y>=0.0                                      ;# positive parameter Y
    assert (not q) or len(q)==n                       ;# correct number of training weights
    assert (not q) or min(q)>=0.0                     ;# non-negative training weights
    if cross>=n:
        cross = 1
        print >> sys.stderr, "\n*** WARNING: Cross-validation fold out of range, " \
              + "switching to leave-one-out cross-validation ***\n"
    if cross>0:
        assert m==n

    # GkNN calculations
    preds = []
    for i in range(0, m):        
        
        # check for cross-validation
        if cross==1:
            myacts = np.delete(tracts, [i])
            mysims = np.delete(sims[i], [i])
            if q: myqs = np.delete(q, [i])
        else:
            myacts = tracts
            mysims = sims[i]
            if q: myqs = q
         
        # make estimates
        if q:            
            # weighted GkNN 
            refs = sorted(zip(myacts, mysims, myqs),
                          key=lambda t: t[1], reverse=True)[0:k]
            [myacts, mysims, myqs] = np.array(zip(*refs))
            maxsim = max(mysims)
            myacts = myacts**x
            mysims = mysims**y
            s1 = sum(myacts * mysims * myqs)
            s2 = sum(mysims * myqs)
        else:
            # flat GkNN
            refs = sorted(zip(myacts, mysims),
                          key=lambda t: t[1], reverse=True)[0:k]
            [myacts, mysims] = np.array(zip(*refs))
            maxsim = max(mysims)
            myacts = myacts**x
            mysims = mysims**y
            s1 = sum(myacts * mysims)
            s2 = sum(mysims)
        
            # store this prediction
            eps = 0.000000001
            if s2>0.0:
                # at lease one non-zero similarity
                result = [(s1/s2)**(1.0/(x+eps)), maxsim]
            else:
                # no non-zero similarities - inactive by default(?)
                result = [0.0, 0.0]
            preds.append(result)
        
    # return predictions
    return preds

if __name__ == '__main__':
    print summary

    
    
    
