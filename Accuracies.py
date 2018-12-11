#!/usr/bin/env python

summary = '''
File: Accuracies.py
Project: QSAR
Date: October 2016; modified as requested by the reviewers July 2018
Author: Ilya Balabin <balabin.i@epa.gov>
Purpose: Calculate ROC curve and the following metrics as functions of cutoff:
  TP             exp>thr & pred>thr      |---
  FP             exp<=thr & pred>thr     |
  TN             exp<=thr & pred<=thr    |
  FN             exp>thr & pred<=thr     |---calculated using sklearn
  sensitivity    TPR = TP/(TP+FN)
  specificity    TNR = TN/(TN+FP)
  bal. accuracy  0.5*(TPR+TNR)
  accuracy       (TP+TN)/(TP+FP+TN+FN)
  precision      PPV = TP/(TP+FP)
  neg.pred.val.  NPV = TN/(TN+FN)
  (recall == sensitivity)
Usage: Call functions from Python
Dependencies: sklearn, IOCommon
'''

from sklearn.metrics import roc_curve, auc, confusion_matrix
from sklearn.preprocessing import binarize
from IOCommon import *

# calculate ROC curve (actual - list of ints, preds - list of floats)
def ROCCurve(actual, preds):

    assert len(actual)==len(preds)
    fpr, tpr, thr = roc_curve(actual, preds)
    #roc_auc = auc(fpr, tpr)
    return [tpr, fpr, thr]

# calculate balanced accuracy (actual = classes; must have both 0 and 1)
def Accuracies(actual, preds, DEBUG=False):

    assert actual.ndim==1   ;# must be 1d NumPy array
    assert preds.ndim==1    ;# must be 1d NumPy array
    assert set(actual)==set([0, 1])   ;# includes both active and inactive
    assert len(preds)==len(actual)
    counter = 0
    accuracies = []

    # infer thresholds from predictions
    #thr = sorted(set(preds.tolist() + [0.0] + [1.0]))
    thr = sorted(set(preds.tolist()))
    preds = preds.reshape(1, -1)   ;# need 2d arrays here
    for th in thr:
        binpreds = binarize(preds, th)[0]
        [[tn, fp], [fn, tp]] = \
              [map(float, row) for row in confusion_matrix(actual, binpreds)]
        if not (tn+fp)*(tp+fn): continue   ;# avoid division by zero

        # sensiticity a.k.a. recall
        sens = float(tp) / (tp + fn)
        
        # specificity
        spec = float(tn) / (tn + fp)
        
        # balanced accuracy
        bacc = 0.5 * (sens + spec)
        
        # accuracy
        acc = float(tp + tn) / (tp + fp + tn + fn)

        # precision a.k.a. positive prediction value
        if tp+fp>0:
            pres = float(tp) / (tp + fp)
        else:
            pres = 0.0
        
        # negative predicted value
        if tn+fn>0:
            npv = float(tn) / (tn + fn)
        else:
            npv = 0.0

        accuracies.append([th, sens, spec, bacc, acc, pres, npv])
        if DEBUG and not counter%100: print counter
        counter +=1
    return accuracies

if __name__ == '__main__':
    print summary


