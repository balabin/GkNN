# GkNN
Generalized kNN model (https://doi.org/10.1186/s13321-018-0300-0)

These scripts implement the generalized k-nearest neighbor model described in the above article as well as three other kNN-based models published earlier. The functions of these scripts are briefly described below.

Accuracies.py: routines for calculating numbers of TP, FP, FN, and TN as well as sensitivity, specificity, balanced accuracy, accuracy, precision, negative predicted value, and ROC AUC and generating the ROC curve for actual and predicted labels..
Fingerprints.py: routines for reading lists of chemicals in the SDF or SMILES formats, calculating all fingerprint types supported by RDKit and Indigo, and calculating pairwise similarities between two sets of fingerprints.
Get_Similarities.py: a wrapper for calculating a similarity matrix between two sets of chemicals.
GkNNModel.py: the generalized k-nearest neighbor model class
IOCommon.py: a wrapper around I/O routines
kNN_arithm.py: the arithmetic kNN model class
kNN_exp.py: the exponential kNN model class
kNN_geom.py: the geometric kNN model class
Run_GkNN.py: a wrapper for running the generalized kNN model
Run_kNN_arithm.py: a wrapper for running the arithmetic kNN model 
Run_kNN_exp.py: a wrapper for running the exponential kNN model
Run_kNN_geom.py: a wrapper for running the geometric kNN model

The estrogen receptor data are available from ftp://newftp.epa.gov/COMPTOX/Sustainable_Chemistry_Data/CERAPP_QSAR_Models/ .
