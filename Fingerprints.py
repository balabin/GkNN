#!/usr/bin/env python

summary = '''
File: Fingerprints.py
Project: QSAR
Date: May 2017
Author: Ilya Balabin <ilya.balabin@gmail.com>
Purpose: Read chemicals, calculate fingerprints and similarities
Usage: Call functions from Python
Dependencies: indigo, rdkit, numpy
'''

import re, sys, os, time
import numpy as np

from rdkit import Avalon
from rdkit import Chem
from rdkit import DataStructs
from rdkit import Geometry
from rdkit import RDConfig
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import Descriptors
from rdkit.Chem import Fingerprints
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Sheridan
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem.EState import Fingerprinter
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D
from rdkit.Chem.Pharm2D import Generate
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint

from indigo import *
indigo = Indigo()

# acceptable fingerprints
indigofps = ['full', 'sim', 'sub', 'sub-res', 'sub-tau']
rdkitfps = ['atompair', 'avalon', 'daylight', 'maccs', 'morgan', 
            'pharm2d', 'topological']
rdkitnonbitfps = ['sheridan', 'topotorsion']
rdkitestatefps = ['estate1', 'estate2']

# acceptable similarity metrics
indigosims = ['euclid-sub', 'tanimoto', 'tversky']
rdkitsims = ['allbit', 'asymmetric', 'asymmetric2', 'braunblanquet',
             'cosine', 'dice',  'kulczynski', 'mcconnaughey', 'russel',
             'rogotgoldberg', 'sokal', 'tanimoto', 'tversky']
rdkitnonbitsims = ['dice', 'tanimoto']
rdkitestatesims = ['pearson']
    
# read molecules
def ReadMolecules(infile, toolkit, format):

    # sanity checks
    assert toolkit=='indigo' \
           or toolkit=='rdkit'
    assert format=='smiles' \
           or format=='sdf'
    
    # reading SMILES file
    if format=='smiles':
        if re.match('^.+\.gz$', infile):
            fin = gzip.open(infile, 'r')
        else:
            fin = open(infile, 'r')
        smiles = [line.strip() for line in fin.readlines()]
        fin.close()
        if toolkit=='indigo':
            results = [indigo.loadMolecule(smile) for smile in smiles]
        elif toolkit=='rdkit':
            results = [Chem.MolFromSmiles(smile) for  smile in smiles]
    
    # reading SDF file
    elif format=="sdf":
        if toolkit=='indigo':
            results = [mol for mol in indigo.iterateSDFile(infile)]
        elif toolkit=='rdkit':
            results = [mol for mol in Chem.SDMolSupplier(infile)]
    
    # return results
    return [mol for mol in results if mol is not None]


# calculate fingerprints ('indigo_*' for Indigo mols; others for RDKit mols)
def Fingerprints(mols, fingerprint):

    # Indigo fingerprints
    if fingerprint in indigofps:
        return [mol.fingerprint(fingerprint) for mol in mols]

    # RDKit fingerprints
    if fingerprint in rdkitfps:
        if fingerprint=="atompair":
            return [Pairs.GetAtomPairFingerprintAsBitVect(mol) for mol in mols]
        elif fingerprint=="avalon":
            return [pyAvalonTools.GetAvalonFP(mol) for mol in mols]
        elif fingerprint=="daylight":
            return [Chem.RDKFingerprint(mol, fpSize=2048) for mol in mols]
        elif fingerprint=="maccs":
            return [MACCSkeys.GenMACCSKeys(mol) for mol in mols]
        elif fingerprint=="morgan":
            return [(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024))
                    for mol in mols]
        elif fingerprint=="pharm2d":
            return [Generate.Gen2DFingerprint(mol, Gobbi_Pharm2D.factory)
                    for mol in mols]
        elif fingerprint=="topological":
            return [FingerprintMols.FingerprintMol(mol) for mol in mols]

    # RDKit non-bit (integer or float) fingerprints
    if fingerprint in rdkitnonbitfps:
        if fingerprint=="sheridan":
            return [Sheridan.GetBPFingerprint(mol) for mol in mols]
        elif fingerprint=="topotorsion":
            return [Torsions.GetTopologicalTorsionFingerprint(mol)
                    for mol in mols]        

    # E-state fingerprints 
    if fingerprint in rdkitestatefps:
        if fingerprint=="estate1":
            return [Fingerprinter.FingerprintMol(mol)[0] for mol in mols]
        elif fingerprint=="estate2":
            return [Fingerprinter.FingerprintMol(mol)[1] for mol in mols]
        
    # unknown fingerprint
    return None

# calculate similarities
def Similarities(fps1, fps2, fingerprint, measure='tanimoto'):

    # Indigo fingerptins
    if fingerprint in indigofps and measure in indigosims:
        return np.array([[indigo.similarity(fp1, fp2, measure)
                          for fp2 in fps2] for fp1 in fps1])

    # RDKit fingerprints
    if fingerprint in rdkitfps and measure in rdkitsims:

        if measure=='allbit':
            return np.array([[DataStructs.FingerprintSimilarity(
                            fp1, fp2, metric=DataStructs.AllBitSimilarity)
                              for fp2 in fps2] for fp1 in fps1])
        elif measure=='asymmetric':
            return np.array([[DataStructs.FingerprintSimilarity(
                            fp1, fp2, metric=DataStructs.AsymmetricSimilarity)
                              for fp2 in fps2] for fp1 in fps1])
        elif measure=='asymmetric2':
            return np.array([[DataStructs.FingerprintSimilarity(
                            fp2, fp1, metric=DataStructs.AsymmetricSimilarity)
                              for fp2 in fps2] for fp1 in fps1])
        elif measure=='braunblanquet':
            return np.array([[DataStructs.FingerprintSimilarity(
                            fp1, fp2, metric=DataStructs.BraunBlanquetSimilarity)
                              for fp2 in fps2] for fp1 in fps1])
        elif measure=='cosine':
            return np.array([[DataStructs.FingerprintSimilarity(
                            fp1, fp2, metric=DataStructs.CosineSimilarity)
                              for fp2 in fps2] for fp1 in fps1])
        elif measure=='dice':
            return np.array([[DataStructs.FingerprintSimilarity(
                            fp1, fp2, metric=DataStructs.DiceSimilarity)
                              for fp2 in fps2] for fp1 in fps1])
        elif measure=='kulczynski':
            return np.array([[DataStructs.FingerprintSimilarity(
                            fp1, fp2, metric=DataStructs.KulczynskiSimilarity)
                              for fp2 in fps2] for fp1 in fps1])
        elif measure=='mcconnaughey':
            return np.array([[DataStructs.FingerprintSimilarity(
                            fp1, fp2, metric=DataStructs.McConnaugheySimilarity)
                              for fp2 in fps2] for fp1 in fps1])
        elif measure=='rogotgoldberg':
            return np.array([[DataStructs.FingerprintSimilarity(
                            fp1, fp2, metric=DataStructs.RogotGoldbergSimilarity)
                              for fp2 in fps2] for fp1 in fps1])
        elif measure=='russel':
            return np.array([[DataStructs.FingerprintSimilarity(
                            fp1, fp2, metric=DataStructs.RusselSimilarity)
                              for fp2 in fps2] for fp1 in fps1])
        elif measure=='sokal':
            return np.array([[DataStructs.FingerprintSimilarity(
                            fp1, fp2, metric=DataStructs.SokalSimilarity)
                              for fp2 in fps2] for fp1 in fps1])
        elif measure=='tanimoto':
            return np.array([[DataStructs.FingerprintSimilarity(
                    fp1, fp2, metric=DataStructs.TanimotoSimilarity)
                       for fp2 in fps2] for fp1 in fps1])
        elif measure=='tversky':
            return np.array([[DataStructs.TverskySimilarity(
                            fp1, fp2, a=0.5, b=0.5)
                              for fp2 in fps2] for fp1 in fps1])

    # RDKit non-bit (integer or float) fingerprints
    if fingerprint in rdkitnonbitfps and measure in rdkitnonbitsims:

        if measure=='dice':
            return np.array([[DataStructs.DiceSimilarity(fp1, fp2)
                              for fp2 in fps2] for fp1 in fps1])
        elif measure=='tanimoto':
            return np.array([[DataStructs.TanimotoSimilarity(fp1, fp2)
                              for fp2 in fps2] for fp1 in fps1])
                            
    # E-state fingerprints
    if fingerprint in rdkitestatefps and measure in rdkitestatesims:

        if measure=='pearson':
            sims = np.array([[np.corrcoef(fp1, fp2)[0][1] for fp2 in fps2]
                    for fp1 in fps1])
            # check for nan's
            if sum([sum(np.isnan(sim)) for sim in sims])==0:
                return sims
            else:
                return None

    # unknown fingerprint and/or similarity measure
    return None

if __name__ == '__main__':
    print summary

