#!/usr/bin/env python

summary = '''
File: IOCommon.py
Project: QSAR
Date: October 2016
Author: Ilya Balabin <ilya.balabin@gmail.com>
Purpose: Routines for reading and writing 2d lists
Usage: Call functions from Python
Dependencies: indigo, rdkit, numpy
'''

import sys, os, re, gzip, fileinput
import numpy as np

# read a table from a (gzipped) file
def TableFromFile(infile):
    try:
        if re.match('^.+\.gz$', infile):
            fin = gzip.open(infile, 'r')
        else:
            fin = open(infile, 'r')
    except:
        print "Can not open %s for reading, stopping" % infile
        return None
    table = []
    for line in fin.readlines():
        if re.match('^\s?#.*$', line) or re.match('^\s?$', line): continue
        if re.match('^\s?$', line) or re.match('^\s?$', line): continue
        table.append(line.strip().split())
    fin.close()
    return table

# write table to file
def TableToFile(outfile, table, overwrite=False):
    if overwrite and os.path.exists(outfile):
        print "Can not overwrite %s, stopping" % outfile
        return False
    try:
        if re.match('^.+\.gz$', outfile):
            fout = gzip.open(outfile, 'w')
        else:
            fout = open(outfile, 'w')
    except:
        print "Can not open %s for writing, stopping" % outfile
        return False
    for line in table:
        for field in line:
            fout.write(str(field) + " \t")
        fout.write("\n")
    fout.close()
    return True


# print out a table
def PrintTable(table):
    for line in table:
        for field in line:
            print "%s \t" % field,
        print ""
    return True

if __name__ == '__main__':
    print summary

