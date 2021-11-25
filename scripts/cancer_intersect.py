#!/usr/bin/env python

import pandas as pd
import os, glob, os.path
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('-t', '--tmp1', help='First temp file', required=True)
parser.add_argument('-p', '--tmp2', help='Second temp file', required=True)
parser.add_argument('-o', '--output', help='Output file', required=True)

args = vars(parser.parse_args())

TF_to_P = args['tmp1']
P_to_TF = args['tmp2']
out = args['output']


intersect_1 = pd.read_csv(TF_to_P, sep="\t", header=None)
intersect_2 = pd.read_csv(P_to_TF, sep="\t", header=None)
df = intersect_2
df["4"] = intersect_1.loc[:,4]
df["5"] = intersect_1.loc[:,5]
df["6"] = intersect_1.loc[:,6]
df.to_csv(out, index=False, header=False, compression="gzip", sep ='\t')