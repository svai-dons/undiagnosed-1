#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import os
from tqdm import tqdm
# !pip install tqdm

# grep -v '##' results/SQ9887_L00.combined.genotyped-snpeff-grch37.75.vcf | head -n1

variants = pd.read_csv('results/SQ9887_L00.combined.genotyped-snpeff-grch37.75.vcf', sep='\t',
                       index_col=None,
                       header=None,
                      comment='#')

variants.columns = ['#CHROM',
 'POS',
 'ID',
 'REF',
 'ALT',
 'QUAL',
 'FILTER',
 'INFO',
 'FORMAT',
 'RESULT']


ANNKEYS=['Allele',
 'Annotation',
 'Annotation_Impact',
 'Gene_Name',
 'Gene_ID',
 'Feature_Type',
 'Feature_ID',
 'Transcript_BioType',
 'Rank',
 'HGVS.c',
 'HGVS.p',
 'cDNA.pos.length',
 'CDS.pos.length',
 'AA.pos.length',
 'Distance',
 'ERRORS']

def parse_info(x):
    outdict={}
    for chunk in x.split(';'):
        kk,vv = chunk.split('=')
        outdict[kk] = vv
    if 'ANN' in outdict:
        annotations = []
        for ann in outdict['ANN'].split(','):
            annotations.append(dict(zip(ANNKEYS,ann.split('|'))))
        outdict['ANN'] = annotations
    return outdict


outfile_ann = 'results/SQ9887_L00.combined.genotyped-snpeff-grch37.75-annotations.vcf'
with open(outfile_ann, 'w+') as fh:
    first = True
    for kk,vv in tqdm(variants.iterrows()):

        tmp_bed = vv.loc['#CHROM':'RESULT']
        tmp_bed = tmp_bed.drop('INFO')
        info = parse_info(vv['INFO'])
        if 'ANN' not in info:
            continue
        tmp_ann = pd.DataFrame(info['ANN'])

        tmp_ann = pd.concat([pd.concat([tmp_bed.to_frame().T]*len(tmp_ann),).reset_index(drop=True),
                              tmp_ann, 
                              ], 1)
        
        tmp_ann.to_csv(fh,
            header=first, index=False,
            mode='a', sep='\t')
        
        if first:
            first=None
