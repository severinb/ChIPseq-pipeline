#!/usr/bin/env python

import component_skeleton.main
import sys, os
from string import *

def execute(cf):

    infile = cf.get_input("in_file") #format: chr start end strand (refined_peaks)
    genome = cf.get_parameter("genome", "string")
    genome_path = cf.get_parameter("genome_path", "string")
    outfile = cf.get_output("Sequences") #format genome_chr_start_end_strand \n sequence

    if genome == "hg18":
        genomeDB = "/import/bc2/data/databases/UCSC/hg18/"
    if genome == "hg19":
        genomeDB = "/import/bc2/data/databases/UCSC/hg19/"
    if genome == "dm3":
        genomeDB = "/import/bc2/data/databases/UCSC/dm3/"
    if genome == "mm9":
        genomeDB = "/import/bc2/data/databases/UCSC/mm9/"

    if genome_path:
        genomeDB = genome_path

    SequenceDict = {}  #chr: [chr, start,end,strand,sequence]

    #parse sequences
    f = open(infile)
    for line in f:
        tmp = line.split()
        try:
            SequenceDict[tmp[0]].append(tmp)
        except KeyError:
            SequenceDict[tmp[0]] = []
            SequenceDict[tmp[0]].append(tmp)

    f.close()

    o = open(outfile, 'w')

    #get sequences by loading one chromosome at a time into memory
    for chrom in SequenceDict:
        print chrom
        c = open(os.path.join(genomeDB, chrom + '.fa'))
        c.readline() #skip >chr on first line
        bases = []
        for line in c:
            bases += list(line.strip())
        
        for peak in SequenceDict[chrom]:
            ID = '_'.join([genome] + peak)
            start = int(peak[1])
            end = int(peak[2])
            o.write('>>' + ID + '\n')
            o.write(''.join(bases[start-1:end]).upper() + '\n')
            
        c.close()

    o.close()

    return 0

component_skeleton.main.main(execute)
