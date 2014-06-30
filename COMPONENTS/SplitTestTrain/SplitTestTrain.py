#!/usr/bin/env python
import component_skeleton.main
import os, re


def splitFile(infile, odd_file, even_file):
    """
    This function splits alignments (from AlignPeaks component) or Sequences (from GetSequences or AlignPeaks) into two halfs (training and test set).
    Both cotain equally high Z scoring (or peak heights) regions. Heights are stored in the sequence identifier line: >>hg19_chr3_20227683_20227835_reg1000069.p2_88.427_+
    """

    infile_dict = {}

    #read input file into dictionary
    for line in open(infile):
        if line.startswith('>>'):
            curr = line
            infile_dict[line] = []
        else:
            infile_dict[curr].append(line)



    odd = open(odd_file,'w')
    even = open(even_file,'w')

    i = 0
    for headline in sorted(infile_dict.keys(), key = lambda k: float(k.split('_')[-2]), reverse=True):
        if i%2 == 1:
            odd.write(headline)
            for line in infile_dict[headline]:
                odd.write(line)
        if i%2 == 0:
            even.write(headline)
            for line in infile_dict[headline]:
                even.write(line)

        i += 1


    odd.close()
    even.close()


def execute(cf):
    infile = cf.get_input("infile") #result file from AlignPeaks
    odd_file = cf.get_output("odd_file")
    even_file = cf.get_output("even_file")  

    splitFile(infile, odd_file, even_file)
    
    return 0

component_skeleton.main.main(execute)
                                                                 
