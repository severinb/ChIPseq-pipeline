#!/usr/bin/env python

import os
from string import *
import component_skeleton.main


def execute(cf):
    infile = cf.get_input("infile")
    outfile = cf.get_output("outfile")
    length = cf.get_parameter("length", "int")

    # chrX    139917566       139917652       reg1000025.p1   3000.681        +
    # chrM    8666    8753    reg1000145.p1   1917.321        +
    # chr1    24620664        24620750        reg1000060.p2   1889.144        +
    # chr9    24346496        24346587        reg1000049.p2   1679.997        +
    # chrM    9026    9113    reg1000145.p2   1102.568        +
    # chr1    24620300        24620387        reg1000060.p1   1034.005        +

    o = open(outfile, 'w')

    for line in open(infile):
        t = line.strip().split()
        m = int((float(t[1]) + float(t[2]))/2.0)
        ml = str(max(0, m - int(length/2.0))) #to not get negative coordinates
        mr = str(m + int(length/2.0))

        o.write('\t'.join([t[0], ml, mr, t[3], t[4], t[5]]) + '\n')


    o.close()
 
    return 0


component_skeleton.main.main(execute)
