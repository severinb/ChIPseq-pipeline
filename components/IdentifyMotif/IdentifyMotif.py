#!/usr/bin/env python
import component_skeleton.main
import os, re
from string import *
from pylab import *
import subprocess
from sequence_logo import generate_sequence_logo

def runComparison(WM, WMlib, outfile, ntop):
    """
    runs Florian's compareWMs.py script.
    """

    proc = subprocess.Popen('./compareWMs.py -wm \'%s\' -wmdir %s -norm -ntop %i > %s' %(WM, WMlib, ntop, outfile),
                            stdout=subprocess.PIPE,
                            stderr= subprocess.PIPE,
                            shell=True
                            )

    stdout_value, stderr_value = proc.communicate()
    print stdout_value
    print stderr_value

    if proc.poll() > 0:
        print '\tstderr:', repr(stderr_value.rstrip())
        return -1
    else:
        return 0


def execute(cf):
    """
    This function runs a greedy algorithm to compute the log-likelihood of sequences given the best combination of input WMs.
    """

    ##Ports and parameters
    WM = cf.get_input("WM")

    outfile = cf.get_output("top_motifs")
    logos_dir = cf.get_output("logos")

    WMlib = cf.get_parameter("WMlibrary", "string")
    ntop = cf.get_parameter("ntop", "int")
    mylogo_path = cf.get_parameter("mylogo_path", "string")

    os.mkdir(logos_dir)

    tmpfile = os.path.join(os.path.split(outfile)[0], 'tmp')
    runComparison(WM, WMlib, tmpfile, ntop)

    o = open(outfile, 'w')
    for l in open(tmpfile):
        t = l.strip().split()

        wm_path = os.path.join(WMlib, t[0])
        logo_path = os.path.join(logos_dir, t[0])
        generate_sequence_logo(wmpath, logo_path)

        dist_percent = float(t[1])*100
        o.write('%s\t%.3f\n' %(t[0], float(t[1])))

    o.close()

    return 0


component_skeleton.main.main(execute)
                                                                 
