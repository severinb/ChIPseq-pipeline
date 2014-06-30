#!/usr/bin/env python
import component_skeleton.main
import os, re
from string import *
from pylab import *
import subprocess


def runComparison(WM, WMlib, outfile, ntop):
    """
    runs Florian's compareWMs.py script.
    """

    proc = subprocess.Popen('./compareWMs.py -wm %s -wmdir %s -norm -ntop %i > %s' %(WM, WMlib, ntop, outfile),
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

    WMlib = cf.get_parameter("WMlibrary", "string")
    ntop = cf.get_parameter("ntop", "int")


    runComparison(WM, WMlib, outfile, ntop)


    return 0


component_skeleton.main.main(execute)
                                                                 
