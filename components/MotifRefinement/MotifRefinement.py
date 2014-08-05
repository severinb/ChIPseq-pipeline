#!/usr/bin/env python
import component_skeleton.main
from motevo_stuff import run_motevo
import os

def execute(cf):
    """
    This component receives a weight matrix and set of sequences (non aligned)
    and by running MotEvo in WMREF mode, refines the input motif to explain the
    sequence data better.
    As output, it generates a new version (refined) motif. 
    """
    inputSequences = cf.get_input("InputSequences")
    wmFile = cf.get_input("WM")
    output_file = cf.get_output("RefinedMotif")
    genome = cf.get_parameter("genome", "string")
    output_dir = os.path.dirname(output_file)
    (siteFilename, priorFilename, paramFilename) = run_motevo(wmFile, inputSequences, \
                                               output_file, output_dir, genome)
    os.system( 'rm %s %s %s' % (siteFilename, priorFilename, paramFilename) )
    

component_skeleton.main.main(execute)
