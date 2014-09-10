#!/usr/bin/env python
import component_skeleton.main
from motevo_stuff import run_motevo
from sequence_logo import generate_sequence_logo
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
    output_file = cf.get_output("refWM")
    output_logo = cf.get_output("Logo")    
    genome = cf.get_parameter("genome", "string")
    weblogo_path = cf.get_parameter("weblogo_path", "string")
    output_dir = os.path.dirname(output_file)
    (siteFilename, priorFilename, paramFilename) = run_motevo(wmFile, inputSequences, \
                                               output_file, output_dir, genome)
    desc = ' ' #'%s refined' % os.path.basename(os.path.dirname(output_file))    
    generate_sequence_logo(output_file, output_logo, weblogo_path, desc)
    # os.system( 'rm %s %s %s' % (siteFilename, priorFilename, paramFilename) )
    return 0

component_skeleton.main.main(execute)
