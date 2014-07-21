#!/usr/bin/env python
import component_skeleton.main
from motevo_stuff import *

def execute(cf):
    """
    This component is responsible to fit the background prior for a motif by using MotEvo,
    and as well fitting the Beta which is going to be used to calculate the enrichment score
    for the motif. 
    """
    interm_dir = cf.get_output("intermediate")    
    inputSequences = cf.get_input("InputSequences")
    wmFile = cf.get_input("WM")
    genome = cf.get_parameter("genome", "string")
    motevo_path = cf.get_parameter("motevo_path", "string")
    mylogo_path = cf.get_parameter("mylogo_path", "string")
    os.mkdir(interm_dir)
    (siteFilename, priorFilename) = run_motevo(motevo_path, wmFile, inputSequences, interm_dir, gemome)
    beta = fit_beta(siteFilename, interm_dir, wmFile)
    os.system('rm %s' % siteFilename)
           
component_skeleton.main.main(execute)
