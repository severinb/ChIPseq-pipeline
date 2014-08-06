#!/usr/bin/env python

import component_skeleton.main
from motevo_stuff import *
from fitting_beta import *
import os, yaml 

def create_training_pool(destDir, realSeq, decoySeq):
    trainingPool = os.path.join(destDir, 'training_pool')    
    cmd = ' '.join([
        'cat',
        realSeq,
        decoySeq,
        '>',
        trainingPool
        ])
    os.system(cmd)
    return trainingPool


def extract_priors(prior_file):
    prior = [float(l.split()[1]) for l in open(prior_file) if re.search('^background', l)].pop()
    return prior


def generate_sequence_logo(motifFile, mylogo_path):
    cmd = ' '.join([
        mylogo_path,
        '-n -a -c -p -Y -F PDF',
        '-f %s' % motifFile,
        ])
    os.system(cmd)
    return '%s.pdf' % motifFile

    
def execute(cf):
    """
    This component is responsible to fit the background prior for a motif by using MotEvo,
    and as well fitting the Beta which is going to be used to calculate the enrichment score
    for the motif. 
    """
    inputSequences = cf.get_input("InputSequences")
    decoySequences = cf.get_input("DecoySequences")
    wmFile = cf.get_input("WM")    
    genome = cf.get_parameter("genome", "string")
    motevo_path = cf.get_parameter("motevo_path", "string")
    out_file = cf.get_output("FittedParameters")
    interm_dir = os.path.join(os.path.dirname(out_file), "intermediate")
    os.system('mkdir %s' % interm_dir)
    trainingPool = create_training_pool(interm_dir, inputSequences, decoySequences)
    (siteFilename, priorFilename) = run_motevo(motevo_path, wmFile, trainingPool, interm_dir, genome)
    prior = extract_priors(priorFilename)
    beta = fit_beta(siteFilename, interm_dir, wmFile)
    # writing the fitted parameters
    output_file = open(out_file, 'w')  
    yaml.dump({'prior':prior, 'beta':round(beta, 10)}, output_file) # we use YAML format for saving the parameters
    output_file.close()
    # cleaning up a bit! 
    # os.system('rm %s' % siteFilename) # we won't need sites file that is generated by MotEvo
    os.system('rm %s' % trainingPool)
    os.system('gzip -r %s ' % interm_dir ) # compressing the intermediate directory
    return 0
    
component_skeleton.main.main(execute)
