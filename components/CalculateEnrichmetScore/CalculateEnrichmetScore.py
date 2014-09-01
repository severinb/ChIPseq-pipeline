#!/usr/bin/env python

import component_skeleton.main
import yaml, os
from motevo_stuff import run_motevo
from enrichment_score import calculate_enrichment_scores 

def create_test_pool(destDir, realSeq, decoySeq):
    testPool = os.path.join(destDir, 'test_pool')    
    cmd = ' '.join([
        'cat',
        realSeq,
        decoySeq,
        '>',
        testPool
        ])
    os.system(cmd)
    return testPool


def execute(cf):
    """
    Calculates the enrichment score for a motif over the test pool.
    This component receives beta (already fitted by another component), and background
    prior. In addition to that it gets a set of sequences and for this set of sequences
    it calculates the enrichmet score. 
    """
    inputSequences = cf.get_input("InputSequences")
    decoySequences = cf.get_input("DecoySequences")
    wmFile = cf.get_input("WM")
    fittedParams = cf.get_input("FittedParams")    
    genome = cf.get_parameter("genome", "string")
    motevo_path = cf.get_parameter("motevo_path", "string")
    enrichment_file = cf.get_output("EnrichmentScores")
    interm_dir = os.path.join(os.path.dirname(enrichment_file), "intermediate")    
    os.system('mkdir %s' % interm_dir)
    testPool = create_test_pool(interm_dir, inputSequences, decoySequences)
    with open(fittedParams) as f:
        params = yaml.load( f )
    siteFilename = run_motevo(motevo_path, wmFile, \
                                               params['prior'], testPool, interm_dir, genome)
    calculate_enrichment_scores(siteFilename, params['beta'], enrichment_file)    
    return 0

component_skeleton.main.main(execute)
