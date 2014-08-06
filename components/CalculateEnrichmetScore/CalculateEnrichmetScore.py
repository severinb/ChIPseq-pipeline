#!/usr/bin/env python

import component_skeleton.main

def create_test_pool(destDir, realSeq, decoySeq):
    testPool = os.path.join(destDir, 'tesg_pool')    
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
    
    fittedParams = cf.get_input("FittedParameters")    
    genome = cf.get_parameter("genome", "string")
    motevo_path = cf.get_parameter("motevo_path", "string")
    
    # out_file = cf.get_output("")
    return 0

component_skeleton.main.main(execute)
