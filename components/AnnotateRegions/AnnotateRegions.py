#!/usr/bin/env python
import component_skeleton.main
import subprocess
import os
from string import *
import datetime, time
from pylab import *

def produceHomerBed(infile, newfile):

    o = open(newfile, 'w')
    ids = []

    for line in open(infile):
        t = line.strip().split()
        chrom = t[0]
        start = t[1]
        end = t[2]
        strand = t[-1]
        ID = t[3]
        out = '\t'.join([chrom, start, end, ID, '', strand])+'\n'
        o.write(out)

    o.close()

    return newfile


def execute(cf):
    """
    This component runs homer's (http://biowhat.ucsd.edu/homer/index.html) annotatePeaks.pl to annotate each given peak with:
        -Nearest TSS:
            -distance to nearest TSS from RefSeq (negative means upstream, positive means downstream)
            -ID of nearest promoter (RefSeq)
            -ID of associated gene (RefSeq, Enseml, Entrez, Unigene)
            -Gene name, alias, description and type
        -Genomic Annotation. Does peak lie in:
            -TSS (-1kb +100bp)
            -TTS (-100bp +1kb)
            -CDS exons
            -5' UTR
            -3' UTR
            -CpG islands
            -repeats
            -Introns
            -Intergenic
        -Gene ontology analysis of nearest genes: enrichment for biological functions of genes.
        -Genome ontology analysis of peaks: genomic regions that are enriched in the input peaks.
    """

    ##Ports and parameters
    regions = cf.get_input("regions") #result file from AlignPeaks spltting from SplitAlignments (typically even half)
    interm = cf.get_output("intermediate")
    genomeOnt = cf.get_output("GenomeOntology")
    geneOnt = cf.get_output("GeneOntology")
    outfile = cf.get_output("regionAnnotations")
    genome = cf.get_parameter("genome", "string")
    homerPATH = cf.get_parameter("homerPATH", "string")
    #perlPATH = cf.get_parameter("perlPATH", "string")

    os.mkdir(interm)

    regions = produceHomerBed(regions, os.path.join(interm, 'regions.homer'))

    #add homer's directory to the PATH
    os.environ['PATH'] = homerPATH + ':' + os.environ['PATH']    

    command = ' '.join(['annotatePeaks.pl',
                        regions, 
                        genome,
                        #'/import/bc2/data/databases/UCSC/hg19',
                        '-go', geneOnt,
                        '-genomeOntology', genomeOnt,
                        '>', outfile])

    print '\nrun homer on %s' %regions
    proc = subprocess.Popen(command,
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


component_skeleton.main.main(execute)
                                                                 
