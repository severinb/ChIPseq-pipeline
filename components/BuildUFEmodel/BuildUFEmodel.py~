#!/usr/bin/env python
import component_skeleton.main
import subprocess
import os, re
from string import *
import datetime
from pylab import *

def getBaseFreqs(f):
    """
    This function computes base frequencies for input sequences
    """

    bd = {}
    bd['A'] = 0.
    bd['C'] = 0.
    bd['G'] = 0.
    bd['T'] = 0.

    for line in open(f):
        if line.startswith('>'):
            continue
        for b in bd:
            bases = list(line.strip())
            bd[b] += bases.count(b)

    print bd
    #normalize to frequencies
    tot = sum(bd.values())
    for i in bd:
        bd[i] /= tot

    print bd
    #take same frequencies for A and T or C and G respectively, because we are using double stranded DNA at the end.
    prec = 3
    ATfreq = round((bd['A']+bd['T'])/2., prec)
    GCfreq = round((bd['C']+bd['G'])/2., prec)

    #round frequencies to some precision and be sure that they sum up to exactly 1
    #round one frequency (AT or GC) down and the other up
    while (ATfreq + GCfreq) != 0.5:
        if (ATfreq + GCfreq) > 0.5:
            if rand(1) >= 0.5:
                ATfreq -= 1.0/(10**prec)
            else:
                GCfreq -= 1.0/(10**prec)
        else:
            if rand(1) >= 0.5:
                ATfreq += 1.0/(10**prec)
            else:
                GCfreq += 1.0/(10**prec)

    return ATfreq, GCfreq


    
def runMotevo(motevo_path, treefile, outfile, ATfreq, GCfreq):
    """
    runs Motevo
    """

    print 'A/T frequency:', ATfreq
    print 'G/C frequency:', GCfreq

    proc = subprocess.Popen(motevo_path + ' %s %s %s %s %s > %s' %(treefile, ATfreq, GCfreq, GCfreq, ATfreq, outfile),
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
    This function/component refines a given weight matrix, then predicts sites for old and refined WM and true sites and on background peaks.
    With this it makes a positive predictive value - sensitivity plot.
    It also makes a logo of the refined WM.
    """

    ##Ports and parameters
    seqs = cf.get_input("Sequences") #result file from AlignPeaks spltting from SplitAlignments (typically even half)
    outfile = cf.get_output("UFEmodel")
    basefreqs = cf.get_output("BaseFrequencies")
    genome = cf.get_parameter("genome", "string")
    motevo_path = cf.get_parameter("motevoUFE_path", "string")


    genome_dict = {}
    genome_dict['hg19'] = '((((hg19:0.032973,rheMac2:0.057695):0.09821,mm9:0.352605):0.020666,(bosTau6:0.186713,(equCab2:0.107726,canFam2:0.150374):0.010431):0.032764):0.156024,monDom5:0.425899);'
    genome_dict['hg18'] = '((((hg18:0.032973,rheMac2:0.057695):0.09821,mm9:0.352605):0.020666,(bosTau3:0.186713,(equCab1:0.107726,canFam2:0.150374):0.010431):0.032764):0.156024,monDom4:0.425899);'
    genome_dict['mm9'] = '((((hg19:0.032973,rheMac2:0.057695):0.09821,mm9:0.352605):0.020666,(bosTau7:0.186713,(equCab2:0.107726,canFam2:0.150374):0.010431):0.032764):0.156024,monDom5:0.425899);'

    if genome == 'dm3':
        dm3model = '/import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/MotEvo_v1.0/UFEmodels/dm3UFEparallel/UFEmodel_dm3'
        print 'Warning!! UFE model for dm3 can not be built because of too many species in phylogenetic tree.'
        print 'Using standard UFE model instead (A/T:): %s'%dm3model
        ATfreq = 0.295
        GCfreq = 0.205
        o = open(basefreqs, 'w')
        o.write('%s\t%s\n' %('AT', ATfreq))
        o.write('%s\t%s\n' %('GC', GCfreq))
        o.close()
        os.system('cp %s %s' %(dm3model, outfile))

        return 0

    #genome_dict['dm3'] = ['((((((dm3:0.059,droSim1:0.075):0.041,(droYak2:0.104,droEre2:0.107):0.054):0.120,droAna3:0.377):0.072,dp4:0.397):0.061,droWil1:0.536):0.020,((droVir3:0.196,droMoj3:0.255):0.073,droGri2:0.291):0.337);', '/import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/MotEvo_v1.0/UFEmodels/dm3UFEparallel/UFEmodel_dm3']


    ATfreq, GCfreq = getBaseFreqs(seqs)

    o = open(basefreqs, 'w')
    o.write('%s\t%s\n' %('AT', ATfreq))
    o.write('%s\t%s\n' %('GC', GCfreq))
    o.close()


    treefile = os.path.join(os.path.split(outfile)[0], 'treefile')
    o = open(treefile, 'w')
    o.write(genome_dict[genome])
    o.close()

    #run motevo to build UFE model
    rval = runMotevo(motevo_path, treefile, outfile, ATfreq, GCfreq)


    return rval


component_skeleton.main.main(execute)
                                                                 
