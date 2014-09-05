#!/usr/bin/env python

import component_skeleton.main
import sys, os
from math import isnan
from pylab import *


def rename_and_copy_WM(wm, outwm, wmname):

    wmlines = open(wm).readlines()

    #NA Logo
    wmlines[1] = 'NA %s\n' %wmname

    o = open(outwm, 'w')
    for i in wmlines:
        o.write(i)


def filterQuarterWMs(wmdict):
    """
    Sometimes WMs just contain 1 everywhere (I produce them in RunMotevo).
    This function filters those out.
    """

    badwms = []

    for wm in wmdict:
        name, mat, matrev = wm2mat(open(wm))
        if mat.shape[0] * mat.shape[1] == len(where(mat == 0.25)[0]):
            badwms.append(wm)

    for wm in badwms:
        del wmdict[wm]

    return wmdict


def wm2mat(wmFile,ps=0.5):
    """ Read swiss regulon style WM file """
    M = False
    ok = False
    lines = wmFile.read().splitlines()
    for line in lines:
        if not (line.startswith("//") or line==""):
            if ok:
                elm = reshape(array(map(float,line.split()[1:5])),(1,4))
                if M is False:
                    M = elm
                else:
                    M = vstack((M,elm))
            else:
                if line[0:2]=='NA':
                    na,name = line.split()
                elif line[0:2]=='PO' or line[0:2]=='P0': # WM start after PO/0 line
                    ok = True
                else:
                    pass
    # automatic detection if WM is already normalized
    if all(M.sum(axis=1)-1.0 < 0.0001):
        M += 0.0001 # to avoid the case that we have completely polarized column
    else: # matrix is in counts; add pseudo count
        M += 0.5

    # convert counts to frequencies
    M = M/M.sum(axis=1)[:,newaxis]
    # also build the reverse complement matrix
    Mrev = vstack((M[::-1,3],M[::-1,2],M[::-1,1],M[::-1,0])).T

    return(name,M,Mrev)


def getSimilarityScore(wm1, wm1rev, wm2, wm2rev):

    s1 = zeros((wm1.shape[0]+wm2.shape[0]-1,))
    s2 = zeros((wm1.shape[0]+wm2.shape[0]-1,))
    for n in arange(wm1.shape[1]): # over A,T,G,C                                                                                      
        s1 += convolve(wm1[:,n],wm2[::-1,n])
        s2 += convolve(wm1[:,n],wm2rev[::-1,n])
    score = vstack((s1,s2))
    idx = argmax(score,axis=0)
    idx2 = argmax(score[idx,arange(score.shape[1])])
    max_score = score[idx,arange(score.shape[1])][idx2] # idx[idx2] = 0 origianl, idx[idx2] = 1, reverse complement matches

    return max_score


def reduceWMs(wmdict, dist_co):
    """
    This function takes a dictionary with WM filepaths as keys and a score (AUCs or likelihoods) as values.
    The function goes down the sorted dictionary (by score), computes which WMs from the same dictionary 
    are similar (distance smaller than dist_co) and removes them from the dictionary.
    It does this as long as there is no WM left in the given wmdict.
    """

    final_wms = [] #list containing names of WMs to give out. WMs that have high likelihood and are not similar to any other in this list

    while len(wmdict.keys()) != 0:

        print 'Remaining candidate WMs:', len(wmdict.keys())
        #reformat wmdict to a matrix
        wmmat1 = []
        for k in wmdict:
            wmmat1.append([k, wmdict[k]])

        wmmat = sorted(wmmat1, key = lambda k: k[1], reverse=True)

        refWM = wmmat[0][0]
        print '-----------------'
        print refWM
        final_wms.append(refWM)
        ignore, rWM, rWMrev = wm2mat(open(refWM))
        queries = array(wmmat).T[0] #wmdict.keys()

        #compute score between reference WM and itself to be able to compute distance later
        rscore = getSimilarityScore(rWM, rWMrev, rWM, rWMrev)

        for wm in queries:

            ignore, qWM, qWMrev = wm2mat(open(wm))

            #compute score between query WM and itself to  be able to compute distance later
            qscore = getSimilarityScore(qWM, qWMrev, qWM, qWMrev)

            #compute score between reference WM and query WM
            rqscore = getSimilarityScore(rWM, rWMrev, qWM, qWMrev)

            # normalize score to a distance between 0 (=both WMs are the same) and 1                                                                  
            dist = 1 - 2*rqscore / (rscore + qscore)

            # compute divergence. Measures how much the longer is away from the shorter! I take this instead of symmetrical distance to detect sub matrices as sub matrices!
            rWMlen = rWM.shape[0]
            qWMlen = qWM.shape[0]
            if rWMlen < qWMlen:
                diverg = 1 - (2*rqscore) / (2*rscore)
            elif qWMlen < rWMlen:
                diverg = 1 - (2*rqscore) / (2*qscore)
            else:
                diverg = dist

            print wm, dist, diverg
            #if diverg <= dist_co:
            if dist <= dist_co:
                del wmdict[wm]
    return final_wms


def convertFloat(x):
    res = float(x)
    if isnan(res):
        raise Exception
    return res 


def execute(cf):
    """
    This component reduces candidate WMs.
    """

    ##Ports and parameters
    infile = cf.get_input("infile")
    outdir = cf.get_output("WMdir")
    log_file = cf.get_output("log_file")
    dist_co = cf.get_parameter("distance_cutoff", "float")
    minscore = cf.get_parameter("minscore", "float")
    
    os.mkdir(outdir)    
    wmdict = {} #filename: [AUC, Likelihood]

    badwms = []

    for i, line in enumerate(open(infile)):
        if i == 0:
            continue
        t = line.strip().split()
        try:
            score = convertFloat(t[1])
            if score > minscore:
                wmdict[t[0]] = score
            else:
                print "Motif %s wasn't included to the initial list (too low score)" % t[0]
        except:
            print "Motif %s wasn't included to the initial list (score probably nan)" % t[0]
            continue 

    wmdict = filterQuarterWMs(wmdict)

    final_wms = reduceWMs(wmdict, dist_co)
    print final_wms

    l = open(log_file, 'w')

    l.write('Passed WMs:\n')
    for i, wm in enumerate(final_wms):
        os.system('cp \'%s\' %s' %(wm, outdir))
        #rename_and_copy_WM(wm, '%s/WM_%i' %(outdir, i+1), 'WM_%i' %(i+1))
        l.write('WM_%i\t%s\n' %(i+1, wm))

    l.close()


    return 0

component_skeleton.main.main(execute)
