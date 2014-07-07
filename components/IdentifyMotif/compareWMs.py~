#!/usr/bin/env python

""" Compute the similarity between a given WM and a WM data base (directory of WM files) """

import sys,os,argparse,textwrap
import numpy as np
from subprocess import call

def wm2mat(wmFile,ps=0.5):
    """ Read swiss regulon style WM file """
    M = False
    ok = False
    lines = wmFile.read().splitlines()
    for line in lines:
        if not (line.startswith("//") or line==""):
            if ok:
                elm = np.reshape(np.array(map(float,line.split()[1:5])),(1,4))
                if M is False:
                    M = elm
                else:
                    M = np.vstack((M,elm))
            else:
                if line[0:2]=='NA':
                    na,name = line.split()
                elif line[0:2]=='PO' or line[0:2]=='P0': # WM start after PO/0 line 
                    ok = True
                else:
                    pass
    # add pseudo count
    M += ps
    # convert counts to frequencies
    M = M/M.sum(axis=1)[:,np.newaxis]
    # also build the reverse complement matrix
    Mrev = np.vstack((M[::-1,3],M[::-1,2],M[::-1,1],M[::-1,0])).T

    return(name,M,Mrev)

def scoreWMs(qWM,rWM,rWMrev):
    # distance: matrix convolution
    s1 = np.zeros((qWM.shape[0]+rWM.shape[0]-1,))
    s2 = np.zeros((qWM.shape[0]+rWM.shape[0]-1,))
    for n in np.arange(qWM.shape[1]): # over A,T,G,C
        s1 += np.convolve(qWM[:,n],rWM[::-1,n])
        s2 += np.convolve(qWM[:,n],rWMrev[::-1,n])
    score = np.vstack((s1,s2))
    idx = np.argmax(score,axis=0)
    idx2 = np.argmax(score[idx,np.arange(score.shape[1])])
    max_score = score[idx,np.arange(score.shape[1])][idx2]
    rev_comp = idx[idx2] # idx[idx2] = 0 original, idx[idx2] = 1, reverse complement matches
    return(max_score,rev_comp)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description=textwrap.dedent(__doc__))
    parser.add_argument('-wmdir', dest='wmdir', action='store', required=True, help='Directory containing reference WMs. If [WM] is omitted calculates the similarity between all pairs of WMs in [WMDIR].')
    parser.add_argument('-wm', dest='wm', action='store', default=None, help='Query WM. Calulates similarity between query WM and all reference WMs in [WMDIR].')
    parser.add_argument('-norm', dest='norm', action='store_true', default=False, help='Normalize scores to distances: D_ij = 1 - 2*S_ij/(S_ii+S_jj)')
    parser.add_argument('-ntop', dest='ntop', action='store', type=int, default=0, help='Only print a sorted list top N matches. Only works if query is single WM!')

    args = parser.parse_args()

    # reference WMs
    refWM = {}
    for wm in os.listdir(args.wmdir):
        WMname,WM,WMrev = wm2mat(open(os.path.join(args.wmdir,wm)))
        refWM[WMname] = (WM,WMrev)
        
    # query WM
    queryWM = {}
    if not args.wm is None:
        qWMname,qWM,qWMrev = wm2mat(open(args.wm))
        qq_max_score,qq_rev_com = scoreWMs(qWM,qWM,qWMrev)
        S = np.zeros((len(refWM),2))
        for r,rwm in enumerate(sorted(refWM.keys())):
            key = str(sorted((qWM,rwm)))
            rWM,rWMrev = refWM[rwm]
            # distance: matrix convolution
            qr_max_score,rev_com = scoreWMs(qWM,rWM,rWMrev)
            S[r,0] = qr_max_score
            rr_max_score,rev_com = scoreWMs(rWM,rWM,rWMrev)
            S[r,1] = rr_max_score
    else:
        queryWM = refWM # all against all  
        scores = {}
        S = np.zeros((len(queryWM),len(refWM)))
        for q,qwm in enumerate(sorted(queryWM.keys())):
            qWM,qWMrev = queryWM[qwm]
            for r,rwm in enumerate(sorted(refWM.keys())):
                key = str(sorted((qwm,rwm)))
                rWM,rWMrev = refWM[rwm]
                if not scores.has_key(key):
                    # distance: matrix convolution
                    max_score,rev_comp = scoreWMs(qWM,rWM,rWMrev)
                    scores[key] = (max_score,wm, rev_comp) 
                    S[q,r] = max_score
                    if S.shape[0] == S.shape[1]: # all against all query 
                        S[r,q] = S[q,r]
                else:
                    pass
    
    # normalize score to a distance between 0 (=both WMs are the same) and 1
    if args.norm:
        if args.wm is None: # all against all query 
            S = 1.0 - 2.0*S / np.transpose(np.tile(S.diagonal(),(S.shape[0],1)).T + np.tile(S.diagonal(),(S.shape[0],1)))
        else:
            S = 1.0 - 2.0*S[:,0]/(S[:,1]+qq_max_score)

    if args.ntop > 0 and not args.wm is None: # single query WM
        if args.norm:
            top = sorted(zip(S.flatten(),sorted(refWM.keys())),reverse=False)
        else:
            top = sorted(zip(S.flatten(),sorted(refWM.keys())),reverse=True)
        for t in range(min(args.ntop,S.shape[0])):
            print "%s\t%g"%(os.path.basename(top[t][1]),top[t][0])
    else:        
        # print output in R matrix form: read in with: d = read.table('output.dat')
        print '\t' + '\t'.join("%s" %os.path.basename(s) for s in sorted(refWM.keys()))
        for q,qwm in enumerate(sorted(queryWM.keys())):
            row = os.path.basename(qwm) + '\t'
            for r,rwm in enumerate(sorted(refWM.keys())):
                row += "%g\t" % S[q,r]
            print row.strip()
