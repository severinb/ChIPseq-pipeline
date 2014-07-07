#!/usr/bin/env python

""" Scan the genome with all WMs in WM_DIR. Input bed file is sorted and from one chromosome only. """

import sys,os
import tables
import gzip
import argparse
import numpy as np

# A, C, G, T, N
#BG_DEFAULT = np.array([0.25, 0.25, 0.25, 0.25, 0.25],dtype=float)
BG_DEFAULT = '0.25,0.25,0.25,0.25,0.25'

#
NUM2SEQ = np.array(['A','C','G','T','N'],dtype='S1')
STRAND = ['+','-']

def wm2mat(wmFile,BG,ps=0.5):
    """ Read swiss regulon style WM file """
    M = False
    ok = False
    lines = wmFile.read().splitlines()
    for line in lines:
        if not (line[0:2]=="//" or line==""):
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

    # add an 'N' column and set it to the mean over A,T,G,C freq, i.e. to 0.25
    M = np.hstack((M,M.mean(axis=1)[:,np.newaxis]))

    # also build the reverse complement matrix
    Mrev = np.vstack((M[::-1,3],M[::-1,2],M[::-1,1],M[::-1,0],M[::-1,4])).T

    # now build log odds
    M = np.log(M/BG)
    Mrev = np.log(Mrev/BG)
    
    return(name,M,Mrev)        
    
if __name__ == '__main__':

    mat = sys.argv[1]
    bed = sys.argv[2]
    out = sys.argv[3]
    cutoff = float(sys.argv[4])
    chrDir = sys.argv[5]
    bg = sys.argv[6]

    jobid = int(os.environ['SGE_TASK_ID'])

    try:
        bedFile = open(bed,'r')    
    except:
        print "Can't read bed input file!"
        sys.exit(1)

    try:
        outFile = open(os.path.join(out, 'outfile_' + str(jobid)), 'w')    
        print outFile
    except:
        print "Can't open output bed file!"
        sys.exit(1)
    
    # genomic background frequencies
    try:
        BG = np.fromstring(bg,dtype=float,sep=',')
    except:
        print "Can't read genomic background frequencies!"
        sys.exit(1)
    if not len(BG) == 5:
        print "Wrong size of genomic background frequencies! Required format: A,C,G,T,N"
        sys.exit(1)

    # weight matrix
    try:
        matFile = open(mat,'r')
        WMname,WM,WMrev = wm2mat(matFile,BG)
        WMlen = WM.shape[0]
        rowidx = np.arange(WM.shape[0]) # 0...len(WM)
    except:
        print "Can't read matrix file!"
        sys.exit(1)

        
    # open table files
    if not os.path.exists(chrDir):
        print "Can't read genome directory!"
        sys.exit()
        
    lines =  bedFile.readlines()


    #select e a specific line from chrominfo file (according to jobid) 
    line = lines[jobid-1]
    elm = line.split()
    print elm
    chrom = elm[0]
    tableFile = tables.openFile(os.path.join(chrDir,chrom+'.h5'),'r')

    start = int(elm[1]) #1
    stop = int(elm[2]) #int(elm[1])
    nseq = tableFile.root.numseq[start:stop]
    for p in xrange(len(nseq)-WMlen):
        colidx = nseq[p:(p+WMlen)]
        score = np.hstack((np.sum(WM[rowidx,colidx]),np.sum(WMrev[rowidx,colidx]))) # forward and reverse score
        if np.any(score >= cutoff):
            maxidx = score.argmax()
            seq = NUM2SEQ[colidx].tostring()
            outFile.write("%s\t%d\t%d\t%s\t%.8f\t%s\n" %(chrom,start+p,start+p+WMlen,seq,score[maxidx],STRAND[maxidx]))

    outFile.close()
    
