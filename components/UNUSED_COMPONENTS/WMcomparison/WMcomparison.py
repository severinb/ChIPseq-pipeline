#!/usr/bin/env python
import component_skeleton.main
import subprocess
import os, re
from string import *
import datetime
from pylab import *

def normalizeWM(wm):
    """
    This function reads in a WM file to a array/matrix.
    This function normalizes WM columns to sum up to 1.
    """

    cols = open(wm).readlines()
    
    h0 = cols[0]
    h1 = cols[1]
    h2 = cols[2]
    t = cols[-1]

    cols = cols[3:-1]

    #make a matrix
    m = 4
    n = len(cols)
    WMat = zeros(m*n).reshape(4,n)

    for i, c in enumerate(cols):
        a1 = c.split()
        a = array(a1[1:5], dtype='float')
        try:
            num = int(a1[0])
        except Exception, e:
            print '\n\n\n\n'
            print wm
            print e

        aN = a/sum(a)

        for j in arange(len(aN)): 
            WMat.T[i][j] = aN[j]

    return WMat    


def complementWM(WM):
    """
    This function just swaps entries of each column between A/T and G/C 
    WM:
      1 2 3 4 5 6 7 8 9 10
    A
    C
    G
    T
    Additionally it reverses the matrix, i.e. first becomes last column and so on
    """

    cWM = array(WM)

    for i in arange(len(WM.T)):
        col = list(WM.T[i])
        cWM.T[-i-1][0] = col[3] #put T val to A spot
        cWM.T[-i-1][3] = col[0] #put A val to T spot
        cWM.T[-i-1][1] = col[2] #put G val to C spot
        cWM.T[-i-1][2] = col[1] #put C val to G spot

    return cWM


def writeWM(WM, interm, name):

    outfile = os.path.join(interm, name)
    o = open(outfile, 'w')
    o.write('//\n')
    o.write('NA\t%s\n' %name)
    o.write('P0\tA\tC\tG\tT\n')

    for pos, c in enumerate(WM.T):
        o.write('\t'.join([str(pos+1)] + map(str, list(c)) + ['\n']))

    o.write('//\n')

    return outfile

def SSD(cA, cB):
    """
    This function computes the sum of squared distances between two columns.
    2 is the maximum score. e.g. [1,0,0,0] and [0,1,0,0]
    """

    return 2 - sum((cA - cB)**2)


def TraceBack(Smat, Pmat):
    """
    this function looks for the max value in Smat. Starting from there and with the help of Pmat it traces back the alignment
    """

    m = Smat.shape[0]
    n = Smat.shape[1]

    maxVal = 0
    posi = 0
    posj = 0
    for i in arange(m):
        for j in arange(n):
            if Smat[i][j] > maxVal:
                maxVal = Smat[i][j]
                posi = i
                posj = j

    print 'maximum Alignment Score: ', maxVal
    print posi, posj

    seq1 = [] #alignment sequence of WM1. This WM is the first row
    seq2 = [] #This WM2 is the first column
    AlignScore = 0

    caution = 0

    while 1:
        caution += 1
        if caution > 100:
            break
        if posi <= 0 and posj <= 0:
            break

        AlignScore += Smat[posi][posj]
        if Pmat[posi][posj] == 0: #diag
            seq1 = ['x'] + seq1
            seq2 = ['x'] + seq2
            posi -= 1
            posj -= 1
        elif Pmat[posi][posj] == 1: #up
            seq1 = ['x'] + seq1
            seq2 = ['-'] + seq2
            posi -= 1
        elif Pmat[posi][posj] == 2: #left
            seq1 = ['-'] + seq1
            seq2 = ['x'] + seq2
            posj -= 1

        print posi, posj


    print 'WM1 ', ''.join(seq1)
    print 'WM2 ', ''.join(seq2)

    L = len(seq1)
    return (maxVal)/L, ' | '.join([seq1[i] + '' + seq2[i] for i in arange(L)]) #normalize the maxVal to account for longer alignments (because scores are always positive). But maybe this is not a good normalization.


def trimAln(aln, score):
    """
    this function takes an alignments string, splits it and counts gaps on the edges of the alignment. 
    If there are gaps, they're cut away and gap scores are also subtracted to not dilute the bp normalized alignment score.
    """

    bps = aln.split(' | ')

    rawScore = score * len(bps)

    #trim start
    while 1:  
        t = list(bps[0])
        if t[0] == '-' or t[1] == '-':
            del bps[0]
            rawScore -= gappen
            continue
        else:
            break

    
    #trim end
    if not len(bps) == 0:        
        while 1:  
            t = list(bps[-1])
            if t[0] == '-' or t[1] == '-':
                del bps[-1]
                rawScore -= gappen
                continue
            else:
                break


    return (rawScore)/len(bps), ' | '.join(bps)


def SWalign(wm1, wm2):
    """
    This function executes smith waterman alignment.
    For wm1 of length 10 and wm2 of length 6:
     X - c1 c2 c3 c4 c5 c6 c7 c8 c9 c10
     - 0 g  g  g  g  g  g  g  g  g  g
    c1 g
    c2 g
    c3 g
    c4 g
    c5 g
    c6 g

    g= gap penalty
    Options are: match, gap, indel
    Alignment of cmi with cnj A(i,j) = max( A(i-1,j-1) + S(i,j) ,  A(i-1,j) + g , A(i,j-1) + g  )
    g = gap penalty (e.g. 0) (#gap penalty. Can't be negative, since scores also can't be negative)
    A(i-1,j) + g means cm(i-1) is aligned with cn(j) and cm(i) is paired with a gap (so gap in cn)
    A(i,j-1) + g means cm(i) is aligned with cn(j-1) and cn(j) is paired with a gap (so gap in cm)
    S(i,j) is the SSD of columns i and j
    """

    m = wm1.shape[1] + 1 #+1 because of adding the gap penalty row and column for initial gap aligning
    n = wm2.shape[1] + 1

    Smat = zeros(m*n).reshape(m,n)  #matrix that contains scores
    Pmat = zeros(m*n).reshape(m,n)  #matrix that contains pointers: 0 means diagonal (aligned), 1 means up (gap in cn), 2 means left (gap in cm)

    print 'non-affine Gap Penalty: ', gappen

    for i in arange(m):
        for j in arange(n):
            if i==0:
                Smat[i][j] = gappen
                Pmat[i][j] = 2 #2 means up, gap in cm 
            elif j==0:
                Smat[i][j] = gappen
                Pmat[i][j] = 1 #1 means left, gap in cn
            else:
                vals = [ Smat[i-1][j-1] + SSD(wm1.T[i-1], wm2.T[j-1])  ,  Smat[i-1][j] + gappen  ,  Smat[i][j-1] + gappen  ]  #wm1.T[i-1] subtract 1 because WMs have one column and row less than Smat and Pmat 
                a = argmax(vals)
                Smat[i][j] = round(vals[a],3)
                Pmat[i][j] = a
    print Smat
    print Pmat
    
    score, aln = TraceBack(Smat, Pmat)

    scoreT, alnT = trimAln(aln, score)

    return scoreT, aln




def DrawLogo(WM, interm, mylogopath):
    """
    This function draws the sequence logo using weblogo. 
    """

    #set the output image inside the wm file
    wmt = open(WM).readlines()
    tempWMfile = os.path.join(interm, os.path.split(WM)[1])
    wmt[1] = 'NA %s\n' %tempWMfile
    wmo = open(tempWMfile, 'w')
    wmo.write(''.join(wmt))
    wmo.close()

    proc = subprocess.Popen('%s -n -a -c -p -Y -F PDF -f %s' %(mylogopath, tempWMfile),
                             stdout=subprocess.PIPE,
                             stderr= subprocess.PIPE,
                             shell=True
                            )

    stdout_value, stderr_value = proc.communicate()
    print stdout_value
    print stderr_value


def AlignMain(wm1,wm2):
    """
    Smith waterman alignment of the WMs with sum of squared distances scoring.

    Note: Also cares about complement. Maybe test all three possibilities: input-input, input-complement, complement-input. I think the last two are the same
    """

    WM1 = normalizeWM(wm1)
    WM2 = normalizeWM(wm2)
    cWM2 = complementWM(WM2)

    print '\n-----------------------\ninput-input\n'
    score1, aln1 = SWalign(WM1, WM2)
    print '\n-----------------------\ninput-comp\n'
    score2, aln2 = SWalign(WM1, cWM2)

    a = argmax([score1,score2])  #variable a: 0 means alignment to non-complement is better than to complement. 1 is the other way round. 
    best = [score1,score2][a] 
    print best,a

    #writeWM(cWM1, 'compJun')
    #writeWM(WM1, 'Junnnn')

    return best, a, [aln1,aln2][a]



def execute(cf):
    WM = cf.get_input("WM") #result file from AlignPeaks spltting from SplitAlignments (typically even half)
    WM_dir = cf.get_input("WM_dir") 

    out_file = cf.get_output("out_file")
    interm = cf.get_output("intermediate")
    logfile = cf.get_output("log_file") 

    mylogo_path = cf.get_parameter("mylogo_path", "string")
    numWMs = cf.get_parameter("numWMs", "int")
    global gappen
    gappen = cf.get_parameter("gap_penalty", "float")


    T1 = datetime.datetime.now()

    os.mkdir(interm)

    o = open(out_file,'w')
    scorelist = []
    for i in os.listdir(WM_dir):
        wmi = os.path.join(WM_dir,i)
        s, a, aln = AlignMain(WM,wmi)
        scorelist.append((wmi,s,a))

        o.write(wmi + '\t%s\t%s\t%s\n' %(s, a, aln))

    o.close()

    sortedscorelist = sorted(scorelist, key= lambda k: k[1], reverse=1)[:5]
    print sortedscorelist
    DrawLogo(WM, interm, mylogo_path)

    for i in sortedscorelist:
        wm = i[0]
        a = i[-1]
        if a == 1:
            wmpath = writeWM(complementWM(normalizeWM(wm)), interm, os.path.split(wm)[1])
            DrawLogo(wmpath, interm, mylogo_path)
        else:
            DrawLogo(wm, interm, mylogo_path)



    T2 = datetime.datetime.now()

    lf = open(logfile, 'w')
    lf.close
                        

    return 0

component_skeleton.main.main(execute)
                                                                 
