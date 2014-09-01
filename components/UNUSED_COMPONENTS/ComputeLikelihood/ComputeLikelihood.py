#!/usr/bin/env python
import component_skeleton.main
import os, re
from string import *
import datetime, time
from pylab import *

def getBaseFreqs(f1,f2):
    """
    This function computes base frequencies for input sequences
    """

    bd = {}
    bd['A'] = 0.
    bd['C'] = 0.
    bd['G'] = 0.
    bd['T'] = 0.
    bd['N'] = 0.

    for f in [f1, f2]:
        for line in open(f1):
            if line.startswith('>'):
                continue
            for b in bd:
                bases = list(line.strip())
                bd[b] += bases.count(b)

    #normalize to frequencies
    tot = sum(bd.values())
    for i in bd:
        bd[i] /= tot


    #put A,T and G,C to the same frequencies
    ATfreq = (bd['A'] + bd['T'])/2.0
    GCfreq = (bd['C'] + bd['G'])/2.0

    bd2 = {}
    bd2['A'] = ATfreq
    bd2['C'] = ATfreq
    bd2['G'] = GCfreq
    bd2['T'] = GCfreq
    bd2['N'] = bd['N']

    return bd2


def wm2mat(wmFile,ps=0.5):
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


    return M, Mrev


def computeLikelihood(seqs, M, Mrev, prior, basedict):

    #A:0, C:1, G:2, T:3
    seq2num = {}
    seq2num['A'] = 0
    seq2num['C'] = 1
    seq2num['G'] = 2
    seq2num['T'] = 3
    seq2num['N'] = 4

    wmlen = len(M)
    rowidx = arange(wmlen)

    Fs = []  #list of the likelihoods of the sequences
    for seq in open(seqs):
        if seq.startswith('>>'):
            continue
        else:
            seqlist = list(seq.strip())

            Fw = zeros(len(seqlist))  #Fw[i] contains the lilkelihood of the sequence with having either a WM ending at i or a bg at i
            for i in arange(len(seqlist)):
                if i+1 - wmlen >= 0: #start with P(s|wm) when there is enough space to put it there. Use i+1 because indexes start at 0
                    if i+1 - wmlen - 1 < 0:
                        Fw1 = 1
                    else:
                        Fw1 = Fw[i+1-wmlen-1]

                    colidx = [seq2num[seqlist[j]] for j in arange(i+1-wmlen, i+1, 1)]
                    Fww = (prod(M[rowidx, colidx]) + prod(Mrev[rowidx, colidx])) * (prior/2) * Fw1

                    #Fww = max(prod(M[rowidx, colidx]), prod(Mrev[rowidx, colidx])) * prior * Fw1
                    #Fww = (prod(Psw) + prod(Pswrev)) * prior * Fw1
                    #Fww = max(prod(Psw), prod(Pswrev)) * prior * Fw1

                else:
                    Fww = 0

                if i-1 < 0:
                    Fw1 = 1
                else:
                    Fw1 = Fw[i-1]

                Psb = basedict[seqlist[i]]
                Fwb = Psb * (1-prior) * Fw1
                #Fwb = 2*Psb * (1-prior) * Fw1
                Fwi = Fww + Fwb

                Fw[i] = Fwi

            Fs.append(Fw[-1])

    LL = sum(log(array(Fs)))

    return LL
    

def WMmain(train_seqs, test_seqs, wm, max_prior, basedict):

    M,Mrev = wm2mat(open(wm))

    #first find optimal prior (pi) just by testing different values
    priorr = arange(0.0, max_prior, 0.001)
    LLs = []
    for prior in priorr:
        L = computeLikelihood(train_seqs, M, Mrev, prior, basedict)
        print prior, L
        LLs.append(L)

    opt_prior = priorr[argmax(LLs)]
    training_LL = max(LLs)
    print opt_prior, training_LL, LLs[0], LLs[-1]

    #now compute likelihood of the test_seqs with opt_prior and with prior=0
    LL0 = computeLikelihood(test_seqs, M, Mrev, 0.0, basedict)
    LLopt = computeLikelihood(test_seqs, M, Mrev, opt_prior, basedict)
    print LLopt, LLopt-LL0

    return opt_prior, training_LL, LLopt, LLopt-LL0


def execute(cf):
    """
    This function/component refines a given weight matrix, then predicts sites for old and refined WM and true sites and on background peaks.
    With this it makes a positive predictive value - sensitivity plot.
    It also makes a logo of the refined WM.
    """

    ##Ports and parameters
    train_seqs = cf.get_input("train_sequences") #training sequences. Typically even_file
    test_seqs = cf.get_input("test_sequences") #test set. Typically odd_file
    WM1 = cf.get_input("WM1")
    WM2 = cf.get_input("WM2")
    WM3 = cf.get_input("WM3")
    WM4 = cf.get_input("WM4")
    WM5 = cf.get_input("WM5")
    WM6 = cf.get_input("WM6")
    WM7 = cf.get_input("WM7")
    WM8 = cf.get_input("WM8")
    WM9 = cf.get_input("WM9")
    WM10 = cf.get_input("WM10")
    WM11 = cf.get_input("WM11")
    WM12 = cf.get_input("WM12")
    WM13 = cf.get_input("WM13")
    WM14 = cf.get_input("WM14")
    WM15 = cf.get_input("WM15")
    WM16 = cf.get_input("WM16")
    WM17 = cf.get_input("WM17")
    WM18 = cf.get_input("WM18")
    WM19 = cf.get_input("WM19")
    WM20 = cf.get_input("WM20")
    WMdir = cf.get_input("WMdir")
    WMdir2 = cf.get_input("WMdir2")

    bestWM = cf.get_output("BestWM")
    log_file = cf.get_output("log_file")

    max_prior = cf.get_parameter("max_prior", "float")

    WMs = [i for i in[WM1, WM2, WM3, WM4, WM5, WM6, WM7, WM8, WM9, WM10, WM11, WM12, WM13, WM14, WM15, WM16, WM17, WM18, WM19, WM20] if i]

    if WMdir:
        WMs += [os.path.join(WMdir, wm) for wm in  os.listdir(WMdir)]

    if WMdir2:
        WMs += [os.path.join(WMdir2, wm) for wm in  os.listdir(WMdir2)]

    basedict = getBaseFreqs(train_seqs, test_seqs)

    #LLdiffs is the log likelihood difference between without WM model and WM model with optimal prior
    LLdiffs = []
    opt_priors = []
    for WM in WMs:
        LL = WMmain(train_seqs, test_seqs, WM, max_prior, basedict)
        LLdiffs.append(LL[3])
        opt_priors.append(LL[0])

    print LLdiffs

    #replace name in WM file with bestWM
    lines = open(WMs[argmax(LLdiffs)]).readlines()
    lines[1] = 'NA BestWM\n'
    bwm = open(bestWM, 'w')
    bwm.write(''.join(lines))
    #os.system('cp %s %s' %(WMs[bestindex], bestWM))


    l = open(log_file, 'w')

    l.write('WM_name\tWM_path\tLL_diff\topt_prior\n')

    names = ['WM_%i\t%s\t%.4f\t%s' %(i+1, WMs[i], LLdiffs[i], opt_priors[i]) for i in arange(len(WMs))]

    l.write('\n'.join(names))
    l.close()


    return 0


component_skeleton.main.main(execute)
                                                                 
