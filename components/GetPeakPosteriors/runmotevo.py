#!/usr/bin/env python

import subprocess
import os, re
from string import *
from pylab import *
import sys
import pickle

def giveMotevoParamFile(genome, inter_dir, tag, ATfreq, GCfreq, emprior, bgorder, bgprior, site_bool, loglik_bool):
    """
    Returns a parameter file for motevo.
    """

    sitefilepath = os.path.join(inter_dir, 'sites_' + tag)
    priorfilepath = os.path.join(inter_dir, 'priors_' + tag)
    loglikfile = os.path.join(inter_dir, 'loglik_' + tag)

    print '\nCreate motevo parameter file %s' %tag
    motevo_params = '\n'.join(['refspecies %s' %genome,
                               'TREE (%s: 1)' %genome,
                               'Mode TFBS',
                               'EMprior %s' %emprior,
                               'priordiff %s' %0.05,
                               'markovorderBG %s' %bgorder,
                               'bgprior %s' %bgprior,
                               'bg A %s' %ATfreq,
                               'bg T %s' %ATfreq,
                               'bg G %s' %GCfreq,
                               'bg C %s' %GCfreq,
                               'restrictparses %s' %0,
                               'priorfile %s' %priorfilepath,
                               'printsiteals %s' %0,
                               'minposterior %f' %0.01])

    if site_bool:
        motevo_params += '\nsitefile %s' %sitefilepath
    if loglik_bool:
        motevo_params += '\nloglikfile %s' %loglikfile


    params_path = os.path.join(inter_dir, 'motevo_TFBS_params_' + tag)
    pf = open(params_path, 'w')
    pf.write(motevo_params)
    pf.close()

    return (params_path, sitefilepath, priorfilepath, loglikfile)    
    


def runMotevo(motevo_path, seqs, params, WM, tag):
    """
    runs Motevo
    """
    
    pwd = os.getcwd()
    os.chdir(os.path.split(params)[0])

    print '\nrun Motevo %s' %tag
    print motevo_path + ' %s %s \"%s\"' %(seqs, params, WM)
    proc = subprocess.Popen(motevo_path + ' %s %s \"%s\"' %(seqs, params, WM),
                            stdout=subprocess.PIPE,
                            stderr= subprocess.PIPE,
                            shell=True
                            )

    stdout_value, stderr_value = proc.communicate()

    os.chdir(pwd)
    
    if proc.poll() > 0:
        print '\tstderr:', repr(stderr_value.rstrip())
        return -1
    else:
        return 0

def predict_sites(train_set, test_set, bg_train_set, WM, ATfreq, GCfreq, interm, genome, motevo_path, bg_prior, markovorder):

    ## 1. Fit prior on training set with EM
    tag = 'fitP'
    params, sites, priors, loglikfile = giveMotevoParamFile(genome, interm, tag, ATfreq, GCfreq, emprior=1, bgorder=markovorder, bgprior=bg_prior, site_bool=False, loglik_bool=False)

    tmpfile1 = os.path.join(interm, 'catted')
    os.system('cat %s %s > %s' %(train_set, bg_train_set, tmpfile1))

    r = runMotevo(motevo_path, tmpfile1, params, WM, tag)
    if r != 0:
        print >> sys.stderr, 'motevo failed ', tag
        sys.exit(1)

    os.system('rm %s' %params)
    os.system('rm %s' %tmpfile1)

    # Read priors and add them into wm files
    # prior file:
    # WM_name final_prior nr_of_sites density
    # /import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/PipeLineSource/TESTRUN/NRF1_Z2/OUTPUT/NRF1_FgBg-runmotevoPG2_1/Logo 0.016554 635.008 0.251863
    # background 0.983446 37724.8 0.748137
    # UFEwm 0 0 0

    prior_d = {}
    for line in open(priors):
        if line.startswith('WM_name'):
            continue
        t = line.strip().split()
        prior_d[t[0]] = float(t[1])

    tmpwm2 = os.path.join(interm, os.path.split(WM)[1] + '_updated')
    o = open(tmpwm2, 'w')
    for line in open(WM):
        if line.startswith('NA'):
            o.write(line)
            t = line.strip().split()
            try:
                pri = prior_d[t[1]]
            except KeyError:
                prior_d[t[1]] = 0.0
                pri = 0.0
            o.write('PW\t%s\n' %(pri))
        elif line.startswith('PW'):
            continue
        else:
            o.write(line)
    o.close()


    try:
        bgprior = prior_d['background']
    except KeyError, e:
        print >> sys.stderr, 'bg_prior Issue: ', tag
        print >> sys.stderr, e
        print >> sys.stderr, prior_d
        for i in open(priors):
            print >> sys.stderr, 'l', i.rstrip()

    os.system('rm %s' %priors)

    ## 2. Compute log-likelihood on test set with optimal prior from training set and without EM
    tag = 'TFBS'
    params, fg_sites, priors, loglikfile = giveMotevoParamFile(genome, interm, tag, ATfreq, GCfreq, emprior=0, bgorder=markovorder, bgprior=bg_prior, site_bool=True, loglik_bool=False)
    runMotevo(motevo_path, test_set, params, tmpwm2, tag)

    os.system('rm %s' %params)
    os.system('rm %s' %priors)
    os.system('rm %s' %tmpwm2)

    return fg_sites


def getDicts(sites, statsfile, regcov_dir):
    """
    sitesDict: here IDs are just IDs of the region, thus no .p1 or .rep1 etc...
    13-22 - 0.573747 BestWM hg19_chr7_915000_916000_reg1000004_9.279_+
    This function returns a dictionary: reg1000004: [(13-22, -, 0.573747), (13-22, - ,0.573747) ...]

    IDstats: store all stats for one region under region ID key
    statsfile: reg1000087.p1   4.393   1.77485636798e-05       409.430 34.822 (height, rmsd, mu, sigma)
    reg1000087: [[4.393, 1.77485636798e-05, 409.430, 34.822], ]

    IDcoords:
    regcovfile: chr10   121356000       121357250       reg1000053      1       0.5
    chr10_121356000_121357250: reg1000053
    """


    sitesDict = {}
    for line in open(sites): #e.g. 1500-1511 + 0.001162 BestWM hg19_chr13_99962000_99964000_reg1002740_14.9861052478_+
        t = line.strip().split()
        regID = t[4].split('_')[-3]
        try:
            sitesDict[regID].append( (t[0], t[1], float(t[2])) )
        except KeyError:
            sitesDict[regID] = []
            sitesDict[regID].append( (t[0], t[1], float(t[2])) )


    IDstats = {} #the peakstats file from SelectPeaks or from NormalizeHieghts contains stats for ALL fitted regions of PeakMerger outfile
    for line in open(statsfile):
        t = line.strip().split()
        try:
            IDstats[t[0].split('.')[0]].append(map(float, t[1:]) + [t[0]])
        except KeyError:
            IDstats[t[0].split('.')[0]] = []
            IDstats[t[0].split('.')[0]].append(map(float, t[1:]) + [t[0]])


    IDcoords = {}
    for fi in os.listdir(regcov_dir):
        f = open(os.path.join(regcov_dir, fi))
        t = f.readline().split()
        f.close()
        IDcoords['_'.join(t[:-3])] = t[-3]


    return sitesDict, IDstats, IDcoords



def main():
    """
    This script runs motevo to predict sites. It also creates dictionaries.
    """

    motevo_path = sys.argv[1]
    seqs = sys.argv[2]
    train_set = sys.argv[3]
    bg_train_set = sys.argv[4]
    WM =  sys.argv[5]

    pickled_sitesdict = sys.argv[6]
    pickled_idstats = sys.argv[7]
    pickled_idcoords = sys.argv[8]

    statsfile = sys.argv[9]
    regcov_dir = sys.argv[10]

    genome = sys.argv[11]
    markovorder = sys.argv[12]
    interm = sys.argv[13]

    ATfreq = 0.25
    GCfreq = 0.25
    bg_prior = 0.99

    sitesfile = predict_sites(train_set, seqs, bg_train_set, WM, ATfreq, GCfreq, interm, genome, motevo_path, bg_prior, markovorder)
    print sitesfile

    #get true sites and a post dict that contains non-refined coordinates as keys and lists of all posteriors as values
    sitesDict, IDstats, IDcoords = getDicts(sitesfile, statsfile, regcov_dir)


    pickle.dump(sitesDict, open(pickled_sitesdict, 'w'))
    pickle.dump(IDstats, open(pickled_idstats, 'w'))
    pickle.dump(IDcoords, open(pickled_idcoords, 'w'))


if __name__ == '__main__':
    main()
