#!/usr/bin/env python
import component_skeleton.main
import os, re
from string import *
import datetime, time
from pylab import *
import subprocess


def giveMotevoParamFile(genome, wmlen, inter_dir, tag, aligned, ufemodel_path, ATfreq, GCfreq, emprior, bgorder, bgprior):
    """
    Returns a parameter file for motevo.
    """

    ##UFE_models from genome_dict are not used anymore
    #UFEmodel_hg19 is UFE model for mammal species
    genome_dict = {}
    genome_dict['hg19'] = ['((((hg19:0.032973,rheMac2:0.057695):0.09821,mm9:0.352605):0.020666,(bosTau6:0.186713,(equCab2:0.107726,canFam2:0.150374):0.010431):0.032764):0.156024,monDom5:0.425899);', '/import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/MotEvo_v1.0/UFEmodels/UFEmodel_hg19']
    genome_dict['hg18'] = ['((((hg18:0.032973,rheMac2:0.057695):0.09821,mm9:0.352605):0.020666,(bosTau3:0.186713,(equCab1:0.107726,canFam2:0.150374):0.010431):0.032764):0.156024,monDom4:0.425899);', '/import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/MotEvo_v1.0/UFEmodels/UFE_mammals']
    #genome_dict['dm3'] = ['((((((dm3:0.059,droSim1:0.075):0.041,(droYak2:0.104,droEre2:0.107):0.054):0.120,droAna3:0.377):0.072,dp4:0.397):0.061,droWil1:0.536):0.020,((droVir3:0.196,droMoj3:0.255):0.073,droGri2:0.291):0.337);', '/import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/MotEvo_v1.0/UFEmodels/UFEmodel_dm3']
    genome_dict['dm3'] = ['((((((dm3:0.059,droSim1:0.075):0.041,(droYak2:0.104,droEre2:0.107):0.054):0.120,droAna3:0.377):0.072,dp4:0.397):0.061,droWil1:0.536):0.020,((droVir3:0.196,droMoj3:0.255):0.073,droGri2:0.291):0.337);', '/import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/MotEvo_v1.0/UFEmodels/dm3UFEparallel/UFEmodel_dm3']
    genome_dict['mm9'] = ['((((hg19:0.032973,rheMac2:0.057695):0.09821,mm9:0.352605):0.020666,(bosTau7:0.186713,(equCab2:0.107726,canFam2:0.150374):0.010431):0.032764):0.156024,monDom5:0.425899);', '/import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/MotEvo_v1.0/UFEmodels/UFEmodel_mm9']


    sitefilepath = os.path.join(inter_dir, 'sites_' + tag)
    priorfilepath = os.path.join(inter_dir, 'priors_' + tag)
    loglikfile = os.path.join(inter_dir, 'loglik_' + tag)


    print '\nCreate motevo parameter file %s' %tag
    print 'aligned', aligned
    if aligned:
        motevo_params = '\n'.join(['refspecies %s' %genome,
                                   'TREE %s' %genome_dict[genome][0],
                                   'Mode TFBS',
                                   'EMprior %s' %emprior,
                                   'priordiff %s' %0.05,
                                   'UFEwmprior %s' %200,
                                   'UFEwmfile %s' %ufemodel_path,
                                   'UFEwmlen %s' %wmlen,
                                   'UFEprint %s' %0,
                                   'markovorderBG %s' %bgorder,
                                   'bgprior %s' %bgprior,
                                   'bg A %s' %ATfreq,
                                   'bg T %s' %ATfreq,
                                   'bg G %s' %GCfreq,
                                   'bg C %s' %GCfreq,
                                   'restrictparses %s' %0,
                                   'sitefile %s' %sitefilepath,
                                   'priorfile %s' %priorfilepath,
                                   'printsiteals %s' %0,
                                   'minposterior %f' %0.0,
                                   'loglikfile %s' %loglikfile])
    else:
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
                                   'sitefile %s' %sitefilepath,
                                   'priorfile %s' %priorfilepath,
                                   'printsiteals %s' %0,
                                   'minposterior %f' %0.0,
                                   'loglikfile %s' %loglikfile])       

    params_path = os.path.join(inter_dir, 'motevo_TFBS_params_' + tag)
    pf = open(params_path, 'w')
    pf.write(motevo_params)
    return (params_path, sitefilepath, priorfilepath, loglikfile)    

    

    
def runMotevo(motevo_path, seqs, params, WM, interm, tag):
    """
    runs Motevo
    """
    
    # pwd = os.getcwd()
    # os.chdir(interm)

    print '\nrun Motevo %s' %tag
    proc = subprocess.Popen(motevo_path + ' %s %s %s' %(seqs, params, WM),
                            stdout=subprocess.PIPE,
                            stderr= subprocess.PIPE,
                            shell=True
                            )

    stdout_value, stderr_value = proc.communicate()
    print stdout_value
    print stderr_value

    # os.chdir(pwd)
    
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
    infile = cf.get_input("infile") #file containing non-similar looking WM paths
    regions = cf.get_input("regions") #training set. Typically even_file
    basefreqs = cf.get_input("BaseFrequencies")
    ufemodel_path = cf.get_input("UFEmodel")
    keepWM = cf.get_input("keepWM")

    WMdir = cf.get_output("WMdir")
    bestWM = cf.get_output("BestWM")
    log_file = cf.get_output("log_file")
    interm = cf.get_output("intermediate")
    outplot = cf.get_output("contribution_plot")

    genome = cf.get_parameter('genome', 'string')
    motevo_path = cf.get_parameter('motevo_path', 'string')
    aligned = cf.get_parameter("aligned", "boolean")
    ll_co = cf.get_parameter('loglik_cutoff', 'float')

    os.mkdir(interm)
    os.mkdir(WMdir)

    # Read stuff in
    f = open(basefreqs)
    ATfreq = float(f.readline().strip().split()[1])
    GCfreq = float(f.readline().strip().split()[1])
    f.close()


    # Passed WMs:
    # WM_0    /import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/PipeLineSource/TESTRUN/NRF1_Z2/OUTPUT/NRF1_FgBg-filterwms/WMdir/WM_6
    # WM_1    /import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/PipeLineSource/TESTRUN/NRF1_Z2/OUTPUT/NRF1_FgBg-trimwm/outdir/trimmedWM_3
    # WM_2    /import/bc2/home/nimwegen/GROUP/WMs/Mammals/CurrentWMs/NRF1.p2
    wmnames = []
    wmpaths = []

    for line in open(infile):
        if line.startswith('WM_'):
            if line.startswith('WM_name'):
                continue
            else:
                t = line.strip().split()
                wmnames.append(t[0])
                wmpaths.append(t[1])


    # Compute stuff: optimal priors and then likelihood of test set
    logliks = []
    badwms = []
    goodwmpaths = []

    figure()
    plotx = []
    ploty = []

    for i, WM in enumerate(wmpaths):

        tmpwm = os.path.join(interm, 'tmpwm')

        os.system('cat %s > %s' %(' '.join(wmpaths[:i+1]), tmpwm))

        wmlen = 12 #just some number. I can't take length of one WM because I have several WMs

        # 1. Fit prior on training set with EM
        tag = 'compLL_%i' %(i+1)
        params, sites, priors, loglikfile = giveMotevoParamFile(genome, wmlen, interm, tag, aligned, ufemodel_path, ATfreq, GCfreq, emprior=1, bgorder=0, bgprior=0.99)
        r = runMotevo(motevo_path, regions, params, tmpwm, interm, tag)
        if r != 0:
            print 'motevo failed ', tag
            sys.exit(1)


        a = loadtxt(loglikfile, usecols=[1])
        ll = sum(a)

        # always take first WM. If an added WM doesn't help more than ll_co, then it gets thrown out
        if i > 0:
            if ll - logliks[-1] > ll_co: 
                plotx.append(len(logliks)+1)
                ploty.append(ll)

                logliks.append(ll)
                goodwmpaths.append(wmpaths[i])
            else:
                badwms.append(wmpaths[i])
                
        else:
            plotx.append(len(logliks)+1)
            ploty.append(ll)

            goodwmpaths.append(wmpaths[i])
            logliks.append(ll)


    # if there is a keepWM defined, also compute likelihood for this one but add it anyway
    if keepWM:
        if not keepWM in goodwmpaths:
            print 'Given WM %s was added to list' %keepWM

            tmpwm = os.path.join(interm, 'tmpwm')

            os.system('cat %s > %s' %(' '.join(goodwmpaths + [keepWM]), tmpwm))

            wmlen = 12 #just some number. I can't take length of one WM because I have several WMs

            # 1. Fit prior on training set with EM
            tag = 'compLL_keepWM' 
            params, sites, priors, loglikfile = giveMotevoParamFile(genome, wmlen, interm, tag, aligned, ufemodel_path, ATfreq, GCfreq, emprior=1, bgorder=0, bgprior=0.99)
            r = runMotevo(motevo_path, regions, params, tmpwm, interm, tag)
            if r != 0:
                print 'motevo failed ', tag
                sys.exit(1)

            a = loadtxt(loglikfile, usecols=[1])
            ll = sum(a)

            plotx.append(len(logliks)+1)
            ploty.append(ll)

            goodwmpaths.append(keepWM)
            logliks.append(ll)


    print plotx, ploty
    plot([0]+plotx, [0]+ploty, 'o-')
    xticks([0]+plotx, [0]+['WM_%i' %(i+1) for i in arange(len(goodwmpaths))])
    #xlabel('WM number')
    ylabel('log-likelihood given WM number')
    #ylim([0, max(ploty)+ max(ploty)*0.05])
    savefig(outplot)
    close()

    l = open(log_file, 'w')

    l.write('Name\tPath\tlog_likelihood\n')

    names = ['WM_%i\t%s\t%.2f' %(i+1, goodwmpaths[i], logliks[i]) for i in arange(len(goodwmpaths))]
    badnames = ['%s' %(i) for i in badwms]

    l.write('\n'.join(names))
    l.write('\n\nNon-contributing WMs:\n')
    l.write('\n'.join(badnames))
    l.close()

    for i in arange(len(goodwmpaths)):
        os.system('cp %s %s' %(goodwmpaths[i], os.path.join(WMdir, 'WM_' + str(i+1))))



    return 0


component_skeleton.main.main(execute)
                                                                 
