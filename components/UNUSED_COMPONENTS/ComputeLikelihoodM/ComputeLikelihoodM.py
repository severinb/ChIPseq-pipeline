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
    #print 'Command: ./motevo %s %s %s' %(seqs, params, WM)
    proc = subprocess.Popen(motevo_path + ' %s %s \"%s\"' %(seqs, params, WM),
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
    train_set = cf.get_input("train_set") #training set. Typically even_file
    test_set = cf.get_input("test_set") #test set. Typically odd_file
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
    basefreqs = cf.get_input("BaseFrequencies")
    ufemodel_path = cf.get_input("UFEmodel")

    bestWM = cf.get_output("BestWM")
    log_file = cf.get_output("log_file")
    interm = cf.get_output("intermediate")

    genome = cf.get_parameter('genome', 'string')
    motevo_path = cf.get_parameter('motevo_path', 'string')
    aligned = cf.get_parameter("aligned", "boolean")

    os.mkdir(interm)



    # Read stuff in
    WMs = [i for i in[WM1, WM2, WM3, WM4, WM5, WM6, WM7, WM8, WM9, WM10, WM11, WM12, WM13, WM14, WM15, WM16, WM17, WM18, WM19, WM20] if i]

    if WMdir:
        WMs += [os.path.join(WMdir, wm) for wm in  os.listdir(WMdir)]

    if WMdir2:
        WMs += [os.path.join(WMdir2, wm) for wm in  os.listdir(WMdir2)]

    f = open(basefreqs)
    ATfreq = float(f.readline().strip().split()[1])
    GCfreq = float(f.readline().strip().split()[1])
    f.close()


    # Compute stuff: optimal priors and then likelihood of test set
    optpriors = []
    logliks = []

    for i, WM in enumerate(WMs):

        wmlen = len(open(WM).readlines())-4

        # 1. Fit prior on training set with EM
        tag = 'fitP_%i' %(i+1)
        params, sites, priors, loglikfile = giveMotevoParamFile(genome, wmlen, interm, tag, aligned, ufemodel_path, ATfreq, GCfreq, emprior=1, bgorder=0, bgprior=0.99)
        r = runMotevo(motevo_path, train_set, params, WM, interm, tag)
        if r != 0:
            print 'motevo failed ', tag
            sys.exit(1)

        # prior file:
        # WM_name final_prior nr_of_sites density
        # /import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/PipeLineSource/TESTRUN/NRF1_Z2/OUTPUT/NRF1_FgBg-runmotevoPG2_1/Logo 0.016554 635.008 0.251863
        # background 0.983446 37724.8 0.748137
        # UFEwm 0 0 0

        optprior = float(open(priors).readlines()[1].split()[1])
        bgprior=(1-optprior)
        print bgprior

        # 2. Compute log-likelihood on test set with optimal prior from training set and without EM
        tag = 'compLL_%i' %(i+1)
        params, sites, priors, loglikfile = giveMotevoParamFile(genome, wmlen, interm, tag, aligned, ufemodel_path, ATfreq, GCfreq, emprior=0, bgorder=0, bgprior=bgprior)
        runMotevo(motevo_path, train_set, params, WM, interm, tag)

        a = loadtxt(loglikfile, usecols=[1])
        ll = sum(a)

        logliks.append(ll)
        optpriors.append(optprior)

    print logliks



    #replace name in WM file with bestWM
    lines = open(WMs[argmax(logliks)]).readlines()
    lines[1] = 'NA BestWM\n'
    bwm = open(bestWM, 'w')
    bwm.write(''.join(lines))


    l = open(log_file, 'w')

    l.write('WM_name\tWM_path\tlog_likelihood\topt_prior\n')

    names = ['WM_%i\t%s\t%.4f\t%s' %(i+1, WMs[i], logliks[i], optpriors[i]) for i in arange(len(WMs))]

    l.write('\n'.join(names))
    l.close()


    return 0


component_skeleton.main.main(execute)
                                                                 
