#!/usr/bin/env python
import component_skeleton.main
import subprocess
import os, re
from string import *
import ROC2
import datetime, time
import signal
from pylab import *
from scipy import interpolate

def giveMotevoParamFile(genome, inter_dir, tag, aligned, ufemodel_path, ATfreq, GCfreq, emprior, bgorder, bgprior, UFEwmprior, site_bool, loglik_bool):
    """
    Returns a parameter file for motevo.
    """

    ##UFE_models from genome_dict are not used anymore
    #UFEmodel_hg19 is UFE model for mammal species
    genome_dict = {}
    genome_dict['hg19'] = ['((((hg19:0.032973,rheMac2:0.057695):0.09821,mm9:0.352605):0.020666,(bosTau6:0.186713,(equCab2:0.107726,canFam2:0.150374):0.010431):0.032764):0.156024,monDom5:0.425899);', '/import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/MotEvo_v1.0/UFEmodels/UFEmodel_hg19']
    genome_dict['hg18'] = ['((((hg18:0.032973,rheMac2:0.057695):0.09821,mm9:0.352605):0.020666,(bosTau3:0.186713,(equCab1:0.107726,canFam2:0.150374):0.010431):0.032764):0.156024,monDom4:0.425899);', '/import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/MotEvo_v1.0/UFEmodels/UFE_mammals']
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
                                   'UFEwmprior %s' %UFEwmprior,
                                   'UFEwmfile %s' %ufemodel_path,
                                   'UFEwmlen 10', #keep this constant, for several WMs it gets difficult anyway.
                                   'UFEprint %s' %0,
                                   'markovorderBG %s' %bgorder,
                                   'bgprior %s' %bgprior,
                                   'bg A %s' %ATfreq,
                                   'bg T %s' %ATfreq,
                                   'bg G %s' %GCfreq,
                                   'bg C %s' %GCfreq,
                                   'restrictparses %s' %0,
                                   'priorfile %s' %priorfilepath,
                                   'printsiteals %s' %0,
                                   'minposterior %f' %0.001])

        if site_bool:
            motevo_params += '\nsitefile %s' %sitefilepath
        if loglik_bool:
            motevo_params += '\nloglikfile %s' %loglikfile

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
                                   'priorfile %s' %priorfilepath,
                                   'printsiteals %s' %0,
                                   'minposterior %f' %0.001])

        if site_bool:
            motevo_params += '\nsitefile %s' %sitefilepath
        if loglik_bool:
            motevo_params += '\nloglikfile %s' %loglikfile


    params_path = os.path.join(inter_dir, 'motevo_TFBS_params_' + tag)
    pf = open(params_path, 'w')
    pf.write(motevo_params)
    pf.close()

    return (params_path, sitefilepath, priorfilepath, loglikfile)    
    

    
def runMotevo(motevo_path, seqs, params, WM, interm, tag):
    """
    runs Motevo
    """
    
    print '\nrun Motevo %s' %tag
    proc = subprocess.Popen(motevo_path + ' %s %s \"%s\"' %(seqs, params, WM),
                            stdout=subprocess.PIPE,
                            stderr= subprocess.PIPE,
                            shell=True
                            )

    stdout_value, stderr_value = proc.communicate()
    print stdout_value
    print >> sys.stderr, stderr_value
    
    if proc.poll() > 0:
        print >> sys.stderr, repr(stderr_value.rstrip())
        return -1
    else:
        return 0


def predict_sites(train_set, test_set, bg_train_set, bg_test_set, tmpwm, ufemodel_path, ATfreq, GCfreq, interm, genome, motevo_path, aligned, tag_i, bg_prior, ufewm_prior):

    ## 1. Fit prior on training set with EM
    tag = 'fitP_' + tag_i
    params, sites, priors, loglikfile = giveMotevoParamFile(genome, interm, tag, aligned, ufemodel_path, ATfreq, GCfreq, emprior=1, bgorder=0, bgprior=bg_prior, UFEwmprior=ufewm_prior, site_bool=False, loglik_bool=False)

    tmpfile1 = os.path.join(interm, 'catted_' + tag_i)
    os.system('cat %s %s > %s' %(train_set, bg_train_set, tmpfile1))

    r = runMotevo(motevo_path, tmpfile1, params, tmpwm, interm, tag)
    if r != 0:
        print >> sys.stderr, 'motevo failed ', tag
        sys.exit(1)

    # os.system('rm %s' %params)
    # os.system('rm %s' %tmpfile1)

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

    #os.system('rm %s' %priors)

    tmpwm2 = os.path.join(interm, tmpwm + '_updated')
    o = open(tmpwm2, 'w')
    for line in open(tmpwm):
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
        UFEwmprior = prior_d['UFEwm']
    except KeyError:
        UFEwmprior = 0

    #print 'UFEwmprior: ', UFEwmprior
    try:
        bgprior = prior_d['background']
    except KeyError, e:
        print >> sys.stderr, 'bg_prior Issue: ', tag
        print >> sys.stderr, e
        print >> sys.stderr, prior_d
        for i in open(priors):
            print >> sys.stderr, 'l', i.rstrip()

    # os.system('rm %s' %priors)

    ## 2. Compute log-likelihood on test set with optimal prior from training set and without EM
    tag = 'fg_sites_' + tag_i
    params, fg_sites, priors, loglikfile = giveMotevoParamFile(genome, interm, tag, aligned, ufemodel_path, ATfreq, GCfreq, emprior=0, bgorder=0, bgprior=bgprior, UFEwmprior=UFEwmprior, site_bool=True, loglik_bool=False)
    runMotevo(motevo_path, test_set, params, tmpwm2, interm, tag)

    # os.system('rm %s' %params)
    # os.system('rm %s' %priors)

    ## 3. Compute number of predicted binding sites on background set with optimal prior from training set and without EM
    tag = 'bg_sites_' + tag_i
    params, bg_sites, priors, loglikfile = giveMotevoParamFile(genome, interm, tag, aligned, ufemodel_path, ATfreq, GCfreq, emprior=0, bgorder=0, bgprior=bgprior, UFEwmprior=UFEwmprior, site_bool=True, loglik_bool=False)
    runMotevo(motevo_path, bg_test_set, params, tmpwm2, interm, tag)


    # os.system('rm %s' %params)
    # os.system('rm %s' %tmpwm2)
    # os.system('rm %s' %priors)


    return fg_sites, bg_sites

    

def ROCanalysis(FGsitefiles, BGsitefiles, ROCplot, interm, peaks, bgPeaks, wm_names):
    """
    This function produces a ROC plot with functions from the ROC.py module.
    """

    if len(FGsitefiles) != len(BGsitefiles):
        print '\nError: Number of foreground and background site files is not the same!\n'
        return

    AUCs = []

    AUCdict = {} #dict with index: [senslist, ppvlist, AUC]

    for i in arange(len(FGsitefiles)):

        FGposteriors = ROC2.get_peak_posteriors(FGsitefiles[i], peaks)
        BGposteriors = ROC2.get_peak_posteriors(BGsitefiles[i], bgPeaks)

        senslist, ppvlist = ROC2.TP_P_FN_N(FGposteriors, BGposteriors)

        #AUC by trapezoid integration
        AUC = dot([senslist[j] - senslist[j+1] for j in arange(len(senslist)-1)], [mean(ppvlist[j:j+1]) for j in arange(len(ppvlist)-1)])
        AUCs.append(AUC)

        AUCdict[i] = [senslist, ppvlist, AUC]

        figure()
        plot(senslist, ppvlist)
        xlabel('sensitivity')
        ylabel('positive predictive value')
        ylim([0,1.1])
        title('AUC: %.3f' %AUC)
        savefig(os.path.join(interm, 'WM_%i.pdf' %(i+1)))
        savefig(os.path.join(interm, 'WM_%i' %(i+1))) #png for website
        close()


    figure()
    for i in sorted(AUCdict, key= lambda k: AUCdict[k][-1], reverse=1)[:4]:
        plot(AUCdict[i][0], AUCdict[i][1], label='%s' %(wm_names[i]))

    xlabel('sensitivity')
    ylabel('positive predictive value')
    ylim([0,1.1])
    legend(bbox_to_anchor=(1.1, 1.1))
    savefig(ROCplot)
    close()


    return AUCs
    

def execute(cf):
    """
    This component takes a bunch of WMs and then, for every single WM:
        - fits a prior on the training pool (training set + shuffled training set)
        - predicts TFBS using the fitted prior on the testing pool (testing set + shuffled testing set)
        - plots a precision-recall curve using the site predictions
    """

    ##Ports and parameters
    train_set = cf.get_input("train_set") #alignments or sequences
    bg_train_set = cf.get_input("bg_train_set") #shuffled alignments or sequences
    test_set = cf.get_input("test_set")
    bg_test_set = cf.get_input("bg_test_set")
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
    WMlist = cf.get_input("WMlist")
    ufemodel_path = cf.get_input("UFEmodel")
    basefreqs = cf.get_input("BaseFrequencies")

    interm = cf.get_output("intermediate")
    ROCplot = cf.get_output("sens_ppv")
    log_file = cf.get_output("log_file")

    genome = cf.get_parameter("genome", "string")
    motevo_path = cf.get_parameter("motevo_path", "string")
    aligned = cf.get_parameter("aligned", "boolean")


    if not aligned:
        ufewm_prior = None
        ufemodel_path = None
        ATfreq = 0.25
        GCfreq = 0.25
    else:
        #read base frequencies
        ufewm_prior = 200
        f = open(basefreqs)
        ATfreq = float(f.readline().strip().split()[1])
        GCfreq = float(f.readline().strip().split()[1])
        f.close()

    os.mkdir(interm)
    
    WMs = [i for i in[WM1, WM2, WM3, WM4, WM5, WM6, WM7, WM8, WM9, WM10, WM11, WM12, WM13, WM14, WM15, WM16, WM17, WM18, WM19, WM20] if i]
    if WMdir:
        WMs += [os.path.join(WMdir, wm) for wm in  os.listdir(WMdir)]
    if WMdir2:
        WMs += [os.path.join(WMdir2, wm) for wm in  os.listdir(WMdir2)]
    if WMlist: #here are all the WMs.
        f = open(WMlist)
        wmlist = f.readlines()[1:]
        f.close()
        for w in wmlist:
            WMs.append(w.split()[0])

    FGsitefiles = []
    BGsitefiles = []
    wm_names = []

    for i, WM in enumerate(WMs):

        wmlen = len(open(WM).readlines())-4
        wm_names.append(os.path.split(WM)[1])
        bg_prior = 0.99

        #get sitefile for foreground and background
        tag = 'wm_%i' %i
        (fg_sites, bg_sites) = predict_sites(train_set, test_set, bg_train_set, bg_test_set, WM, ufemodel_path, ATfreq, GCfreq, interm, genome, motevo_path, aligned, tag, bg_prior, ufewm_prior)

        # (fgParams, fgSites) = giveMotevoParamFile('TFBS', genome, wmlen, interm, tag, aligned, ufemodel_path, ATfreq, GCfreq)
        # runMotevo(motevo_path, peaks, fgParams, WM, interm, tag)
        FGsitefiles.append(fg_sites)

        # (bgParams, bgSites) = giveMotevoParamFile('TFBS', genome, wmlen, interm, tag, aligned, ufemodel_path, ATfreq, GCfreq)
        # runMotevo(motevo_path, bgPeaks, bgParams, WM, interm, tag)
        BGsitefiles.append(bg_sites)


    #ROC analysis
    print 'create ROC plot'
    AUCs = ROCanalysis(FGsitefiles, BGsitefiles, ROCplot, interm, test_set, bg_test_set, wm_names)


    l = open(log_file, 'w')

    names = ['WM_%i\t%s\t%.4f' %(i+1, WMs[i], AUCs[i]) for i in arange(len(WMs))]

    l.write('WM_name\tWM_path\tAUC\n')
    l.write('\n'.join(names))
    l.close()

    #remove sites files because they are using a lot of space
    for i in FGsitefiles + BGsitefiles:
        os.system('rm %s' %i)


    return 0


component_skeleton.main.main(execute)
                                                                 
