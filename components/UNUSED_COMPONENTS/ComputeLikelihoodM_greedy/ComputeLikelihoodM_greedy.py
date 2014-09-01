#!/usr/bin/env python
import component_skeleton.main
import os, re
from string import *
import datetime, time
from pylab import *
import subprocess


def giveMotevoParamFile(genome, wmlen, inter_dir, tag, aligned, ufemodel_path, ATfreq, GCfreq, emprior, bgorder, bgprior, UFEwmprior):
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
    
    print '\nrun Motevo %s' %tag
    proc = subprocess.Popen(motevo_path + ' %s %s \"%s\"' %(seqs, params, WM),
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
    This function runs a greedy algorithm to compute the log-likelihood of sequences given the best combination of input WMs.
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
    keepWM = cf.get_input("keepWM")
    WMdir = cf.get_input("WMdir")
    WMdir2 = cf.get_input("WMdir2")
    basefreqs = cf.get_input("BaseFrequencies")
    ufemodel_path = cf.get_input("UFEmodel")

    loglik_all_motifs = cf.get_output("loglik_all_motifs")
    loglik_contributing_motifs = cf.get_output("loglik_contributing_motifs")
    loglik_combined = cf.get_output("loglik_combined")
    interm = cf.get_output("intermediate")
    outplot = cf.get_output("loglik_plot")
    WMoutdir = cf.get_output("WMoutdir")

    genome = cf.get_parameter('genome', 'string')
    motevo_path = cf.get_parameter('motevo_path', 'string')
    aligned = cf.get_parameter("aligned", "boolean")
    slim = cf.get_parameter("slim", "boolean")
    loglik_co = cf.get_parameter("loglik_cutoff", "float") #this log-likelihood cut-off of course depends on the number of input sequences... 20 should be used for 500 seqs...

    os.mkdir(interm)
    os.mkdir(WMoutdir)

    # adapt loglik_co to number of input sequences: the one given is thought for 500 sequences... scale it linearly
    N = 0
    for l in open(test_set):
        if l.startswith('>>'):
            N+=1
        else:
            continue

    loglik_co /= 500.
    loglik_co *= N

    loglik_co = max(5.0, loglik_co) #too low cut-offs result in too many matrices...

    print 'log_likelihood cut-off: %s (%i sequences)' %(loglik_co, N)

    ## Read stuff in
    WMs = [i for i in[WM1, WM2, WM3, WM4, WM5, WM6, WM7, WM8, WM9, WM10, WM11, WM12, WM13, WM14, WM15, WM16, WM17, WM18, WM19, WM20] if i]

    if WMdir:
        WMs += [os.path.join(WMdir, wm) for wm in  os.listdir(WMdir)]

    if WMdir2:
        WMs += [os.path.join(WMdir2, wm) for wm in  os.listdir(WMdir2)]

    if keepWM: #add keepWM. If keepWM gets thrown out because it's not contributing, it still gets added later on so that one can see in any case how this WM performs
        WMs.append(keepWM)

    f = open(basefreqs)
    ATfreq = float(f.readline().strip().split()[1])
    GCfreq = float(f.readline().strip().split()[1])
    f.close()


    ## start greedy algorithm:

    final_wms_indxs = []
    final_lls = []

    init_lls = []
    init_priors = []
    init = True

    tmpwm = os.path.join(interm, 'tmpwm')

    badindxs = [] # list of indices in WMs that do not contribute to the log-lik at any arbitrary step. 

    # cut-off for log-likelihood. Best would be zero, but then it takes long. Here a WM has to improve log-likelihood by at least 10
    # I think a motif can't increase log-likelihood more when combined with other motifs than just by itself. (That's why taking bad WMs out doesn't hurt)

    j = 0

    while 1:
        # Compute stuff: optimal priors and then likelihood of test set

        logliks = []
        WMindxs = []

        j += 1

        for i, WM in enumerate(WMs):

            if i in final_wms_indxs + badindxs:
                continue

            WMindxs.append(i)

            os.system('cat \"%s\" > %s' %(WM, tmpwm))
            for wmidx in final_wms_indxs:
                os.system('cat \"%s\" >> %s' %(WMs[wmidx], tmpwm))

            wmlen = len(open(WM).readlines())-4

            # 1. Fit prior on training set with EM
            tag = 'fitP_%i_%i' %(j, i+1)
            params, sites, priors, loglikfile = giveMotevoParamFile(genome, wmlen, interm, tag, aligned, ufemodel_path, ATfreq, GCfreq, emprior=1, bgorder=0, bgprior=0.99, UFEwmprior=200)
            r = runMotevo(motevo_path, train_set, params, tmpwm, interm, tag)
            if r != 0:
                print 'motevo failed ', tag
                sys.exit(1)

            os.system('rm %s' %sites)

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

            tmpwm2 = os.path.join(interm, 'tmpwm2')
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
                else:
                    o.write(line)

            o.close()

            try:
                UFEwmprior = prior_d['UFEwm']
            except KeyError:
                UFEwmprior = 0

            print 'UFEwmprior: ', UFEwmprior
            bgprior = prior_d['background']

            # 2. Compute log-likelihood on test set with optimal prior from training set and without EM
            tag = 'compLL_%i_%i' %(j, i+1)
            params, sites, priors, loglikfile = giveMotevoParamFile(genome, wmlen, interm, tag, aligned, ufemodel_path, ATfreq, GCfreq, emprior=0, bgorder=0, bgprior=bgprior, UFEwmprior=UFEwmprior)
            runMotevo(motevo_path, test_set, params, tmpwm2, interm, tag)

            os.system('rm %s' %sites)

            a = loadtxt(loglikfile, usecols=[1])
            ll = sum(a)

            logliks.append(ll)

            if init:
                init_lls.append(ll)
                init_priors.append([prior_d[k] for k in prior_d if k != 'UFEwm' and k!= 'background'][0])


        if init:
            init = False

            if slim:
                break

        indxs_logliks = zip(WMindxs, logliks)
        sorted_indxs_logliks = sorted(indxs_logliks, key= lambda k: k[1], reverse=True)

        print sorted_indxs_logliks
        print final_lls

        if len(final_lls) >= 1:
            if sorted_indxs_logliks[0][1] - final_lls[-1] <= loglik_co:
                print sorted_indxs_logliks[0][1],  final_lls[-1]
                print 'no loglik improvement'
                break
        else:
            if sorted_indxs_logliks[0][1] <= loglik_co:
                print 'Warning: no motif from given library contributes! One motif added by default.\n'
                final_wms_indxs.append(sorted_indxs_logliks[0][0])
                final_lls.append(sorted_indxs_logliks[0][1])
                break

        # sort out non-contributing motifs (they have to contribute at least a little at every step)
        for wmll in sorted_indxs_logliks:
            if len(final_lls) == 0: #if single WMs give log-lik smaller than 0
                if wmll[1] <= loglik_co:
                    badindxs.append(wmll[0])
            else:
                if wmll[1] - final_lls[-1] <= loglik_co: #if WM groups do not improve log-lik from before
                    badindxs.append(wmll[0])

        final_wms_indxs.append(sorted_indxs_logliks[0][0])
        final_lls.append(sorted_indxs_logliks[0][1])

        # check this first, because already after the first or second run all motifs can be in badindxs and nothing in sorted_indxs_logliks. So the likelihood convergence check would fail...
        if len(final_wms_indxs) + len(badindxs) == len(WMs): # at the end all WMs are either contributing or not. In case the log-likelihood doesn't improve anymore even earlier, the break statement above will kick in
            print 'len(final+bad) == len(WMs)'
            break

        if j == len(WMs): #This 
            print 'j == len(WMs)'
            break

        print badindxs


    print WMs
    print final_wms_indxs
    print final_lls
    print init_lls
    print init_priors

    goodwmpaths = [WMs[i] for i in final_wms_indxs]

    keepWMcontributes = True
    # if there is a keepWM defined, also compute likelihood for this one but add it anyway (same function as above)
    if keepWM:
        if not keepWM in goodwmpaths:
            print 'Given WM %s was added to the list even though it is not contributing' %keepWM

            # now compute likelihood given contributing WMs plus keepWM
            tmpwm = os.path.join(interm, 'tmpwm')

            os.system('cat %s > %s' %(' '.join(goodwmpaths + [keepWM]), tmpwm))

            wmlen = 12 #just some number. I can't take length of one WM because I have several WMs

            # 1. Fit prior on training set with EM
            tag = 'fitP_keepWM'
            params, sites, priors, loglikfile = giveMotevoParamFile(genome, wmlen, interm, tag, aligned, ufemodel_path, ATfreq, GCfreq, emprior=1, bgorder=0, bgprior=0.99, UFEwmprior=200)
            r = runMotevo(motevo_path, train_set, params, tmpwm, interm, tag)
            if r != 0:
                print 'motevo failed ', tag
                sys.exit(1)

            os.system('rm %s' %sites)

            # Read priors and add them into wm files
            prior_d = {}
            for line in open(priors):
                if line.startswith('WM_name'):
                    continue
                t = line.strip().split()
                prior_d[t[0]] = float(t[1])

            tmpwm2 = os.path.join(interm, 'tmpwm2')
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
                else:
                    o.write(line)

            o.close()

            try:
                UFEwmprior = prior_d['UFEwm']
            except KeyError:
                UFEwmprior = 0

            print 'UFEwmprior: ', UFEwmprior
            bgprior = prior_d['background']

            # 2. Compute log-likelihood on test set with optimal prior from training set and without EM
            tag = 'compLL_keepWM'
            params, sites, priors, loglikfile = giveMotevoParamFile(genome, wmlen, interm, tag, aligned, ufemodel_path, ATfreq, GCfreq, emprior=0, bgorder=0, bgprior=bgprior, UFEwmprior=UFEwmprior)
            runMotevo(motevo_path, test_set, params, tmpwm2, interm, tag)

            os.system('rm %s' %sites)

            a = loadtxt(loglikfile, usecols=[1])
            ll = sum(a)

            goodwmpaths.append(keepWM)
            final_lls.append(ll)
            final_wms_indxs.append(len(WMs)-1) #keepWM is always in last position of WMs list

            keepWMcontributes = False
        else:
            keepWMcontributes = True


    plot(arange(len(goodwmpaths)+1), [0] + final_lls, 'r-')
    plot(arange(len(goodwmpaths)+1), [0] + final_lls, 'ko')


    # Get the names for the plot. If WM file name was not generated by me (i.e. something like WM_1 or WM_2), then take this name. This is the case when I process known WMs.
    if not keepWMcontributes:
        names = ['']
        for i in arange(len(goodwmpaths)):
            wmfile_name = os.path.split(goodwmpaths[i])[1]
            if not wmfile_name.startswith('WM_'):
                names.append(wmfile_name)
            else:
                names.append('WM_%i' %(i+1))
        names[-1] += '\n(given)'

    else:
        names = ['']
        for i in arange(len(goodwmpaths)):
            wmfile_name = os.path.split(goodwmpaths[i])[1]
            if not wmfile_name.startswith('WM_'):
                if goodwmpaths[i] == keepWM:
                    names.append('%s\n(given)' %(wmfile_name))
                else:
                    names.append(wmfile_name)
            else:
                if goodwmpaths[i] == keepWM:
                    names.append('WM_%i\n(given)' %(i+1))
                else:
                    names.append('WM_%i' %(i+1))


    xticks(arange(len(goodwmpaths)+1), names, rotation=30)
    ylabel('log-likelihood')
    savefig(outplot)
    savefig(outplot.rstrip('.pdf'))


    l = open(loglik_all_motifs, 'w')

    l.write('WM_path\tlog_likelihood\topt_prior\n')

    names = ['%s\t%.4f\t%s' %(WMs[i], init_lls[i], init_priors[i]) for i in arange(len(WMs))]

    l.write('\n'.join(names))
    l.close()


    # here write WM paths of the WMs that were copied to the output directory, so that the naming is consistent afterwards.
    # dictionary to associate old WM paths with new ones. Copy WMs to out_dir
    old_new_name_dict = {}
    for i in arange(len(goodwmpaths)):
        os.system('cp \"%s\" %s/WM_%i' %(goodwmpaths[i], WMoutdir, i+1))
        old_new_name_dict[goodwmpaths[i]] = '%s/WM_%i' %(WMoutdir, i+1)


    lc = open(loglik_combined, 'w')

    lc.write('WM_name\tWM_path\tlog_likelihood\n')
    names = ['WM_%i\t%s\t%s' %(i+1, old_new_name_dict[goodwmpaths[i]], final_lls[i]) for i in arange(len(goodwmpaths)) ]
    lc.write('\n'.join(names))
    lc.close()



    l = open(loglik_contributing_motifs, 'w')

    l.write('WM_name\tWM_path\tlog_likelihood\topt_prior\n')

    #sort WMs here not by there order of contribution but by the log-likelihood they achieve by them selves (i.e. init_lls)
    sorted_final_wms_indxs = sorted(final_wms_indxs, key = lambda idx: init_lls[idx], reverse=True)
    names = ['WM_%i\t%s\t%.4f\t%s' %(j+1, old_new_name_dict[WMs[i]], init_lls[i], init_priors[i]) for j, i in enumerate(sorted_final_wms_indxs)]

    l.write('\n'.join(names))
    l.close()



    return 0


component_skeleton.main.main(execute)
                                                                 
