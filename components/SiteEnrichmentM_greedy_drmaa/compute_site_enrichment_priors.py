#!/usr/bin/env python

import os, re, sys
from string import *
from pylab import *
import subprocess

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


    #sitefilepath = os.path.join(inter_dir, 'sites_' + tag)
    sitefilepath = '/scratch/motevo_sites_' + tag + '_%i_%s' %(randint(10000000), datetime.datetime.now().strftime("%Y%m%d%H%M%S%f"))
    #priorfilepath = os.path.join(inter_dir, 'priors_' + tag)
    priorfilepath = '/scratch/priors_' + tag + '_%i_%s' %(randint(10000000), datetime.datetime.now().strftime("%Y%m%d%H%M%S%f"))
    #loglikfile = os.path.join(inter_dir, 'loglik_' + tag)
    loglikfile = '/scratch/loglik_' + tag + '_%i_%s' %(randint(10000000), datetime.datetime.now().strftime("%Y%m%d%H%M%S%f"))


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


    #params_path = os.path.join(inter_dir, 'motevo_TFBS_params_' + tag)
    params_path = '/scratch/motevo_TFBS_params_' + tag + '_%i_%s' %(randint(10000000), datetime.datetime.now().strftime("%Y%m%d%H%M%S%f"))
    try:
        pf = open(params_path, 'w')
    except IOError,e:
        print >> sys.stderr, 'params_path Issue:'
        print >> sys.stderr, e
        os.system('uname -n 1>&2')
        os.system('ls -laht /scratch 1>&2')
        os.system('df -h /scratch 1>&2')
        os.system('rm %s' %params_path)
        sys.exit(1)

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
    print >> sys.stderr, stderr_value
    
    if proc.poll() > 0:
        print >> sys.stderr, repr(stderr_value.rstrip())
        return -1
    else:
        return 0


def compute_site_enrichment(train_set, test_set, bg_train_set, background_set, tmpwm, ufemodel_path, ATfreq, GCfreq, interm, genome, motevo_path, aligned, tag_i, input_seqs_num, background_seqs_num, pseudocount, bg_prior, ufewm_prior):

    ## 1. Fit prior on training set with EM
    tag = 'fitP_' + tag_i
    params, sites, priors, loglikfile = giveMotevoParamFile(genome, interm, tag, aligned, ufemodel_path, ATfreq, GCfreq, emprior=1, bgorder=0, bgprior=bg_prior, UFEwmprior=ufewm_prior, site_bool=False, loglik_bool=False)

    #tmpfile1 = os.path.join(interm, 'catted_' + tag_i)
    tmpfile1 = '/scratch/catted_' + tag_i + '_%i_%s' %(randint(10000000), datetime.datetime.now().strftime("%Y%m%d%H%M%S%f"))
    os.system('cat %s %s > %s' %(train_set, bg_train_set, tmpfile1))

    r = runMotevo(motevo_path, tmpfile1, params, tmpwm, interm, tag)
    if r != 0:
        print >> sys.stderr, 'motevo failed ', tag
        sys.exit(1)

    #os.system('rm %s' %sites)
    os.system('rm %s' %params)
    os.system('rm %s' %tmpfile1)
    #os.system('rm %s' %loglikfile)

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

    #tmpwm2 = os.path.join(interm, tmpwm + '_updated')
    tmpwm2 = '/scratch/' + os.path.split(tmpwm)[1] + '_updated_%i_%s' %(randint(10000000), datetime.datetime.now().strftime("%Y%m%d%H%M%S%f"))
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

    os.system('rm %s' %priors)

    ## 2. Compute log-likelihood on test set with optimal prior from training set and without EM
    tag = 'compLL_' + tag_i
    params, sites, priors, loglikfile = giveMotevoParamFile(genome, interm, tag, aligned, ufemodel_path, ATfreq, GCfreq, emprior=0, bgorder=0, bgprior=bgprior, UFEwmprior=UFEwmprior, site_bool=True, loglik_bool=True)
    runMotevo(motevo_path, test_set, params, tmpwm2, interm, tag)

    posts = loadtxt(sites, usecols=[2])
    seqs = loadtxt(sites, usecols=[4], dtype=str)
    # test whether posts and seqs are really arrays (lists) and not just empty or single floats..
    try:
        len(posts)
    except TypeError:
        posts = array([posts])
        seqs = array([seqs])

    unique_seqs = unique(seqs)
    seqs_n = [] # expected number of binding sites for test set peaks
    for seq in unique_seqs:
        seqs_n.append(sum(posts[where(seqs==seq)]))

    seqs_n += list(zeros(input_seqs_num - len(seqs_n))) #has to be of the same length for every motif!
    seqs_n = array(seqs_n)
    sumseqs = sum(seqs_n)
    #pseudocount = 0.01
    seqs_n += pseudocount #pseudo count to not get -inf in the logs for zero sites sequences.

    try:
        a = loadtxt(loglikfile, usecols=[1])
    except ValueError, e:
        print >> sys.stderr, 'loglik Issue: ', tag
        print >> sys.stderr, e
        for i in open(loglikfile):
            print >> sys.stderr, 'l', i.rstrip()

    ll = sum(a)

    os.system('rm %s' %loglikfile)
    os.system('rm %s' %sites)
    os.system('rm %s' %params)
    os.system('rm %s' %priors)

    ## 3. Compute number of predicted binding sites on background set with optimal prior from training set and without EM
    tag = 'compBG_' + tag_i
    params, sites, priors, loglikfile = giveMotevoParamFile(genome, interm, tag, aligned, ufemodel_path, ATfreq, GCfreq, emprior=0, bgorder=0, bgprior=bgprior, UFEwmprior=UFEwmprior, site_bool=True, loglik_bool=True)
    runMotevo(motevo_path, background_set, params, tmpwm2, interm, tag)

    try:
        posts = loadtxt(sites, usecols=[2])
        N = sum(posts)
    except ValueError, e:
        print >> sys.stderr, 'sites Issue: ', tag
        print >> sys.stderr, e
        for i in open(sites):
            print >> sys.stderr, 'l', i.rstrip()

    try:
        a = loadtxt(loglikfile, usecols=[1])
    except ValueError, e:
        print >> sys.stderr, 'loglik Issue: ', tag
        print >> sys.stderr, e
        for i in open(loglikfile):
            print >> sys.stderr, 'l', i.rstrip()

    llbg = sum(a)
    llnorm = ll-(0.1*llbg)

    os.system('rm %s' %loglikfile)
    os.system('rm %s' %sites)
    os.system('rm %s' %params)
    os.system('rm %s' %tmpwm2)
    os.system('rm %s' %priors)

    site_enrichment = sum(log(seqs_n)) - (len(seqs_n) * log(N + (background_seqs_num*pseudocount) + sum(seqs_n)))
    #site_enrichment = llnorm

    #also print how site numbers are distributed over peaks:
    bin_edges = list(linspace(pseudocount, 2, 10))
    if max(seqs_n) > 2.:
        bin_edges += [max(seqs_n)]
    figure()
    h = hist(seqs_n, bins=bin_edges)
    close()
    print 'Number of peaks: ', len(seqs_n)
    print 'Number of sites in peaks: ', sumseqs
    print 'Number of sites in peaks (with pseudocounts): ', sum(seqs_n)
    print 'Distribution of sites in peaks ([site_counts], [site_bins]): ', h[0], h[1]
    print 'Number of sites in background, bg_pseudocount: ', N, background_seqs_num*pseudocount
    print 'Site enrichment: ', site_enrichment
    print ''
    print 'll: ', ll
    print 'llbg: ', llbg
    print 'llnorm: ', llnorm

    return site_enrichment, prior_d


if __name__ == '__main__':

    if len(sys.argv) != 20:
        print '\nUsage: python prog.py tagsroot wmroot outfileroot train_set test_set background_set ufemodel_path ATfreq GCfreq interm genome motevo_path aligned input_seqs_number background_seqs_num pseudocount\n'
        sys.exit(0)

    tagsroot = sys.argv[1] #a file containing tags of WM indices. filesfileroot.SGE_TASK_ID is actual file
    wmroot = sys.argv[2] #wmroot.tag (from tagsroot) will be the actual WM to use
    outfileroot = sys.argv[3] #files that contains site enrichment for every line in filesfileroot. outfileroot.SGE_TASK_ID is used.
    train_set = sys.argv[4]
    test_set = sys.argv[5]
    bg_train_set = sys.argv[6]
    background_set = sys.argv[7]
    ufemodel_path = sys.argv[8]
    ATfreq = float(sys.argv[9])
    GCfreq = float(sys.argv[10])
    interm = sys.argv[11]
    genome = sys.argv[12]
    motevo_path = sys.argv[13]
    aligned = int(sys.argv[14])
    input_seqs_num = int(sys.argv[15])
    background_seqs_num = int(sys.argv[16])
    pseudocount = float(sys.argv[17])
    bg_prior = sys.argv[18]
    ufewm_prior = sys.argv[19]

    taskid = os.environ['SGE_TASK_ID']
    print taskid
    filesfile = tagsroot + '.%s' %taskid
    outfile = outfileroot + '.%s' %taskid

    print filesfile, outfile

    o = open(outfile, 'w')
    for line in open(filesfile):

        tag = line.strip()
        tmpwm = wmroot + '.' + tag

        site_enrichment, prior_d = compute_site_enrichment(train_set, test_set, bg_train_set, background_set, tmpwm, ufemodel_path, ATfreq, GCfreq, interm, genome, motevo_path, aligned, tag, input_seqs_num, background_seqs_num, pseudocount, bg_prior, ufewm_prior)

        print  [prior_d[k] for k in prior_d if k != 'UFEwm' and k!= 'background']
        o.write('%s\t%s\t' %(tag, site_enrichment))
        for tf in prior_d:
            o.write('%s %s;' %(tf, prior_d[tf]))
        o.write('\n')


        os.system('rm %s' %tmpwm)

    o.close()

