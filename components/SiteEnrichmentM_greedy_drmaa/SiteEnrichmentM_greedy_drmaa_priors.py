#!/usr/bin/env python

import component_skeleton.main
import os, re
from string import *
from pylab import *
import drmaa
#import compute_site_enrichment

def execute(cf):
    """
    This function runs a greedy algorithm to compute site enrichment of sequences given the best combination of input WMs.
    Site enrichment is the probability of pulling down the peak sequences in a pool of shuffled sequences.
    P(s) = ns/N, P(s) is the probability of pulling down sequence s. ns is the number of binding sites on s and N is the number of binding sites in the background pool.
    site enrichment is then (ns/N)**S, where S is the number of sequences.
    Here it is in log-space, i.e. = sum(log(ns/N))
    """

    ##Ports and parameters
    train_set = cf.get_input("train_set") #training set. Typically even_file
    test_set = cf.get_input("test_set") #test set. Typically odd_file
    bg_train_set = cf.get_input("bg_train_set") #training set. Typically shuffled even_file
    background_set = cf.get_input("background_set") #typically shuffled test set (or random genomic regions)
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
    WMdir3 = cf.get_input("WMdir3")
    basefreqs = cf.get_input("BaseFrequencies")
    ufemodel_path = cf.get_input("UFEmodel")
    instance_name = cf.get_metadata("instanceName")

    loglik_all_motifs = cf.get_output("loglik_all_motifs")
    loglik_top_motifs = cf.get_output("loglik_top_motifs")
    loglik_contributing_motifs = cf.get_output("loglik_contributing_motifs")
    loglik_combined = cf.get_output("loglik_combined")
    outplot = cf.get_output("loglik_plot")

    interm = os.path.join(os.path.split(outplot)[0], "intermediate")

    genome = cf.get_parameter('genome', 'string')
    motevo_path = cf.get_parameter('motevo_path', 'string')
    aligned = cf.get_parameter("aligned", "boolean")
    slim = cf.get_parameter("slim", "boolean")
    loglik_co = cf.get_parameter("site_enrichment_cutoff", "float") #this log-likelihood cut-off of course depends on the number of input sequences... 20 should be used for 500 seqs...
    pseudocount = cf.get_parameter("pseudocount", "float")
    project_leader = cf.get_parameter('project_leader', 'string')
    queue_name = cf.get_parameter('queue_name', 'string')
    top_motifs = cf.get_parameter('top_motifs', 'int')

    os.mkdir(interm)

    if aligned:
        aligned = 1
    else:
        aligned = 0

    # adapt loglik_co to number of input sequences: the one given is thought for 500 sequences... scale it linearly
    input_seqs_num = 0
    for l in open(test_set):
        if l.startswith('>>'):
            input_seqs_num+=1
        else:
            continue

    loglik_co = (loglik_co/500.)*input_seqs_num

    background_seqs_num = 0
    for l in open(background_set):
        if l.startswith('>>'):
            background_seqs_num+=1
        else:
            continue

    # randomly, we would expect the following site enrichment: input_seqs_num * log(1/background_seqs_num)
    # One could thus subtract this number from the site enrichment of each motif...
    no_enrichment = input_seqs_num * log(1./background_seqs_num)

    print 'site enrichment cut-off: %s (%i sequences)' %(loglik_co, input_seqs_num)
    print 'no enrichment: %s' %(no_enrichment)

    ## Read stuff in
    WMs = [i for i in[WM1, WM2, WM3, WM4, WM5, WM6, WM7, WM8, WM9, WM10, WM11, WM12, WM13, WM14, WM15, WM16, WM17, WM18, WM19, WM20] if i]

    if WMdir:
        WMs += [os.path.join(WMdir, wm) for wm in  os.listdir(WMdir)]

    if WMdir2:
        WMs += [os.path.join(WMdir2, wm) for wm in  os.listdir(WMdir2)]

    if WMdir3:
        WMs += [os.path.join(WMdir3, wm) for wm in  os.listdir(WMdir3)]

    if keepWM: #add keepWM. If keepWM gets thrown out because it's not contributing, it still gets added later on so that one can see in any case how this WM performs
        WMs.append(keepWM)

    if basefreqs:
        f = open(basefreqs)
        ATfreq = f.readline().strip().split()[1]
        GCfreq = f.readline().strip().split()[1]
        f.close()
    else:
        ATfreq = '0.25'
        GCfreq = '0.25'

    ## start greedy algorithm:
    final_wms_indxs = []
    final_lls = [] #should be final_site_enrichment or so, generally, everything that is called loglik or so here should be site_enrichment

    init = True

    tmpwm = os.path.join(interm, 'tmpwm')
    outfileroot = os.path.join(interm, 'outfile')
    tagroot = os.path.join(interm, 'tags')

    badindxs = [] # list of indices in WMs that do not contribute to the log-lik at any arbitrary step. 

    final_set_prior_d = {} #just used in the first round
    final_set_prior_d['background'] = '0.99'
    # final_set_prior_d['UFEwm'] = '200'
    final_set_prior_d['UFEwm'] = '0'

    j = 0

    while 1:

        j += 1

        #create files with tags
        tasks_per_job = 5
        jobidx = 0
        closed_flag = True
        tasks_count = 0

        for i, WM in enumerate(WMs):
            if i in final_wms_indxs + badindxs:
                continue

            if closed_flag:
                jobidx += 1
                tf = open(tagroot + '.%s' %jobidx, 'w')
                closed_flag = False

            tag = '%i_%i' %(j, i+1)
            tf.write('%s\n' %(tag))

            tmpwm_tag = tmpwm + '.' + tag
            # create WM for the tag
            if not init:
                o = open(tmpwm_tag, 'w')
                for line in open(WM):
                    if line.startswith('NA'):
                        o.write(line)
                        t = line.strip().split()
                        pri = 0.0001 #start with a not too big prior. If I took 0 it, motevo'd give me too high prior for the new WM
                        o.write('PW\t%s\n' %(pri))
                    else:
                        o.write(line.strip())
                        o.write('\n')
                for wmidx in final_wms_indxs:
                    for line in open(WMs[wmidx]):
                        if line.startswith('NA'):
                            o.write(line)
                            t = line.strip().split()
                            pri = final_set_prior_d[t[1]]
                            o.write('PW\t%s\n' %(pri))
                        else:
                            o.write(line.strip())
                            o.write('\n')
                o.close()
            else:
                os.system('cat \"%s\" > %s' %(WM, tmpwm_tag))
                for wmidx in final_wms_indxs:
                    os.system('cat \"%s\" >> %s' %(WMs[wmidx], tmpwm_tag))

            tasks_count += 1
            if tasks_count == tasks_per_job:
                tf.close()
                closed_flag = True
                tasks_count = 0

        tf.close()


        JOB_PARAM = '-q %s -P %s -e %s/job.stderr -o %s/job.stdout -j n -w n -N %s-CSE.%i -cwd -V -b y' %(queue_name, project_leader, os.path.split(interm)[0], os.path.split(interm)[0], instance_name, j)

        s = drmaa.Session()
        s.initialize()

        jt = s.createJobTemplate()
        jt.nativeSpecification = JOB_PARAM
        jt.remoteCommand = './compute_site_enrichment_priors.py'
        jt.args = [tagroot, tmpwm, outfileroot, train_set, test_set, bg_train_set, background_set, ufemodel_path, ATfreq, GCfreq, interm, genome, motevo_path, str(aligned), str(input_seqs_num), str(background_seqs_num), str(pseudocount), final_set_prior_d['background'], final_set_prior_d['UFEwm']]

        jobs = s.runBulkJobs(jt,1,jobidx,1)

        print 'submitted', jobs[0]

        s.synchronize(jobs, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)

        s.deleteJobTemplate(jt)
        s.exit()


        WMindx_site_enrichment_prior = []
        for jidx in arange(1, jobidx + 1, 1):

            outfile = outfileroot + '.' + str(jidx)

            for l in open(outfile):
                t = l.strip().split('\t')
                tag = t[0]
                site_enrichment = float(t[1])
                prior_d = {}
                for pl in t[2].rstrip(';').split(';'):
                    plt = pl.split()
                    prior_d[plt[0]] = plt[1]

                i_tag = int(tag.split('_')[1]) - 1

                site_enrichment -= no_enrichment #normalize by what one would expect by randomly fishing for sequences
                WMindx_site_enrichment_prior.append((i_tag, site_enrichment, prior_d))


        WMindxs = []
        site_enrichments = []
        priors = [] # this list is only meaningful for the first round, when just one matrix is used
        prior_ds = []
        for elm in sorted(WMindx_site_enrichment_prior, key = lambda k: k[0]):
            WMindxs.append(elm[0])
            site_enrichments.append(elm[1])
            prior_d = elm[2]
            priors.append([prior_d[k] for k in prior_d if k != 'UFEwm' and k!= 'background'][0])
            prior_ds.append(prior_d)

        if init:
            init_lls = list(site_enrichments)
            init_priors = list(priors)
            init = False

            if slim:
                break


        indxs_logliks = zip(WMindxs, site_enrichments, prior_ds)
        sorted_indxs_logliks = sorted(indxs_logliks, key= lambda k: k[1], reverse=True)


        print 'sorted indxs: ', sorted_indxs_logliks
        print 'final_lls: ', final_lls


        if len(final_lls) >= 1:
            if (sorted_indxs_logliks[0][1] - final_lls[-1]) <= loglik_co:
                print sorted_indxs_logliks[0][1],  final_lls[-1]
                print 'no site enrichment improvement'
                break

            # sort out non-contributing motifs (they have to contribute at least a little at every step)
            # for wmll in sorted_indxs_logliks:
            #     if (wmll[1] - final_lls[-1]) <= loglik_co: #if WM groups do not improve log-lik from before
            #         badindxs.append(wmll[0])

            final_wms_indxs.append(sorted_indxs_logliks[0][0])
            final_lls.append(sorted_indxs_logliks[0][1])
            final_set_prior_d = sorted_indxs_logliks[0][2]

            print 'bad WMs: ', badindxs
            print 'total WMs %i - bad WMs %i - contributing WMs %i = WMs %i' %(len(WMs), len(badindxs), len(final_wms_indxs), len(WMs)-len(badindxs)-len(final_wms_indxs))

            # check this first, because already after the first or second run all motifs can be in badindxs and nothing in sorted_indxs_logliks. So the likelihood convergence check would fail...
            if len(final_wms_indxs) + len(badindxs) == len(WMs): # at the end all WMs are either contributing or not. In case the log-likelihood doesn't improve anymore even earlier, the break statement above will kick in
                print 'len(final+bad) == len(WMs)'
                break

            if j == len(WMs): #This 
                print 'j == len(WMs)'
                break


        else: #for the first round always add the first matrix independent of whether it is bad or good...
            final_wms_indxs.append(sorted_indxs_logliks[0][0])
            final_lls.append(sorted_indxs_logliks[0][1])
            final_set_prior_d = sorted_indxs_logliks[0][2]

            # sort out non-contributing motifs (they have to improve at least a little compared to no enrichment...)
            # do not do this, just to give all matrices a chance to help the best one.
            # for wmll in sorted_indxs_logliks:
            #     if wmll[1] <= loglik_co:
            #         badindxs.append(wmll[0])
            continue


    print WMs
    print final_wms_indxs
    print final_lls
    print init_lls
    print init_priors

    goodwmpaths = [WMs[i] for i in final_wms_indxs]
    goodwmnames = [os.path.split(i)[1] for i in goodwmpaths]

    keepWMcontributes = True
    # if there is a keepWM defined, also compute likelihood for this one but add it anyway (same function as above)
    if keepWM:
        if not os.path.split(keepWM)[1] in goodwmnames:
            print 'Given WM %s was added to the list even though it is not contributing' %keepWM

            tmpwm = os.path.join(interm, 'tmpwm')
            outfileroot = os.path.join(interm, 'outfile_keepWM')
            tagroot = os.path.join(interm, 'tags_keepWM')

            tag = 'keepWM'
            #create files with tags
            tf = open(tagroot + '.1', 'w')
            tf.write(tag + '\n')
            tf.close()

            tmpwm_tag = tmpwm + '.' + tag
            # create WM for the tag
            os.system('cat %s > %s' %(' '.join(goodwmpaths + [keepWM]), tmpwm_tag))


            JOB_PARAM = '-q %s -P %s -e %s/job.stderr -o %s/job.stdout -j n -w n -N %s-CSE_keepWM -cwd -V -b y' %(queue_name, project_leader, os.path.split(interm)[0], os.path.split(interm)[0], instance_name)
    
            s = drmaa.Session()
            s.initialize()

            jt = s.createJobTemplate()
            jt.nativeSpecification = JOB_PARAM
            jt.remoteCommand = './compute_site_enrichment_priors.py'
            jt.args = [tagroot, tmpwm, outfileroot, train_set, test_set, bg_train_set, background_set, ufemodel_path, ATfreq, GCfreq, interm, genome, motevo_path, str(aligned), str(input_seqs_num), str(background_seqs_num), str(pseudocount), '0.99', '200']

            jobs = s.runBulkJobs(jt,1,1,1)
            print 'submitted', jobs[0]

            s.synchronize(jobs, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)
            s.deleteJobTemplate(jt)
            s.exit()

            outfile = outfileroot + '.1'
            for l in open(outfile):
                t = l.strip().split()
                tag = t[0]
                site_enrichment = float(t[1])
                prior = t[2]
                print prior
                site_enrichment -= no_enrichment

            goodwmpaths.append(keepWM)
            final_lls.append(site_enrichment)
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
            names.append(wmfile_name)
        names[-1] += '\n(given)'

    else:
        names = ['']
        for i in arange(len(goodwmpaths)):
            wmfile_name = os.path.split(goodwmpaths[i])[1]
            if goodwmpaths[i] == keepWM:
                names.append('%s\n(given)' %(wmfile_name))
            else:
                names.append(wmfile_name)


    xticks(arange(len(goodwmpaths)+1), names, rotation=90)
    ylabel('Enrichment Score')
    tight_layout()
    savefig(outplot)
    savefig(outplot.rstrip('.pdf'))


    print len(WMs)
    print len(init_lls)
    print len(init_priors)

    l = open(loglik_all_motifs, 'w')

    l.write('WM_path\tenrichment_score\topt_prior\n')

    names = ['%s\t%.4f\t%s' %(WMs[i], init_lls[i], init_priors[i]) for i in arange(len(WMs))]

    l.write('\n'.join(names))
    l.close()


    l = open(loglik_top_motifs, 'w')
    WM_ll_p = zip(WMs, init_lls, init_priors)
    WM_ll_p_sorted = sorted(WM_ll_p, key = lambda k: k[1], reverse = 1)

    l.write('WM_path\tenrichment_score\topt_prior\n')

    names = ['%s\t%.4f\t%s' %(WM_ll_p_sorted[i][0], WM_ll_p_sorted[i][1], WM_ll_p_sorted[i][2]) for i in arange(len(WM_ll_p)) if i < top_motifs]

    l.write('\n'.join(names))
    l.close()



    lc = open(loglik_combined, 'w')

    lc.write('WM_path\tenrichment_score\n')
    names = ['%s\t%.1f' %(goodwmpaths[i], final_lls[i]) for i in arange(len(goodwmpaths)) ]
    lc.write('\n'.join(names))
    lc.close()


    l = open(loglik_contributing_motifs, 'w')

    l.write('WM_path\tenrichment_score\topt_prior\n')

    #sort WMs here not by there order of contribution but by the log-likelihood they achieve by themselves (i.e. init_lls)
    sorted_final_wms_indxs = sorted(final_wms_indxs, key = lambda idx: init_lls[idx], reverse=True)
    names = ['%s\t%.1f\t%s' %(WMs[i], init_lls[i], init_priors[i]) for i in sorted_final_wms_indxs]

    l.write('\n'.join(names))
    l.close()

    # clean up: remove tags file chunks
    os.system('rm -r %s' %interm)

    return 0


component_skeleton.main.main(execute)
                                                                 
