#!/usr/bin/env python

import component_skeleton.main
import subprocess
import os, re
from string import *
import datetime, time
from pylab import *
import sys
import drmaa

def findFraglen(readfiles):

    fraglens = []
    rlist = readfiles.split()
    for rf in rlist:
        ##get Fragment length from the log of FragmentLength component
        fraglenRoot = os.path.join(os.path.split(rf)[0], '..', 'intermediate')
        interdir = os.listdir(fraglenRoot)
        for f in interdir:
            m = re.match('\S*.res$', os.path.join(fraglenRoot, f))
            if m:
                fraglenF = m.group()
                break
        fraglen = int(float(open(fraglenF).read().strip().split()[-1]))
        fraglens.append(fraglen)
        print '\nFound: %s with fragment length: %s\n' %(rf,fraglen)

    meanFraglen = mean(fraglens)
    maxFraglen= max(fraglens)
    print '\nFragment Length: %s' %fraglens

    return meanFraglen


def getBaseFreqs(f):
    """           
    This function computes base frequencies for input sequences
    """

    bd = {}
    bd['A'] = 0.
    bd['C'] = 0.
    bd['G'] = 0.
    bd['T'] = 0.

    for line in open(f):
        if line.startswith('>'):
            continue
        for b in bd:
            bases = list(line.strip())
            bd[b] += bases.count(b)

    print bd
    #normalize to frequencies
    tot = sum(bd.values())
    for i in bd:
        bd[i] /= tot

    print bd
    #take same frequencies for A and T or C and G respectively, because we are using double stranded DNA at the end.

    prec = 3
    ATfreq = round((bd['A']+bd['T'])/2., prec)
    GCfreq = round((bd['C']+bd['G'])/2., prec)

    #round frequencies to some precision and be sure that they sum up to exactly 1
    #round one frequency (AT or GC) down and the other up

    while (ATfreq + GCfreq) != 0.5:
        if (ATfreq + GCfreq) > 0.5:
            if rand(1) >= 0.5:
                ATfreq -= 1.0/(10**prec)
            else:
                GCfreq -= 1.0/(10**prec)
        else:
            if rand(1) >= 0.5:
                ATfreq += 1.0/(10**prec)
            else:
                GCfreq += 1.0/(10**prec)

    return ATfreq, GCfreq



def runMotevo(motevo_path, seqs, train_set, bg_train_set, WM, interm, project_leader, pickled_sitesdict, pickled_idstats, pickled_idcoords, statsfile, regcov_dir, instance_name, queue_name, genome, markovorder):
    """
    runs Motevo
    """

    #here call motevo with drmaa
    JOB_PARAM = '-q %s -P %s -e %s/motevo_job.stderr -o %s/motevo_job.stdout -j n -w n -N %s-predSites -cwd -V -b y' %(queue_name, project_leader, os.path.split(interm)[0], os.path.split(interm)[0], instance_name)

    s = drmaa.Session()
    s.initialize()

    jt = s.createJobTemplate()
    jt.nativeSpecification = JOB_PARAM

    jt.remoteCommand = './runmotevo_updated.py'

    jt.args = [motevo_path, seqs, train_set, bg_train_set, WM, pickled_sitesdict, pickled_idstats, pickled_idcoords, statsfile, regcov_dir, genome, str(markovorder), interm]

    job = s.runJob(jt)

    print 'submitted', job

    retval = s.wait(job, drmaa.Session.TIMEOUT_WAIT_FOREVER)

    print "motevo job finished with status: ", retval.hasExited

    s.deleteJobTemplate(jt)
    s.exit()


    os.system('cat %s/motevo_job.stderr' %os.path.split(interm)[0])
    os.system('cat %s/motevo_job.stdout' %os.path.split(interm)[0])

    return retval.hasExited


def splitRegCovs(regcov_dir, interm, fpj):
    """
    This function creates input files for drmaa. It creates files with fpj regcov file paths of.
    """

    count = 0
    i = 0

    froot = os.path.join(interm, 'regcovfile')

    for f in os.listdir(regcov_dir):
        if i >= fpj:
            i = 0
            o.close()
        if i == 0:
            count += 1
            outfile = froot + '.%i' %count
            o = open(outfile, 'w')

        name = os.path.join(regcov_dir, f)
        o.write(name + '\n')
        i += 1

    return froot, count


def submitJobs(regcov_root, pickled_sd, pickled_ids, pickled_idc, plotdir, fraglen, minpost, peakstats_root, tfbsstats_root, count, project_leader, instance_name, queue_name):

    JOB_PARAM = '-q %s -P %s -e %s/job.stderr -o %s/job.stdout -j n -w n -N %s-combPost -cwd -V -b y' %(queue_name, project_leader, os.path.split(plotdir)[0], os.path.split(plotdir)[0], instance_name)

    s = drmaa.Session()
    s.initialize()

    jt = s.createJobTemplate()
    jt.nativeSpecification = JOB_PARAM

    jt.remoteCommand = './combinePosteriors_faster.py'

    jt.args = [pickled_ids, pickled_idc, pickled_sd, regcov_root, tfbsstats_root, peakstats_root, plotdir, str(fraglen), str(minpost)]

    jobs = s.runBulkJobs(jt,1,count,1)

    print count
    print 'submitted', jobs[0]

    s.synchronize(jobs, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)

    s.deleteJobTemplate(jt)
    s.exit()


def combineFiles(peakstats_root, tfbsstats_root, count, peakstats, TFBSstats):

    op = open(peakstats, 'w')
    op.write('#chrom\tstart\tend\tpeakID\tzscore\tquality\tsummed_posterior\n')
    ot = open(TFBSstats, 'w')
    ot.write('#chrom\tstart\tend\tpeakID\tdistance\tposterior\tTFBS_coverage\n')

    for i in arange(1, count+1, 1):
        for pl in open(peakstats_root + '.%i' %i):
            op.write(pl)

        for tl in open(tfbsstats_root + '.%i' %i):
            ot.write(tl)

    op.close()
    ot.close()


def execute(cf):
    """
    This component gives true regions (determined by a posterior cut-off on TFBS).
    It produces some plots: 
        -histogram of region posteriors (one with summed posteriors and one with maximum TFBS posterior per region)
        -plots peak coverage (from RegionCoverage) plots with TFBSs (above 0.5 posterior cut-off)
    """

    ##Ports and parameters
    regions = cf.get_input("regions") #sequences of candidate regions
    train_set = cf.get_input("train_set") #train_set and bg_train_set are used to estimate a prior for the WM that is then used to predict sites on the regions.
    bg_train_set = cf.get_input("bg_train_set")
    regcov_dir = cf.get_input("RegCov_dir")
    WM = cf.get_input("WM") 
    statsfile = cf.get_input("statsfile")

    peakstats = cf.get_output("peakstats")
    TFBSstats = cf.get_output("TFBSstats")
    interm = cf.get_output("intermediate")
    log_file = cf.get_output("log_file")


    plotdir = os.path.join(os.path.split(interm)[0], 'peak_plots') #cf.get_output("peak_plots")

    genome = cf.get_parameter("genome", "string")
    minpost = cf.get_parameter("minposterior", "float")
    motevo_path = cf.get_parameter("motevo_path", "string")
    markovorder = cf.get_parameter("markovorder", "int")
    read_files = cf.get_parameter("read_files", "string")
    fpj = cf.get_parameter("files_per_job", "float")
    project_leader = cf.get_parameter("project_leader", "string")
    queue_name_motevo = cf.get_parameter("queue_name_motevo", "string")
    queue_name = cf.get_parameter("queue_name", "string")

    instance_name = cf.get_metadata("instanceName")

    T1 = datetime.datetime.now()

    ##Main function
    os.mkdir(interm)
    os.mkdir(plotdir)

    fraglen = findFraglen(read_files)

    #get parameter file and predicted sites for best WM
    pickled_sd = os.path.join(interm, 'sitesDict')
    pickled_ids = os.path.join(interm, 'IDstats')
    pickled_idc = os.path.join(interm, 'IDcoords')

    retval = runMotevo(motevo_path, regions, train_set, bg_train_set, WM, interm, project_leader, pickled_sd, pickled_ids, pickled_idc, statsfile, regcov_dir, instance_name, queue_name_motevo, genome, markovorder)
    if not retval:
        return 1

    T2 = datetime.datetime.now()

    #create Plots and build peakstats and TFBSstats files
    regcov_root, count = splitRegCovs(regcov_dir, interm, fpj)

    peakstats_root = os.path.join(interm, 'peakstatsfile')
    tfbsstats_root = os.path.join(interm, 'tfbsstatsfile')

    submitJobs(regcov_root, pickled_sd, pickled_ids, pickled_idc, plotdir, fraglen, minpost, peakstats_root, tfbsstats_root, count, project_leader, instance_name, queue_name)

    combineFiles(peakstats_root, tfbsstats_root, count, peakstats, TFBSstats)


    # make archive:
    pwd = os.getcwd()
    os.system('cd %s && tar -czf ../%s.tar.gz . && cd %s' %(plotdir, os.path.split(plotdir)[1], pwd))
    os.system('rm -r %s' %plotdir)


    # clean up: remove pickled dictionaries and file chunks
    os.system('rm %s %s %s' %(pickled_sd, pickled_ids, pickled_idc))
    for fname in os.listdir(interm):
        f = os.path.join(interm, fname)
        if fname.startswith('peakstatsfile') or fname.startswith('tfbsstatsfile') or fname.startswith('regcovfile'):
            os.system('rm %s' %f)


    T3 = datetime.datetime.now()


    timetext = '\n'.join(['Running time:',
                          '\t-Predicting sites on given regions and loading data into dictionaries: %s' %(T2-T1),
                          '\t-Plotting coverage profiles with TFBSs and computing peak posteriors: %s' %(T3-T2),
                          ])


    lf = open(log_file, 'w')
    lf.write(timetext)
    lf.close()

    print 'Running time: %s' %str(T3-T1)

    return 0


component_skeleton.main.main(execute)
                                                                 
