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


def giveMotevoParamFile(genome, wmlen, inter_dir):
    """
    Returns a parameter file for motevo.
    """

    sitefilepath = os.path.join(inter_dir, 'sites')
    priorfilepath = os.path.join(inter_dir, 'priors')


    print '\nCreate motevo parameter file'
    motevo_params = '\n'.join(['refspecies %s' %genome,
                               'TREE (%s: 1)' %genome,
                               'Mode TFBS',
                               'EMprior %s' %1,
                               'priordiff %s' %0.01,
                               'markovorderBG %s' %1,
                               'bgprior %s' %0.99,
                               'restrictparses %s' %0,
                               'sitefile %s' %sitefilepath,
                               'priorfile %s' %priorfilepath,
                               'printsiteals %s' %0,
                               'minposterior %f' %0.01])            

    params_path = os.path.join(inter_dir, 'motevo_TFBS_params')
    pf = open(params_path, 'w')
    pf.write(motevo_params)

    return (params_path, sitefilepath)    

    
def runMotevo(motevo_path, alignments, params, WM, interm, project_leader):
    """
    runs Motevo
    """

    pwd = os.getcwd()
    os.chdir(interm)

    #here call motevo with drmaa
    JOB_PARAM = '-q fs_long -P %s -e %s/motevo_job.stderr -o %s/motevo_job.stdout -j n -N combPost -cwd -V -b y' %(project_leader, interm, interm)

    s = drmaa.Session()
    s.initialize()

    jt = s.createJobTemplate()
    jt.nativeSpecification = JOB_PARAM

    jt.remoteCommand = motevo_path

    jt.args = [alignments, params, WM]

    job = s.runJob(jt)

    print 'submitted', job

    s.wait(job, drmaa.Session.TIMEOUT_WAIT_FOREVER)

    s.deleteJobTemplate(jt)
    s.exit()

    os.chdir(pwd)


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
            fh.close()
        if i == 0:
            count += 1
            outfile = froot + '.%i' %count
            o = open(outfile, 'w')

        name = os.path.join(regcov_dir, f)
        o.write(name + '\n')
        i += 1

    return froot, count


def submitJobs(regcov_root, pickled_sd, pickled_ids, pickled_idc, plotdir, fraglen, minpost, peakstats_root, tfbsstats_root, count, project_leader):

    JOB_PARAM = '-q fs_long -P %s -e %s/job.stderr -o %s/job.stdout -j n -N combPost -cwd -V -b y' %(project_leader, os.path.split(regcov_root)[0], os.path.split(regcov_root)[0])

    s = drmaa.Session()
    s.initialize()

    jt = s.createJobTemplate()
    jt.nativeSpecification = JOB_PARAM

    jt.remoteCommand = './combinePosteriors.py'

    jt.args = [pickled_ids, pickled_idc, pickled_sd, regcov_root, tfbsstats_root, peakstats_root, plotdir, str(fraglen), str(minpost)]

    jobs = s.runBulkJobs(jt,1,count,1)

    print count
    print 'submitted', jobs[1]

    s.synchronize(jobs, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)

    s.deleteJobTemplate(jt)
    s.exit()


def combineFiles(peakstats_root, tfbsstats_root, count, peakstats, TFBSstats):

    op = open(peakstats, 'w')
    ot = open(TFBSstats, 'w')

    for i in arange(1, count+1, 1):
        for pl in open(peakstats_root + '%i' %i):
            op.write(pl)

        for tl in open(tfbsstats_root + '%i' %i):
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
    regcov_dir = cf.get_input("RegCov_dir")
    WM = cf.get_input("WM") 
    statsfile = cf.get_input("statsfile")

    peakstats = cf.get_output("peakstats")
    TFBSstats = cf.get_output("TFBSstats")
    interm = cf.get_output("intermediate")
    log_file = cf.get_output("log_file")
    plotdir = cf.get_output("peak_plots")

    genome = cf.get_parameter("genome", "string")
    minpost = cf.get_parameter("minposterior", "float")
    motevo_path = cf.get_parameter("motevo_path", "string")
    read_files = cf.get_parameter("read_files", "string")
    fpj = cf.get_parameter("files_per_job", "float")
    project_leader = cf.get_parameter("project_leader", "string")

    T1 = datetime.datetime.now()

    ##Main function
    os.mkdir(interm)
    os.mkdir(plotdir)

    fraglen = findFraglen(read_files)

    wmlen = len(open(WM).readlines())-4

    #get parameter file and predicted sites for best WM
    (params, sites) = giveMotevoParamFile(genome, wmlen, interm)
    runMotevo(motevo_path, regions, params, WM, interm, project_leader)

    T2 = datetime.datetime.now()

    #get true sites and a post dict that contains non-refined coordinates as keys and lists of all posteriors as values
    sitesDict, IDstats, IDcoords = getDicts(sites, statsfile, regcov_dir)

    #pickle dictionaries
    pickled_sd = os.path.join(interm, 'sitesDict')
    pickle.dump(sitesDict, pickled_sd)

    pickled_ids = os.path.join(interm, 'IDstats')
    pickle.dump(IDstats, pickled_ids)

    pickled_idc = os.path.join(interm, 'IDcoords')
    pickle.dump(IDcoords, pickled_idc)


    T3 = datetime.datetime.now()

    #create Plots and build peakstats and TFBSstats files
    regcov_root, count = splitRegCovs(regcov_dir, interm, fpj)

    peakstats_root = os.path.join(interm, 'peakstatsfile')
    tfbsstats_root = os.path.join(interm, 'tfbsstatsfile')

    submitJobs(regcov_root, pickled_sd, pickled_ids, pickled_idc, plotdir, fraglen, minpost, peakstats_root, tfbsstats_root, count, project_leader)

    combineFiles(peakstats_root, tfbsstats_root, count, peakstats, TFBSstats)


    T4 = datetime.datetime.now()


    timetext = '\n'.join(['Running time:',
                          '\t-Predicting sites on given regions: %s' %(T2-T1),
                          '\t-Loading data into dictionaries: %s' %(T3-T2),
                          '\t-Plotting coverage profiles with TFBSs and computing peak posteriors: %s' %(T4-T3),
                          ])

    lf = open(log_file, 'w')
    lf.write(text + '\n')
    lf.write(timetext)
    lf.close()


    return 0


component_skeleton.main.main(execute)
                                                                 
