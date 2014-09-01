#!/usr/bin/env python
import component_skeleton.main
import os, re
from string import *
import datetime
from pylab import *
import drmaa

# A, C, G, T, N
#BG_DEFAULT = np.array([0.25, 0.25, 0.25, 0.25, 0.25],dtype=float)
BG_DEFAULT = '0.25,0.25,0.25,0.25,0.25'
hg19_DEFAULT_DIR = '/import/bc2/home/nimwegen/geier/data/hg19_genome'
mm9_DEFAULT_DIR = '/import/bc2/home/nimwegen/geier/data/mm9_genome'
dm3_DEFAULT_DIR = '/import/bc2/home/nimwegen/geier/data/dm3_genome'


def makeBedFile(chrominfo, interm, winlen=10000000):
    """
    This function makes a bed file that contains windows of about winlen. This way site predictions can be split into several jobs...
    """

    o = open(os.path.join(interm, 'bedfile'), 'w')

    for line in open(chrominfo):
        chrom, end = line.strip().split()
        start = 1
        stop = False
        while not stop:
            end1 = start + winlen - 1
            if end1 >= int(end):
                stop = True
            if (int(end) - end1) <= 1000000: #this assures that at the end the real end is printed and not end1 and also that not a very short window is printed at the end
                end1 = int(end)
                stop = True
            o.write('\t'.join([chrom, str(start), str(end1)]) + '\n')
            start += winlen

    o.close()

    return os.path.join(interm, 'bedfile')



def predictSites(mat, chrominfo, interm, cutoff, chromdir, bg, outfile, project_leader):

    #give an index (SGE_TASK_ID) to each job which specifies which line to read from chrominfo and how to name the output file
    chromnum = len(open(chrominfo).readlines())
    print 'chromnum: ', chromnum

    JOB_PARAM = '-q fs_long -P %s -e %s/job.stderr -o %s/job.stdout -j n -N GWWMS -cwd -V -b y' %(project_leader, interm, interm)

    s = drmaa.Session()
    s.initialize()

    jt = s.createJobTemplate()
    jt.nativeSpecification = JOB_PARAM
    jt.remoteCommand = './runWM.py'
    print mat, chrominfo, interm, cutoff, chromdir, bg
    jt.args = [mat, chrominfo, interm, str(cutoff), chromdir, bg]

    jobs = s.runBulkJobs(jt,1,chromnum,1)

    print 'submitted', jobs

    s.synchronize(jobs, drmaa.Session.TIMEOUT_WAIT_FOREVER, False)
    for curjob in jobs:
        print 'Collecting job ' + curjob
        retval = s.wait(curjob, drmaa.Session.TIMEOUT_WAIT_FOREVER)
        print 'Job: ' + str(retval.jobId) + ' finished with status ' + str(retval.hasExited)

    s.deleteJobTemplate(jt)
    s.exit()

    o = open(outfile, 'w')
    for i in range(1,chromnum+1,1):
        f = open(os.path.join(interm, 'outfile_' + str(i)))
        for line in f:
            o.write(line)
        f.close()
        os.system('rm %s' %os.path.join(interm, 'outfile_' + str(i)))

    o.close()
        

def execute(cf):
    """
    This component calls the program that computes WM scores over whole genome
    """

    ##Ports and parameters
    mat = cf.get_input("WM")

    interm = cf.get_output("intermediate")
    log_file = cf.get_output("log_file")
    outfile = cf.get_output("WMscores")
    
    chrominfo = cf.get_parameter("chrominfo", "string") #a bed file with chromosome regions to look for binding sites (chr end)
    cutoff = cf.get_parameter("WMscore_cutoff", "float")
    chromdir = cf.get_parameter("hd5f_chromdir", "string") #directory with chromosomes in hdf5 format
    bg = cf.get_parameter("bg_frequencies", "string")
    project_leader = cf.get_parameter("project_leader", "string")


    os.system('echo $LD_LIBRARY_PATH')
    os.mkdir(interm)

    T1 = datetime.datetime.now()

    bedfile = makeBedFile(chrominfo, interm)
    predictSites(mat, bedfile, interm, cutoff, chromdir, bg, outfile, project_leader)

    T2 = datetime.datetime.now()

    lf = open(log_file, 'w')
    lf.write('Running time: %s' %str(T2-T1))
    lf.close()

    return 0


component_skeleton.main.main(execute)
                                                                 
