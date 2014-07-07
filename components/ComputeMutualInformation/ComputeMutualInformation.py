#!/usr/bin/env python
import component_skeleton.main
import os, re
from string import *
import datetime
from pylab import *
import drmaa
import pickle

def submitClassifications(infile, interm, project_leader, jobnum=100):
    """
    This function submits an array job with drmaa to classify WM scores and Z scores by cut-offs
    """
    # pickle infile so that classify.py can load it faster
    #don't do it because otherwise this wrapper component needs too much memory and it can't be submitted to the queue!
    # pickled_infile = os.path.join(interm, 'pickled_infile')
    # ia = loadtxt(infile, usecols=[1,2])
    # pickle.dump(ia, open(pickled_infile, 'w'))

    ##make files with jobnum cut-off pairs
    #find ranges for Z and WM score cut-off
    minz = 1000
    maxz = -1000
    minwm = 1000
    maxwm = -1000

    for l in open(infile):
        t = l.strip().split()
        if float(t[1]) > maxz:
            maxz = float(t[1])
        if float(t[1]) < minz:
            minz = float(t[1])
        if float(t[2]) > maxwm:
            maxwm = float(t[2])
        if float(t[2]) < minwm:
            minwm = float(t[2])
        
    maxz /= 2. #take half because maximum is more likely in these regions
    maxwm /= 2.
    print maxz, maxwm

    #constrain maxima, otherwise it will just take too long sometimes which would be unnecessary.
    maxz = min(maxz, 20)
    maxwm = min(maxwm, 250)
    print maxz, maxwm

    infileroot = os.path.join(interm, 'cutoffs.')

    n = 1
    filenum = 1
    for i in arange(minz, maxz, 1):
        for j in arange(minwm, maxwm, 1):
            if n == 1:
                f = open(infileroot + str(filenum), 'w')
                filenum += 1
            f.write('%s\t%s\n' %(i, j))
            n += 1
            if n == jobnum:
                f.close()
                n = 1

    #to be sure that file handle is really closed
    f.close()


    ##submit jobs
    outfileroot = os.path.join(interm, 'result.')

    JOB_PARAM = '-q fs_long@@x3755 -P %s -e %s/job.stderr -o %s/job.stdout -j n -N CZWM -cwd -V -b y' %(project_leader, os.path.split(interm)[0], os.path.split(interm)[0])

    s = drmaa.Session()
    s.initialize()

    jt = s.createJobTemplate()
    jt.nativeSpecification = JOB_PARAM
    jt.remoteCommand = './classify.py'
    jt.args = [infile, infileroot, outfileroot]

    jobs = s.runBulkJobs(jt, 1, filenum-1, 1)

    print 'submitted', jobs

    s.synchronize(jobs, drmaa.Session.TIMEOUT_WAIT_FOREVER, False)
    for curjob in jobs:
        print 'Collecting job ' + curjob
        retval = s.wait(curjob, drmaa.Session.TIMEOUT_WAIT_FOREVER)
        print 'Job: ' + str(retval.jobId) + ' finished with status ' + str(retval.hasExited)

    s.deleteJobTemplate(jt)
    s.exit()

    ##read in the results
    bestmi = -10 #best mutual information
    bestz = 0 #Z score cut-off at best mutual information
    bestwm = 0 #WM score cut-off at best mutual information
    bestnmi = 0

    results = [] #data matrix
    currz = 10000.0
    for i in arange(1, filenum, 1):
        for line in open(outfileroot+ str(i)):
            t = line.strip().split()
            if float(t[2]) > bestmi:
                bestmi = float(t[2])
                bestz = float(t[0])
                bestwm = float(t[1])
                bestnmi = float(t[3])
            #make a matrix
            if float(t[0]) != currz:
                currz = float(t[0])
                results.append([])
                results[-1].append(float(t[2]))
            else:
                results[-1].append(float(t[2]))


    return array(results), bestmi, bestz, bestwm, bestnmi, arange(minz, maxz, 1), arange(minwm, maxwm, 1)



def execute(cf):
    """
    This component calls the program that computes WM scores over whole genome
    """

    ##Ports and parameters
    infile = cf.get_input("infile") #file containing Z and WM scores

    interm = cf.get_output("intermediate")
    log_file = cf.get_output("log_file")
    contourplot = cf.get_output("contourplot")

    os.mkdir(interm)

    #check whether infile is empty (can happen when there were no WM scores above cut-off)
    i = 0
    for l in open(infile):
        i += 1
        if i > 2:
            break

    if i == 0: #Thus file is empty
        #create dummy output
        figure()
        xlabel('WM-score')
        ylabel('Z-score')
        savefig(contourplot)
        close()

        lf = open(log_file, 'w')
        lf.write('Mutual information can not get computed because input file is empty. This is most likely due to no WM scores above cut-off.\n')
        lf.close()

        return 0
        

    T1 = datetime.datetime.now()

    results, bestmi, bestz, bestwm, bestnmi, zrange, wmrange = submitClassifications(infile, interm, 'project_nimwegen')

    T2 = datetime.datetime.now()

    #pickle.dump(results, open(os.path.join(interm, 'results.mat'), 'w'))

    CS = contour(wmrange, zrange, results)
    clabel(CS, inline=1, fontsize=10)
    xlabel('WM-score')
    ylabel('Z-score')
    title('Mutual Information Optimum: %.3f bits (Z cut-off %.1f and WM score cut-off %.1f)' %(bestnmi, bestz, bestwm))
    savefig(contourplot)
    close()

    T3 = datetime.datetime.now()

    lf = open(log_file, 'w')
    lf.write('\n'.join(['Mutual Information (at Z-score cut-off %s and WM-score cut-off %s): %s percent' %(bestz, bestwm, bestnmi*100),
                        'Non-normalized mutual information: %s bits' %(bestmi),
                        'Running times:',
                        '\t-Computing mutual information: %s' %str(T2-T1),
                        '\t-Plotting: %s' %str(T3-T2),
                        '\t-Overall: %s' %str(T3-T1)
                        ]))
    lf.close()


    return 0


component_skeleton.main.main(execute)
                                                                 
