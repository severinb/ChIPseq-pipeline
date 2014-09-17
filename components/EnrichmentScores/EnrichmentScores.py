#!/usr/bin/env python

import component_skeleton.main
import os
import subprocess
import re 
import drmaa
from math import ceil
import numpy as np
import matplotlib.pyplot as plt
import datetime

NUMBER_OF_COMPUTATION_NODES = 120


def createScratchDirectory(outfile):
    scratch_dir = os.path.join(os.path.dirname(outfile), "scratch")
    os.system('mkdir %s' % scratch_dir)
    return scratch_dir


def listOfAllWMs(WMdirectories, scratchDir):
    WmFilename = os.path.join(scratchDir, 'WMs')
    WMlist = []
    with open(WmFilename, 'w') as WMs:    
        for WMdirectory in WMdirectories:
            if not WMdirectory:
                continue
            for wm in os.listdir(u'%s' % WMdirectory):
                WMs.write( '%s\n' % os.path.join(WMdirectory, wm) )
                WMlist.append( os.path.join(WMdirectory, wm) )
    return WmFilename, WMlist


def createSequencePool(InputSequences, DecoySequences, scratchDir, resFile):
    resFilename = os.path.join(scratchDir, resFile)
    cmd = ' '.join(['cat',
        InputSequences,
        DecoySequences])
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    with open(resFilename, 'w') as outf:
        for line in proc.stdout:
            outf.write(line)
    return resFilename


def createJobTemplate(TrainingPool, TestSequences, WMs, scratchDir, genome, NUMBER_OF_MOTIFS_PER_JOB):
    prog = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'calculate_enrichment_score.py')
    jobFileContent = '\n'.join([
        '#! /bin/bash',
        'WMFILE=%s' % WMs,
        'for (( i=$SGE_TASK_ID ; i<($SGE_TASK_ID+%d) ; i++ ))' %NUMBER_OF_MOTIFS_PER_JOB,
        'do',
        '   printf "Doing jobID: %s\n" "$i"',
        '   WM=$(sed -n -e "$i p" $WMFILE)',
        '   printf "Doing WM: %s\n" "$WM"',
        '   python {prog} -w "$WM" \\'.format(prog=prog),
        '   -t \'{trainseq}\' \\'.format(trainseq=TrainingPool),
        '   -s \'{testseq}\' \\'.format(testseq=TestSequences),        
        '   -o \'{scratch}\' -g \'{genome}\' '.format(scratch=scratchDir, genome=genome),
        'done'
        ])
    shellFilename = os.path.join(scratchDir, 'command.sh')
    with open(shellFilename, 'w') as outf:
        outf.write(jobFileContent)
    return shellFilename



def init_job_template(jt, path, args, as_bulk_job, scratchDir, jobName, queue_type, project_leader):
    stderr = os.path.join(scratchDir, 'stderr')
    stdout = os.path.join(scratchDir, 'stdout')    
    JOB_PARAM = '-q %s -P %s -e %s -o %s -b y ' % (queue_type, project_leader, stderr, stdout)
    JOB_PARAM += '-N %s' % jobName
    env = {'PATH': '/bin:usr/bin:/import/bc2/home/nimwegen/GROUP/local/bin'}
    jt.jobEnvironment = env
    jt.workingDirectory = drmaa.JobTemplate.HOME_DIRECTORY
    jt.remoteCommand = '/bin/bash'
    jt.args = [path]
    jt.nativeSpecification = JOB_PARAM
    return jt


def runningDrmaaJob(job_path, scratchDir, jobName, queue_type, project_leader, NUMBER_OF_MOTIFS_PER_JOB=1, NUMBER_OF_JOBS=1):
    s=drmaa.Session()
    s.initialize()
    jt=init_job_template(s.createJobTemplate(), job_path, [], True, scratchDir, jobName, queue_type, project_leader)
    all_jobids = []
    if NUMBER_OF_JOBS > 1:
        all_jobids = s.runBulkJobs(jt, 1, NUMBER_OF_JOBS, NUMBER_OF_MOTIFS_PER_JOB)
        s.synchronize(all_jobids, drmaa.Session.TIMEOUT_WAIT_FOREVER, False)
    else:
        all_jobids = s.runJob(jt)
        retval = s.wait(all_jobids, drmaa.Session.TIMEOUT_WAIT_FOREVER)
    s.deleteJobTemplate(jt)        
    s.exit()
    return all_jobids
    

def concatenateResults(scratchDir, resFilename, col):
    sortedWMs = []
    files = [os.path.join(scratchDir, f) \
             for f in os.listdir(scratchDir) if re.search('\.results$', f)]
    resFileUnsorted = resFilename + '.unsorted'
    with open(resFileUnsorted, 'w') as outf:
        for a_file in files:
            with open(a_file) as inf:
                outf.write(inf.readline())
    ## adding a header line that can be used by anduril
    with open(resFilename, 'w') as outf:
        outf.write('WM_path\tenrichment_score\tstdev\tLL_ratio\tbeta\tbg_prior\n')
    ## sorting according to the enrichment score column
    cmd = 'sort -gr -k %d %s >> %s' % (col, resFileUnsorted, resFilename)
    os.system(cmd)
    with open(resFilename) as outf:
        header = True
        for line in outf:
            if header:
                header = False
                continue
            sortedWMs.append( line.split() )

    for a_file in files:
        os.system( "rm '%s'" % a_file )
    os.system( "rm '%s'" % resFileUnsorted )

    return sortedWMs


def cleaningUpTmpFiles(scratchDir):
    cmd = "rm -fr '%s'" % scratchDir
    os.system(cmd)
    return 0


def createWMcombinedFile(topWM, restWMs, scratchDir):
    WmFilename = os.path.join(scratchDir, 'WMs')
    with open(WmFilename, 'w') as outf:
        for a_wm in restWMs:
            outf.write('%s %s\n' % (topWM, a_wm[0]))
    return WmFilename


def findTopWMinWMs(WMs, topWM):
    query = topWM.split(' ')[-1]
    for index in xrange(len(WMs)):
        if WMs[index][0] == query:
            return index
    return -1
    

def combinedMotifs(trainingPool, testPool, WMs, jobName, scratchDir, GENOME, NUMBER_OF_MOTIFS_PER_JOB, queue_type, project_leader):
    index = 1
    numberOfForegroundSeq = len([line for line in open(testPool) if re.search('_reg\d+', line)])
    convergence_criterion = np.log(10.)
    print 'Convergence criterion for finding complementary motifs is %f (at least 10 fold increase in log-likelihood)' %convergence_criterion
    topWM = WMs[0][0]
    topEnrichmentScoreFirstRound = float(WMs[0][1])
    top_LL_FirstRound = float(WMs[0][3])
    WMs.remove(WMs[0])
    topEnrichmentScoreSecondRound = 0.
    enrichmentScoresEachRound = [topEnrichmentScoreFirstRound]
    LL_RatioEachRound = [top_LL_FirstRound]
    while True:        
        WmFile = createWMcombinedFile(topWM, WMs, scratchDir)
        shellCommand = createJobTemplate(trainingPool, testPool, WmFile, \
                                         scratchDir, GENOME, NUMBER_OF_MOTIFS_PER_JOB)
        runningDrmaaJob(shellCommand, scratchDir, jobName, queue_type, project_leader, \
                        NUMBER_OF_MOTIFS_PER_JOB=NUMBER_OF_MOTIFS_PER_JOB, \
                        NUMBER_OF_JOBS=len(WMs))
        ## make the last result file, sorted by the average enrichment score
        outfile = os.path.join(os.path.dirname(scratchDir), 'EnrichmentScores_%d' % (index+1))
        sortedWMs = concatenateResults(scratchDir, outfile, index+2)
        topEnrichmentScoreSecondRound = float(sortedWMs[0][index+1])
        top_LL_SecondRound = float(sortedWMs[0][index+3])
        if (top_LL_SecondRound - top_LL_FirstRound) < convergence_criterion: # convergence criterion: exp( nr_of_fg_seqs * enrichment_score_diff ) < 100 --> enrichment_score diff should increase at least by log(100)/nr_of_fg_seqs (0.009 for 00 fg seqs)
            break
        topWM = ' '.join(sortedWMs[0][:(index+1)])
        removeIndex = findTopWMinWMs(WMs, topWM)
        if not removeIndex == -1:
            WMs.remove(WMs[removeIndex])
        else:
            raise Exception        
        enrichmentScoresEachRound.append(topEnrichmentScoreSecondRound)        
        LL_RatioEachRound.append(top_LL_SecondRound)
        topEnrichmentScoreFirstRound = topEnrichmentScoreSecondRound
        top_LL_FirstRound = top_LL_SecondRound
        index += 1
        if index > 8:
            break
    return topWM, enrichmentScoresEachRound, LL_RatioEachRound
        

def execute(cf):
    """
    It receives two set of motifs: the de novo motifs, and the databse motifs.
    For each of these motifs, it fits the parameters background prior and beta.
    Using the fitted parameters, it calculates the average and sd enrichment scores
    for each of the motifs.
    Final result is a file that its line holds the value of prior, beta, mean
    enrichment score, and standard deviation enrichment score for each motif. 
    """
    TrainingInputSequences = cf.get_input("TrainingSequences")
    TrainingDecoySequences = cf.get_input("TrainingDecoySequences")
    TestSequences = cf.get_input("TestSequences")
    TestDecoySequences = cf.get_input("TestDecoySequences")
    DenovoWMs = cf.get_input("DenovoWMs")
    DatabaseWMs = cf.get_parameter("DatabaseWMs")
    GENOME = cf.get_parameter('genome', 'string')
    CombinedMotifs = cf.get_parameter('CombinedMotifs', 'boolean')
    top_wms = cf.get_parameter('top_wms', 'int')
    queue_type = cf.get_parameter('queue_type', 'string')
    project_leader = cf.get_parameter('project_leader', 'string')
    outfile = cf.get_output("EnrichmentScores")
    outfile_tops = cf.get_output("EnrichmentScores_tops") #This file contains the top WMs of the single WM run. The number of top WMs is controlled by the top_wms parameter.
    
    T1 = datetime.datetime.now()

    print "Calculating Enrichment Scores for: "
    print DenovoWMs
    print DatabaseWMs
    ## The scratch directory serves as a temporary space for holding files
    scratchDir = createScratchDirectory(outfile)
    ## create a file that lists all the input WMs    
    WmFile, WMs = listOfAllWMs([DenovoWMs, DatabaseWMs], scratchDir)    
    print "Enrichment Scores: There are in total %d WMs" % len(WMs)
    ## create the training pool that contains both real and decoy (shuffled) sequences

    trainingPool = createSequencePool(TrainingInputSequences, TrainingDecoySequences, \
                                      scratchDir, 'trainingPool')
    testPool = createSequencePool(TestSequences, TestDecoySequences, \
                                  scratchDir, 'testPool')

    ## createJobTemplate for the array job (runs for every motif the fitting and enrichment score program)    
    jobName = os.path.basename(os.path.dirname(outfile))
    NUMBER_OF_MOTIFS_PER_JOB = max(int(ceil(len(WMs) / NUMBER_OF_COMPUTATION_NODES)), 4) #minimum of 4 motifs per job.
    shellCommand = createJobTemplate(trainingPool, testPool, WmFile, scratchDir, GENOME, NUMBER_OF_MOTIFS_PER_JOB)
    runningDrmaaJob(shellCommand, scratchDir, jobName, queue_type, project_leader, \
                    NUMBER_OF_MOTIFS_PER_JOB=NUMBER_OF_MOTIFS_PER_JOB, \
                    NUMBER_OF_JOBS=len(WMs))
    
    ## make the last result file, sorted by the average enrichment score
    # sortedWMs = [wm.split() for wm in \
    #              open('/import/bc2/home/nimwegen/omidi/Projects/ChIPseq-pipeline/example2/OUTPUT/IRF3_FgBg-enrichmentScores_all_motifs/EnrichmentScores')][:50]
    # scratchDir = '/import/bc2/home/nimwegen/omidi/Projects/ChIPseq-pipeline/example2/OUTPUT/IRF3_FgBg-enrichmentScores_all_motifs/scratch'
    sortedWMs = concatenateResults(scratchDir, outfile, 2)
    with open(outfile_tops, 'w') as outf_tops:
        with open(outfile) as outf:
            i = 0
            for line in outf:
                if i > top_wms:
                    break
                outf_tops.write(line)
                i += 1
 
    ## cleaning up the scratch directory
    # cleaningUpTmpFiles(scratchDir)
    if CombinedMotifs:
        topWM, enrichmentScoresEachRound, LL_ratioEachRound = combinedMotifs(trainingPool, testPool, \
                                                         sortedWMs, jobName, scratchDir, GENOME, \
                                                          NUMBER_OF_MOTIFS_PER_JOB, queue_type, project_leader)
        topWM_list = topWM.split()
        with open(outfile, 'w') as outf:
            outf.write('WM_path\tenrichment_score\tLL_ratio\n')
            for i in range(len(topWM_list)):
                outf.write(topWM_list[i] + '\t' + str(enrichmentScoresEachRound[i]) + '\t' + str(LL_ratioEachRound[i]) + '\n')

        plt.plot(range(len(topWM_list)+1), [0] + enrichmentScoresEachRound, 'r-')
        plt.plot(range(len(topWM_list)+1), [0] + enrichmentScoresEachRound, 'ko')
        plt.xlabel("Motifs")
        plt.ylabel("Enrichment Score")
        plt.xticks(range(len(topWM_list)+1), [''] + [os.path.split(n)[1] for n in topWM_list], rotation=45)
        plt.tight_layout()
        plt.savefig(outfile + '.png')

    T2 = datetime.datetime.now()
    print 'Running time: %s' %str(T2-T1)

    return 0


component_skeleton.main.main(execute)
