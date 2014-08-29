#!/usr/bin/env python

import component_skeleton.main
import os
import subprocess
import re 
import drmaa


def createScratchDirectory(outfile):
    scratch_dir = os.path.join(os.path.dirname(outfile), "scratch")
    os.system('mkdir %s' % scratch_dir)
    return scratch_dir


def listOfAllWMs(denovoWMs, databaseWMs, scratchDir):
    WmFilename = os.path.join(scratchDir, 'WMs')
    with open(WmFilename, 'w') as WMs:
        for wm in os.listdir(u'%s' % denovoWMs):
            if re.search('^denovo_WM_\d+$', wm):  # to make sure to only include the denovo motifs
                WMs.write(os.path.join(denovoWMs, wm) + '\n')
        for wm in os.listdir(databaseWMs):
            WMs.write(os.path.join(databaseWMs, wm) + '\n')
    return WmFilename


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


def createJobTemplate(TrainingPool, TestSequences, WMs, scratchDir, genome):
    prog = '/import/bc2/home/nimwegen/omidi/Projects/ChIPseq-pipeline/components/EnrichmentScores/calculate_enrichment_score.py'
    jobFileContent = '\n'.join([
        '#! /bin/bash',
        'WMFILE=%s' % WMs,
        'WM=$(sed -n -e "$SGE_TASK_ID p" $WMFILE)',
        'python {prog} -w $WM \\'.format(prog=prog),
        '-t {trainseq} \\'.format(trainseq=TrainingPool),
        '-s {testseq} \\'.format(testseq=TestSequences),        
        '-o {scratch} -g {genome} '.format(scratch=scratchDir, genome=genome),
        ])
    shellFilename = os.path.join(scratchDir, 'command.sh')
    with open(shellFilename, 'w') as outf:
        outf.write(jobFileContent)
    return shellFilename



def init_job_template(jt, path, args, as_bulk_job, scratchDir, jobName):
    stderr = os.path.join(scratchDir, 'stderr')
    stdout = os.path.join(scratchDir, 'stdout')    
    JOB_PARAM = '-q fs_long -P project_nimwegen -e %s -o %s -b y ' % (stderr, stdout)
    JOB_PARAM += '-N %s' % jobName
    env = {'PATH': '/bin:usr/bin:/import/bc2/home/nimwegen/GROUP/local/bin'}
    jt.workingDirectory = drmaa.JobTemplate.HOME_DIRECTORY
    jt.remoteCommand = '/bin/bash'
    jt.args = [path]
    jt.jobEnvironment = env
    jt.nativeSpecification = JOB_PARAM
    return jt


def runningDrmaaJob(job_path, scratchDir, jobName, NUMBER_OF_JOBS=1):
    s=drmaa.Session()
    s.initialize()
    jt=init_job_template(s.createJobTemplate(), job_path, [], True, scratchDir, jobName)
    all_jobids = []
    if NUMBER_OF_JOBS > 1:
        all_jobids = s.runBulkJobs(jt, 1, NUMBER_OF_JOBS, 1)
        s.synchronize(all_jobids, drmaa.Session.TIMEOUT_WAIT_FOREVER, False)
    else:
        all_jobids = s.runJob(jt)
        retval = s.wait(all_jobids, drmaa.Session.TIMEOUT_WAIT_FOREVER)
    s.deleteJobTemplate(jt)        
    s.exit()    
    return all_jobids
    

def concatenateResults(scratchDir, resFilename):
    cmd = ' '.join([
        'cat',
        os.path.join(scratchDir, '*.results'),
        '|',
        'sort -gr -k 2'
        ])
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    with open(resFilename, 'w') as outf:
        for line in proc.stdout:
            outf.write(line)
    return resFilename


def cleaningUpTmpFiles(scratchDir):
    cmd = 'rm -fr %s' % scratchDir
    os.system(cmd)
    return 0
            

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
    outfile = cf.get_output("EnrichmentScores")
    print "Calculating Enrichment Scores for: "
    print DenovoWMs
    print DatabaseWMs
    ## The scratch directory serves as a temporary space for holding files
    scratchDir = createScratchDirectory(outfile)
    ## create a file that lists all the input WMs    
    WMs = listOfAllWMs(DenovoWMs, DatabaseWMs, scratchDir)
    print "Enrichment Scores: There are in total %d WMs" % len(WMs)
    ## create the training pool that contains both real and decoy (shuffled) sequences
    trainingPool = createSequencePool(TrainingInputSequences, TrainingDecoySequences, \
                                      scratchDir, 'trainingPool')
    testPool = createSequencePool(TestSequences, TestDecoySequences, \
                                  scratchDir, 'testPool')
    ## createJobTemplate for the array job (runs for every motif the fitting and enrichment score program)
    jobName = os.path.basename(os.path.dirname(outfile))
    shellCommand = createJobTemplate(trainingPool, testPool, WMs, scratchDir, GENOME)
    runningDrmaaJob(shellCommand, scratchDir, jobName, NUMBER_OF_JOBS=len(WMs))
    ## make the last result file, sorted by the average enrichment score
    concatenateResults(scratchDir, outfile)
    ## cleaning up the scratch directory
    # cleaningUpTmpFiles(scratchDir)
    return 0


component_skeleton.main.main(execute)
