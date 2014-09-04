#!/usr/bin/env python

import component_skeleton.main
import os
import subprocess
import re 
import drmaa
from math import ceil

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
    prog = '/import/bc2/home/nimwegen/omidi/Projects/ChIPseq-pipeline/components/EnrichmentScores/calculate_enrichment_score.py'
    jobFileContent = '\n'.join([
        '#! /bin/bash',
        'WMFILE=%s' % WMs,
        'for (( i=$SGE_TASK_ID ; i<($SGE_TASK_ID+%d) ; i++ ))' % NUMBER_OF_MOTIFS_PER_JOB,
        'do',
        '   WM=$(sed -n -e "$i p" $WMFILE)',
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


def runningDrmaaJob(job_path, scratchDir, jobName, NUMBER_OF_MOTIFS_PER_JOB=1, NUMBER_OF_JOBS=1):
    s=drmaa.Session()
    s.initialize()
    jt=init_job_template(s.createJobTemplate(), job_path, [], True, scratchDir, jobName)
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
    ## sorting according to the enrichment score column
    cmd = 'sort -gr -k %d %s > %s' % (col, resFileUnsorted, resFilename)
    os.system(cmd)
    with open(resFilename) as outf:
        for line in outf:
            sortedWMs.append( line.split() )
    for a_file in files:
        os.system( 'rm %s' % a_file )
    return sortedWMs


def cleaningUpTmpFiles(scratchDir):
    cmd = 'rm -fr %s' % scratchDir
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
    

def combinedMotifs(trainingPool, testPool, WMs, jobName, scratchDir, GENOME, NUMBER_OF_MOTIFS_PER_JOB):
    index = 1
    topWM = WMs[0][0]
    topEnrichmentScoreFirstRound = float(WMs[0][1])
    WMs.remove(WMs[0])
    topEnrichmentScoreSecondRound = 0.
    enrichmentScoresEachRound = [topEnrichmentScoreFirstRound]
    while True:        
        WmFile = createWMcombinedFile(topWM, WMs, scratchDir)
        shellCommand = createJobTemplate(trainingPool, testPool, WmFile, \
                                         scratchDir, GENOME, NUMBER_OF_MOTIFS_PER_JOB)
        runningDrmaaJob(shellCommand, scratchDir, jobName, \
                        NUMBER_OF_MOTIFS_PER_JOB=NUMBER_OF_MOTIFS_PER_JOB, \
                        NUMBER_OF_JOBS=len(WMs))
        ## make the last result file, sorted by the average enrichment score
        outfile = os.path.join(os.path.dirname(scratchDir), 'EnrichmentScores_%d' % (index+1))
        sortedWMs = concatenateResults(scratchDir, outfile, index+2)
        print sortedWMs
        topEnrichmentScoreSecondRound = float(sortedWMs[0][index+1])
        if (topEnrichmentScoreSecondRound - topEnrichmentScoreFirstRound) < 0.005:
            break
        topWM = ' '.join(sortedWMs[0][:(index+1)])
        removeIndex = findTopWMinWMs(WMs, topWM)
        if not removeIndex == -1:
            WMs.remove(WMs[removeIndex])
        else:
            raise Exception        
        enrichmentScoresEachRound.append(topEnrichmentScoreSecondRound)        
        topEnrichmentScoreFirstRound = topEnrichmentScoreSecondRound
        index += 1
        if index > 8:
            break
    return topWM, enrichmentScoresEachRound
        

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
    outfile = cf.get_output("EnrichmentScores")
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
    NUMBER_OF_MOTIFS_PER_JOB = int(ceil(len(WMs) / NUMBER_OF_COMPUTATION_NODES))
    shellCommand = createJobTemplate(trainingPool, testPool, WmFile, scratchDir, GENOME, NUMBER_OF_MOTIFS_PER_JOB)
    runningDrmaaJob(shellCommand, scratchDir, jobName, \
                    NUMBER_OF_MOTIFS_PER_JOB=NUMBER_OF_MOTIFS_PER_JOB, \
                    NUMBER_OF_JOBS=len(WMs))
    
    ## make the last result file, sorted by the average enrichment score
    # sortedWMs = [wm.split() for wm in \
    #              open('/import/bc2/home/nimwegen/omidi/Projects/ChIPseq-pipeline/example2/OUTPUT/IRF3_FgBg-enrichmentScores_all_motifs/EnrichmentScores')][:50]
    # scratchDir = '/import/bc2/home/nimwegen/omidi/Projects/ChIPseq-pipeline/example2/OUTPUT/IRF3_FgBg-enrichmentScores_all_motifs/scratch'
    sortedWMs = concatenateResults(scratchDir, outfile, 2)    
    ## cleaning up the scratch directory
    # cleaningUpTmpFiles(scratchDir)
    if CombinedMotifs:
        topWM, enrichmentScoresEachRound = combinedMotifs(trainingPool, testPool, \
                                                         sortedWMs, jobName, scratchDir, GENOME, \
                                                          NUMBER_OF_MOTIFS_PER_JOB)
        with open(os.path.join(os.path.dirname(outfile), "topWMs"), 'w') as outf:
            for WmName in topWM.split():
                outf.write(WmName + '\n')
        with open(os.path.join(os.path.dirname(outfile), "EnrichmentScores_all_rounds"), 'w') as outf:
            for score in enrichmentScoresEachRound:
                outf.write(str(score) + '\n')
        
    return 0


component_skeleton.main.main(execute)
