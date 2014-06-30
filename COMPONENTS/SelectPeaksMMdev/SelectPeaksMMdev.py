#!/usr/bin/env python

import os, re
from string import *
from pylab import *
import component_skeleton.main
import datetime
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

    return maxFraglen


def createJobFiles(indir, interm, fpj):
    """
    This function creates files (each containing fpj coverage file paths) and naming them file.SGE_TASK_ID (=1-totfiles/fpj) 
    """

    count = 0 #counts the number of files used
    i = 0

    fname = 'covfiles'

    for f in os.listdir(indir):
        if i >= fpj:
            i = 0
            fh.close()
        if i == 0:
            count += 1
            outfile = os.path.join(interm, fname + '.%i' %count)
            fh = open(outfile, 'w')

        name = os.path.join(indir, f)
        fh.write(name + '\n')
        i += 1

    return count, os.path.join(interm, fname)



def submitJobs(count, infileroot, outfileroot, plotdir, fraglen, order, project_leader):

    JOB_PARAM = '-q fs_long@@x3755 -P %s -e %s/job.stderr -o %s/job.stdout -j n -N MMdev -cwd -b y' %(project_leader, os.path.split(infileroot)[0], os.path.split(infileroot)[0])
    
    s = drmaa.Session()
    s.initialize()

    jt = s.createJobTemplate()
    jt.nativeSpecification = JOB_PARAM
    #jt.remoteCommand = '/import/bc2/home/nimwegen/GROUP/local/bin/python /import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/COMPONENTS/COMPONENTS.v5/SelectPeaksMM/MMEM.py'
    #jt.remoteCommand = '/import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/COMPONENTS/COMPONENTS.v6/SelectPeaksMM/MMEM.py'
    jt.remoteCommand = '/import/bc2/home/nimwegen/GROUP/hseq_pipeline/severin/Anduril/Pipeline/COMPONENTS/COMPONENTS.v6/SelectPeaksMMdev/MMEM_peaks_select2_newmerge_newextract_1Gauss_1file.py'
    jt.args = [infileroot, outfileroot, plotdir, str(fraglen), str(order)]

    #jt.joinFiles = True

    jobs = s.runBulkJobs(jt,1,count,1)

    #print 'submitted', jobs

    s.synchronize(jobs, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)

    s.deleteJobTemplate(jt)
    s.exit()


def combineFiles(count, fileroot, interm, allpeaks):

    tmpfile = os.path.join(interm, 'tmpfile')

    o = open(tmpfile, 'w')

    for i in arange(1, count+1, 1):
        f = open(fileroot+'.%i' %i)
        o.write(f.read())
        f.close()

    o.close()

    os.system('sort -k6gnr %s > %s' %(tmpfile, allpeaks))


def execute(cf):

    indir = cf.get_input("in_dir") #covfiles directory

    outdir = cf.get_output("out_dir") #selected covfiles
    baddir = cf.get_output("bad_dir") #filtered out covfiles
    plotdir = cf.get_output("peak_plots")

    logfile = cf.get_output("log_file")
    peakscores = cf.get_output("peakscores")
    allpeaks = cf.get_output("allpeaks")
    interm = cf.get_output("intermediate")

    readfiles = cf.get_parameter("FGfiles_string", "string")
    fraglen = cf.get_parameter("FragmentLength", "int")
    fpj = cf.get_parameter("files_per_job", "int")
    order = cf.get_parameter("order", "int")
    project_leader = cf.get_parameter("project_leader", "string")
    toppeaks = cf.get_parameter("topPeaks", "int")

    os.mkdir(outdir)
    os.mkdir(baddir)
    os.mkdir(interm)
    os.mkdir(plotdir)


    T1 = datetime.datetime.now()

    ##try to find fragment lengths with readfiles
    if fraglen < 0 :
        if readfiles == '':
            print '\nError: No length cut-offs or read files to get maximum length cut-off are given.\n'
            return 1
        else:
            fraglen = findFraglen(readfiles)

    count, filesfileroot = createJobFiles(indir, interm, fpj)

    T2 = datetime.datetime.now()

    statsfileroot = os.path.join(interm, 'stats')

    submitJobs(count, filesfileroot, statsfileroot, plotdir, fraglen, order, project_leader)

    T3 = datetime.datetime.now()

    combineFiles(count, statsfileroot, interm, allpeaks)


    p = open(peakscores, 'w')
    for i, line in enumerate(open(allpeaks)):
        if i >= toppeaks:
            t = line.strip().split()
            lastheight = float(t[5])
            break
        else:
            p.write(line)

    p.close()

    #plot stuff
    a = loadtxt(allpeaks, usecols=[1,2,5,6])
    sigmas = (a.T[1]-a.T[0])/2.
    heights = a.T[2]+0.01 #add pseudocount
    quals = a.T[3]

    figure()
    plot(log10(heights), arange(1,len(heights)+1,1), '.')
    lh = log10(lastheight)
    plot([lh, lh], [0, len(heights)], label='lowest given out peak: %s' %lastheight)
    xlabel('log10(heights)')
    ylabel('number of peaks with up to peak height')
    legend()
    savefig(os.path.join(os.path.split(interm)[0], 'height_revcum.pdf'))
    close()


    figure()
    hist(sigmas, 50)
    savefig(os.path.join(os.path.split(interm)[0], 'sigmahist.pdf'))
    close()

    figure()
    hist(heights, 50)
    savefig(os.path.join(os.path.split(interm)[0], 'heightshist.pdf'))
    close()


    figure()
    plot(sigmas, heights, '.')
    xlabel('sigmas')
    ylabel('heights')
    savefig(os.path.join(os.path.split(interm)[0], 'sigmaheights.pdf'))
    close()

    figure()
    plot(sigmas, quals, '.')
    xlabel('sigmas')
    ylabel('quals')
    savefig(os.path.join(os.path.split(interm)[0], 'sigmaquals.pdf'))
    close()

    figure()
    plot(heights, quals, '.')
    xlabel('heights')
    ylabel('quals')
    savefig(os.path.join(os.path.split(interm)[0], 'heightsquals.pdf'))
    close()

    figure()
    hist(quals, 50)
    savefig(os.path.join(os.path.split(interm)[0], 'qualhist.pdf'))
    close()

    text = '\n'.join(['About %i input files. %i files per job submitted (%i jobs).' %(fpj*count, fpj, count),
                      'Found mean fragment length %.1f. Constrained sigma to X-X' %fraglen,
                      'Running time for:',
                      '\tFinding fragment length and splitting up files: %s' %str(T2-T1),
                      '\tFitting mixture models: %s' %str(T3-T2),
                      '\tOverall: %s' %str(T3-T1)])

    l=open(logfile, 'w')
    l.write(text)
    l.close()


    return 0

component_skeleton.main.main(execute)
