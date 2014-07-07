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

    return meanFraglen


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



def submitJobs(count, infileroot, outfileroot, plotdir, fraglen, order, project_leader, width, instance_name):

    JOB_PARAM = '-q long -e %s/job.stderr -o %s/job.stdout -j n -w n -N %s -cwd -V -b y' %(os.path.split(plotdir)[0], os.path.split(plotdir)[0], instance_name)
    
    s = drmaa.Session()
    s.initialize()

    jt = s.createJobTemplate()
    jt.nativeSpecification = JOB_PARAM

    jt.remoteCommand = './MMEM2.py'

    jt.args = [infileroot, outfileroot, plotdir, str(fraglen), str(order), str(width)]

    jobs = s.runBulkJobs(jt,1,count,1)

    print 'submitted', jobs[0]

    s.synchronize(jobs, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)

    s.deleteJobTemplate(jt)
    s.exit()


def combineFiles(count, fileroot, interm, allpeaks, statsfile, rmsd_co):

    tmpfile = os.path.join(interm, 'tmpfile')
    badqualpeaksfile = os.path.join(os.path.split(interm)[0], 'badqualpeaks')

    o = open(tmpfile, 'w')
    s = open(statsfile, 'w')
    bqf = open(badqualpeaksfile, 'w')

    for i in arange(1, count+1, 1):
        for line in open(fileroot+'.%i' %i):
            t = line.strip().split()
            if float(t[5]) > rmsd_co:
                bqf.write('%s\t%s\t%s\t%s\t%s\t+\n' %(t[0], t[1], t[2], t[3], t[4]))
            else:
                o.write('%s\t%s\t%s\t%s\t%s\t+\n' %(t[0], t[1], t[2], t[3], t[4]))
            s.write('%s\t%s\t%s\t%s\t%s\t%s\n' %(t[3], t[4], t[5], t[6], t[7], t[8]))

    o.close()
    s.close()
    bqf.close()

    os.system('sort -k5gr %s > %s' %(tmpfile, allpeaks))
    os.system('rm %s' %tmpfile)

    return badqualpeaksfile


def execute(cf):

    indir = cf.get_input("in_dir") #covfiles directory


    logfile = cf.get_output("log_file")
    outfile = cf.get_output("outfile") #contains 1000 top peaks (by height). Peaks above rmsd_co were filtered out. sorted by height
    allpeaks = cf.get_output("allpeaks") #contains all peaks that made rmsd_co. sorted by height
    statsfile = cf.get_output("peakstats") #contains statistics of all found peaks. Ones that made rmsd_co and the ones that didn't

    interm = os.path.join(os.path.split(logfile)[0], 'intermediate') 
    plotdir = os.path.join(os.path.split(logfile)[0], 'peak_plots') 


    #plots:
    height_sigma_scatter = cf.get_output("height_sigma_scatter")
    height_rmsd_scatter = cf.get_output("height_rmsd_scatter")
    sigma_rmsd_scatter = cf.get_output("sigma_rmsd_scatter")
    height_revcum = cf.get_output("height_revcum")
    rmsd_hist = cf.get_output("rmsd_hist")
    sigma_hist = cf.get_output("sigma_hist")
    height_hist = cf.get_output("height_hist")

    #parameters:
    readfiles = cf.get_parameter("FGfiles_string", "string")
    fraglen = cf.get_parameter("FragmentLength", "int")
    fpj = cf.get_parameter("files_per_job", "int")
    order = cf.get_parameter("order", "int")
    project_leader = cf.get_parameter("project_leader", "string")
    toppeaks = cf.get_parameter("topPeaks", "int")
    width = cf.get_parameter("widthFactor", "float") #widthFactor * sigma: this is added to each side of mu to define a peak.
    rmsd_co = cf.get_parameter("RMSD_cutoff", "float") #RMSD cut-off

    instance_name = cf.get_metadata("instanceName")

    os.mkdir(interm)
    os.mkdir(plotdir)

    T1 = datetime.datetime.now()

    ##try to find fragment lengths with readfiles
    if fraglen == -1 :
        if readfiles == '':
            print '\nError: No length cut-offs or read files to get maximum length cut-off are given.\n'
            return 1
        else:
            fraglen = findFraglen(readfiles)

    count, filesfileroot = createJobFiles(indir, interm, fpj)

    T2 = datetime.datetime.now()

    statsfileroot = os.path.join(interm, 'stats')

    submitJobs(count, filesfileroot, statsfileroot, plotdir, fraglen, order, project_leader, width, instance_name)

    T3 = datetime.datetime.now()

    badpeaksfile = combineFiles(count, statsfileroot, interm, allpeaks, statsfile, rmsd_co)

    bpnum = 0
    for i in open(badpeaksfile):
        bpnum += 1


    p = open(outfile, 'w')
    for i, line in enumerate(open(allpeaks)):
        t = line.strip().split()
        lastheight = float(t[4])
        if i >= toppeaks:
            break
        else:
            p.write(line)

    p.close()


    #plot stuff
    a = loadtxt(statsfile, usecols=[1,2,4])
    sigmas = a.T[2]
    heights = a.T[0] #add pseudocount
    quals = a.T[1]

    # if there is just one peak called heights, sigmas and quals will be just numbers. Thus test this and make arrays in that case.
    try:
        len(heights)
    except TypeError:
        sigmas = array([sigmas])
        heights = array([heights])
        quals = array([quals])

    figure()
    plot(log10(sorted(heights, reverse=True)), arange(1,len(heights)+1,1), '.', rasterized=True)
    lh = log10(lastheight)
    plot([lh, lh], [0, len(heights)], label='lowest peak used for motif finding: %s' %lastheight)
    xlabel('log10(heights)')
    ylabel('number of peaks with up to peak height')
    legend()
    savefig(height_revcum)
    close()


    figure()
    hist(sigmas, 300, histtype='step')
    title('Peakwidth (sigma) histogram')
    savefig(sigma_hist)
    close()

    figure()
    hist(log10(heights), 300, histtype='step')
    title('Peakheight (log10-space) histogram')
    savefig(height_hist)
    close()


    figure()
    plot(sigmas, heights, '.', rasterized=True)
    xlabel('peak width (sigma)')
    ylabel('peak height')
    savefig(height_sigma_scatter)
    close()

    log10quals = log10(quals)

    figure()
    plot(sigmas, log10quals, '.', rasterized=True)
    xlabel('peak width (sigma)')
    ylabel('log10(peak quality)')
    savefig(sigma_rmsd_scatter)
    close()

    figure()
    plot(heights, log10quals, '.', rasterized=True)
    xlabel('peak height')
    ylabel('log10(peak quality)')
    savefig(height_rmsd_scatter)
    close()

    figure()
    hist(log10quals[isfinite(log10quals)], 300, histtype='step')
    title('log10(Peak quality (RMSD)) histogram')
    savefig(rmsd_hist)
    close()


    text = '\n'.join(['About %i input regions.' %(fpj*count),
                      'Found mean fragment length %.1f. Constrained sigma to %s-%s' %(fraglen, 0.5*fraglen-16, 0.5*fraglen + 94)
                      #'%i peaks were filtered out due to RMSD cut-off of %s' %(bpnum, rmsd_co)
                      #'Running time for:',
                      #'\tFinding fragment length and splitting up files: %s' %str(T2-T1),
                      #'\tFitting mixture models: %s' %str(T3-T2),
                      #'\tOverall: %s' %str(T3-T1)
                      ])

    l=open(logfile, 'w')
    l.write(text)
    l.close()

    # clean up:
    # tar and zip peak_plots
    pwd = os.getcwd()
    os.system('cd %s && tar -czf ../%s.tar.gz . && cd %s' %(plotdir, os.path.split(plotdir)[1], pwd))
    os.system('rm -r %s' %plotdir)

    #remove intermediate dir
    os.system('rm -r %s' %interm)



    return 0

component_skeleton.main.main(execute)
