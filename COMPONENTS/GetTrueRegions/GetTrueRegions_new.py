#!/usr/bin/env python

import component_skeleton.main
import subprocess
import os, re
from string import *
import datetime, time
from pylab import *
from scipy.stats import gaussian_kde
import sys

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

    
def runMotevo(motevo_path, alignments, params, WM, interm):
    """
    runs Motevo
    """
    
    pwd = os.getcwd()
    os.chdir(interm)

    print '\nrun Motevo'
    proc = subprocess.Popen(motevo_path + ' %s %s %s' %(alignments, params, WM),
                            stdout=subprocess.PIPE,
                            stderr= subprocess.PIPE,
                            shell=True
                            )

    stdout_value, stderr_value = proc.communicate()
    print stdout_value
    print stderr_value
    os.chdir(pwd)
    
    if proc.poll() > 0:
        print '\tstderr:', repr(stderr_value.rstrip())
        return -1
    else:
        return 0


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



def getPeakPlotsTFBS(regcov_dir, sitesDict, IDstats, IDcoords, plotdir, fraglen, minpost, peakstats, TFBSstats):
    """
    -This function maps predicted sites inside the provided regions (stored in sitesDict) to peaks inside these regions (peaks are stored in IDstats).
     Peaks in IDstats were called by SelectPeaksMM on the same provided regions.
    -The function writes two files:
       -peakstats: 1 line per peak: all peak statistics (like in IDstats) plus the peak posterior and weighted peak posterior
       -TFBSstats: 1 line per TFBS: location of TFBS, posterior, distance to nearest peak, coverage at TFBS
    -The function plots the coverage profiles plus the predicted sites that are above minpost

    Dictionary formats:
    -sitesDict:
     reg1000004: [(13-22, -, 0.573747), (13-22, - ,0.573747) ...] 
    -IDstats:
     reg1000087: [[4.393, 1.77485636798e-05, 409.430, 34.822], ]
    -IDcoords:
     chr10_121356000_121357250: reg1000053 
    """

    peakheightPost = {} #dictionary that sums up for each peak (key) the posteriors of nearest TFBSs
    peakheightPostNonWeight = {}   #sums up posteriors and has rmsd as key (as peakheightPost dict)

    to = open(TFBSstats, 'w')
    to.write('#chrom\tstart\tend\tID\tdistance\tposterior\tweighted_posterior\tTFBS_coverage\n')

    for fname in os.listdir(regcov_dir):

        a = loadtxt(os.path.join(regcov_dir, fname), usecols=[4,5])

        regionstart = int(fname.split('_')[-2])
        regionend = int(fname.split('_')[-1])
        chrom = '_'.join(fname.split('_')[:-2])

        ##plot all the predicted TFBSs
        regID = IDcoords[fname] #e.g. reg1000012

        #IDstats is a dictionary containing all Gaussian peaks from mixture modeling in one region.
        IDstats_regID = IDstats[regID]
        mus = [i[2] for i in IDstats_regID] #all mus of all Gaussians in that region
        sigmas = [i[3] for i in IDstats_regID] #all sigmas of all Gaussians in that region
        heights = array([i[0] for i in IDstats_regID]) #all heights of all Gaussians in that region
        quals = array([i[1] for i in IDstats_regID]) #all qualities of all Gaussians in that region
        peakids = array([i[4] for i in IDstats_regID]) #all peak IDs of all Gaussians in that region

        #initialize peak height dictionaries so that all peaks (also the one that have no predicted TFBS inside) are accounted for in the statistics later
        for i in arange(len(mus)):
            dictkey = '%s_%s_%s_%s_%s' %(peakids[i], heights[i], quals[i], mus[i], sigmas[i])
            peakheightPost[dictkey] = 0.0
            peakheightPostNonWeight[dictkey] = 0.0


        #Get sites for this region. If there are none, continue
        try:
            sites = sitesDict[regID]
        except KeyError:
            continue


        height = int(max(a.T[1])) + 10 #just used for plotting

        #only plots profiles when there are also sites
        figure()
        plot(a.T[0], a.T[1], label='Coverage')


        for TFBS in sites:
            if not TFBS[-1] >= 0.0001:
                continue
            else:
                s = int(TFBS[0].split('-')[0]) #start of TFBS
                e = int(TFBS[0].split('-')[1])

                possibledists = array(mus) - mean([s,e]) #see which mu - TFBS distance is smallest
                heightweighteddists = possibledists #* (1/heights) #weight distance height such that peaks do not get assigned to some very low peak.
                whichmu = argmin(abs(heightweighteddists))

                if TFBS[-1] >= minpost: #take everything for plots and statistics, but just plot TFBSs that are above minpost
                    plot([s, s], [0, height], label='%i, post %.2f' %(s, TFBS[-1]))
                    plot([e, e], [0, height], label='%i, post %.2f' %(e, TFBS[-1]))
                    plot([mus[whichmu], mean([s,e])], [height, height])

                to.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(chrom, regionstart + s, regionstart + e, peakids[whichmu], possibledists[whichmu], TFBS[-1], TFBS[-1] * (1./(abs(possibledists[whichmu]) + 1)), a.T[1][int(mean([s,e]))]))

                #keep everything in a dictionary that is used to compute peak posteriors
                if abs(possibledists[whichmu]) <= fraglen: #just add those posteriors to the nearest peak that really lie under the peak!
                    #one dictionary with regID_height_rmsd_mu_sigma as key and summed, weighted posterior as value. The regID is just there to make it unique
                    dictkey = '%s_%s_%s_%s_%s' %(peakids[whichmu], heights[whichmu], quals[whichmu], mus[whichmu], sigmas[whichmu])
                    try:
                        peakheightPost[dictkey] += TFBS[-1] * (1./(abs(possibledists[whichmu]) + 1)) #add pseudo count of 1 to not get high summed posteriors
                        peakheightPostNonWeight[dictkey] += TFBS[-1]
                    except KeyError: #this exception should never be taken, since I initialized the dictionary above
                        peakheightPost[dictkey] = TFBS[-1] * (1./(abs(possibledists[whichmu]) + 1))
                        peakheightPostNonWeight[dictkey] = TFBS[-1]


        xlabel('Position')
        ylabel('Coverage')
        title(fname)
        legend(prop={'size':8})
        savefig(os.path.join(plotdir, fname) + '.png')
        close()

    to.close()


    po = open(peakstats, 'w')
    po.write('#chrom\tstart\tend\tpeakID\theight\tquality\tsummed_posterior\tsummed_weighted_posterior\n')

    #print peakheightPost
    peakheights = []
    peakposts = [] #this one contains summed posteriors weighted by distance from peak center
    peakquals = []
    peakposts_nonweight = []
    for key in peakheightPost:
        peakID = key.split('_')[0]
        regID = peakID.split('.')[0]
        #find key for regID in ID dict
        for key1 in IDcoords:
            if IDcoords[key1] == regID:
                coords= key1

        t = coords.split('_')
        chrom = '_'.join(t[:-2])
        start = int(t[-2])
        end = int(t[-1])

        height = float(key.split('_')[1])
        qual = float(key.split('_')[2])
        mu = float(key.split('_')[3])
        sigma = float(key.split('_')[4])

        peakheights.append(float(height))
        peakquals.append(float(qual))
        peakposts.append(peakheightPost[key])
        peakposts_nonweight.append(peakheightPostNonWeight[key])

        po.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(chrom, int(start + mu - sigma), int(start + mu + sigma), peakID, height, qual, peakheightPostNonWeight[key], peakheightPost[key]))


    po.close()



def make_boxplot(xarr, yarr, binnum, outplot, xlab, ylab):

    #make box plots of posterior distributions stratified by height 
    binedges = linspace(min(array(xarr)), max(array(xarr)), binnum)

    datamatrix = [[] for i in binedges]

    for i in arange(len(yarr)):
        yval = yarr[i]
        xval = xarr[i]

        index = where(binedges <= xval)[0][-1] #include smaller edge into bin
        datamatrix[index].append(yval)

    # print '-------------------'
    # print binedges
    #create violin plots on an axis
    figure()
    bp = True
    violin = True

    dist = max(binedges)-min(binedges)
    w = (1.0/(2*binnum))*dist

    if violin:
        for d,p in zip(datamatrix, map(float, [mean([binedges[i], binedges[i+1]]) for i in arange(len(binedges) - 1)])):
            if len(d) <= 1: #gaussian_kde function doesn't accept one value lists and empty lists
                continue
            try:
                k = gaussian_kde(d) #calculates the kernel density
                m = k.dataset.min() #lower bound of violin
                M = k.dataset.max() #upper bound of violin
                x = arange(m,M,(M-m)/100.) # support for violin
                v = k.evaluate(x) #violin profile (density curve)
                v = v/v.max()*w #scaling the violin to the available space
                fill_betweenx(x,p,v+p,facecolor='y',alpha=0.3)
                fill_betweenx(x,p,-v+p,facecolor='y',alpha=0.3)
            except Exception:
                print d, p, 'Couldn\'t make a violin plot'
    if bp:
        boxplot(datamatrix, positions=map(float, [round(mean([binedges[i], binedges[i+1]]), 2) for i in arange(len(binedges) -1)]), widths = w) #ones(len(binedges)) + 4)

    xlim([min(binedges), max(binedges)])
    xticks(rotation=30)
    xlabel(xlab)
    ylabel(ylab)
    savefig(outplot)
    close()



def plotTFBSdist(peakstats, TFBSstats, minpost, histplot, scatterplot):
    """
    This function plots a histogram of the distance of a peak to it's nearest TFBS
    """

    peakDict = {} #peakid: peakcenter (mu)                                                                                                                                                                                                   
    for line in open(peakstats):
        if line.startswith('#'):
            continue
        else:
            t = line.split()
            center = mean([int(t[2]),int(t[1])])
            peakid = t[3]

            peakDict[peakid] = center


    TFBSdict = {}
    postDict = {}
    for line in open(TFBSstats):
        if line.startswith('#'):
            continue
        else:
            t = line.split()
            center = mean([int(t[2]),int(t[1])])
            peakid = t[3]

            if float(t[5]) >= minpost: #only take TFBSs into account that are above minpost and thus assumed to be real
                try:
                    TFBSdict[peakid].append(center)
                    postDict[peakid].append(float(t[5]))
                except KeyError:
                    TFBSdict[peakid] = [center]
                    postDict[peakid] = [float(t[5])]


    dists = []
    posts = []

    noTFBS = 0

    for peakid in peakDict:
        try:
            sites = array(TFBSdict[peakid])
        except KeyError:
            noTFBS += 1
            continue

        sites -= peakDict[peakid] #subtract peak center from site position

        idx = argmin(sqrt(sites**2))
        dists.append(sites[idx]*len(sites))
        posts.append(postDict[peakid][idx])


    print 'len(dists)', len(dists)
    print 'empty peaks:', noTFBS
    print len(dists) + noTFBS

    figure()
    hist(dists, 100)
    savefig(histplot)
    close()

    figure()
    plot(dists, posts, '.')
    savefig(scatterplot)
    close()



def plotStats(peakstats, plotlist, minpost):


    # plotlist = [height_post_scatter, height_weightedpost_scatter, quality_post_scatter, quality_weightedpost_scatter, height_post_violin, height_weightedpost_violin, 
    #             quality_post_violin, quality_weightedpost_violin, post_hist, weightedpost_hist, post_cumulative]

    a = loadtxt(peakstats, usecols = [4,5,6,7], skiprows=1) #4: height, 5: quality, 6: peak posterior, 7: weighted peak posterior
    peakheights = a.T[0]
    peakposts = a.T[3] #this one contains summed posteriors weighted by distance from peak center
    peakquals = a.T[1]
    peakposts_nonweight = a.T[2]


    #make box plots of posterior distributions stratified by height 
    make_boxplot(peakheights, peakposts_nonweight, 20, plotlist[4], 'peak height', 'peak posterior')
    make_boxplot(peakheights, peakposts, 20, plotlist[5], 'peak height', 'weighted peak posterior')

    make_boxplot(log10(array(peakquals)), peakposts_nonweight, 40, plotlist[6], 'log10(peak quality (RMSD))', 'peak posterior')
    make_boxplot(log10(array(peakquals)), peakposts, 40, plotlist[7], 'log10(peak quality (RMSD))', 'weighted peak posterior')


    ##plot height post scatter and also compute PCA
    #compute PCA to get correlation between peak posterior and peak height
    # eigvals, eigvecs = eig(cov(peakheights, peakposts_nonweight))
    # fov1 = max(eigvals)/sum(eigvals)
    # ind = argmax(eigvals)
    # baseslope = eigvecs[:,ind][1] / eigvecs[:,ind][0]

    # #also compute FOV for shuffled vectors to get a feeling for what this number means
    # aheights = array(peakheights)
    # aposts = array(peakposts_nonweight)
    # shuffle(aheights)
    # shuffle(aposts)
    # aeigvals, aeigvecs = eig(cov(aheights, aposts))
    # afov1 = max(aeigvals)/sum(aeigvals)

    # scalefac = mean([max(peakheights), max(peakposts_nonweight)])

    # Compute FOV of first principal component of the variance normalized covariance matrix. 
    # normalized cov-matrix: [[var1/(sig1*sig1), cov/(sig1*sig2)],[cov/(sig1*sig2), var2/(sig2*sig2)]], which results in: [[1, cov/(sig1*sig2)], [cov/(sig1*sig2), 1]]
    # solving the characteristic polynomial yields: lambda = 1 +- sqrt(cov**2/(var1*var2))
    # Thus FOV = (1 + sqrt(cov**2/(var1*var2)) )/2

    covmat = cov(peakheights, peakposts_nonweight)

    fov1 = (1 + sqrt( (covmat[0][1]**2)/(covmat[0][0] * covmat[1][1]) ))/2

    figure()
    plot(peakheights, peakposts_nonweight, '.')
    #plot([0, eigvecs[:,ind][0] * scalefac], [0,eigvecs[:,ind][1] *scalefac], label='slope = %s' %baseslope)
    xlabel('peak height')
    ylabel('peak posterior')
    legend()
    savefig(plotlist[0])
    close()



    figure()
    plot(peakheights, peakposts, '.')
    xlabel('peak height')
    ylabel('weighted peak posterior')
    savefig(plotlist[1])
    close()

    figure()
    plot(log10(array(peakquals)), peakposts_nonweight, '.')
    xlabel('log10(peak quality (RMSD))')
    ylabel('peak posterior')
    savefig(plotlist[2])
    close()

    figure()
    plot(log10(array(peakquals)), peakposts, '.')
    xlabel('log10(peak quality (RMSD))')
    ylabel('weighted peak posterior')
    savefig(plotlist[3])
    close()

    figure()
    hist(peakposts_nonweight, 100)
    xlabel('peak posterior')
    ylabel('number of peaks')
    savefig(plotlist[8])
    close()

    figure()
    hist(peakposts, 100)
    xlabel('weighted peak posterior')
    ylabel('number of peaks')
    savefig(plotlist[9])
    close()

    figure()
    plot(sorted(peakposts_nonweight), arange(1,len(peakposts_nonweight)+1,1))
    plot([minpost, minpost], [1, len(peakposts_nonweight)], label= 'peak posterior cut-off')
    xscale('log')
    yscale('log')
    xlabel('peak posterior')
    ylabel('number of peaks with up to peak posterior')
    legend(loc='lower right')
    savefig(plotlist[10])
    close()


    return fov1


def computeExpectedCoverage(TFBSstats, plotfile):
    """
    -This function plots a violin plot of TFBS posterior versus coverage at that position
    -Secondly it computes:
     sum(c(i)*w(i)) / ( sum(w(i)) * mean(c) )
     This is a measure for how centered to the peaks TFBSs are.
    """

    a = loadtxt(TFBSstats, usecols=[5,7], skiprows=1) #5: posterior, 6: coverage
    posts = a.T[0]
    TFBScovs = a.T[1]
    make_boxplot(TFBScovs, posts, 20, plotfile, 'coverage/height at TFBS', 'posterior of TFBS')

    expcov = sum(a.T[0]*a.T[1])/(sum(a.T[0]) * mean(a.T[1]))

    return expcov


def execute(cf):
    """
    This component gives true regions (determined by a posterior cut-off on TFBS).
    It produces some plots: 
        -histogram of region posteriors (one with summed posteriors and one with maximum TFBS posterior per region)
        -plots peak coverage (from RegionCoverage) plots with TFBSs (above 0.5 posterior cut-off)
    """

    ##Ports and parameters
    regions = cf.get_input("regions") #alignments or sequences of candidate regions
    regcov_dir = cf.get_input("RegCov_dir")
    WM = cf.get_input("WM") 
    statsfile = cf.get_input("statsfile")

    peakstats = cf.get_output("peakstats")
    TFBSstats = cf.get_output("TFBSstats")
    interm = cf.get_output("intermediate")
    log_file = cf.get_output("log_file")
    plotdir = cf.get_output("peak_plots")

    #plots:
    TFBS_peakcenter_dist_hist = cf.get_output("TFBS_peakcenter_dist_hist")
    TFBS_post_peakcenter_dist_scatter = cf.get_output("TFBS_post_peakcenter_dist_scatter")
    height_post_scatter = cf.get_output("height_post_scatter")
    height_weightedpost_scatter = cf.get_output("height_weightedpost_scatter")
    quality_post_scatter = cf.get_output("quality_post_scatter")
    quality_weightedpost_scatter = cf.get_output("quality_weightedpost_scatter")
    height_post_violin = cf.get_output("height_post_violin")
    height_weightedpost_violin = cf.get_output("height_weightedpost_violin")
    quality_post_violin = cf.get_output("quality_post_violin")
    quality_weightedpost_violin = cf.get_output("quality_weightedpost_violin")
    height_posterior_matches_scatter = cf.get_output("height_posterior_matches_scatter")
    TFBSheight_TFBSpost_scatter = cf.get_output("TFBSheight_TFBSpost_scatter")
    post_hist = cf.get_output("post_hist")
    weightedpost_hist = cf.get_output("weightedpost_hist")
    post_cumulative = cf.get_output("post_cumulative")

    plotlist = [height_post_scatter, height_weightedpost_scatter, quality_post_scatter, quality_weightedpost_scatter, height_post_violin,
                height_weightedpost_violin, quality_post_violin, quality_weightedpost_violin, post_hist, weightedpost_hist, post_cumulative]

    genome = cf.get_parameter("genome", "string")
    minpost = cf.get_parameter("minposterior", "float")
    motevo_path = cf.get_parameter("motevo_path", "string")
    aligned = cf.get_parameter("aligned", "boolean")
    read_files = cf.get_parameter("read_files", "string")

    T1 = datetime.datetime.now()

    ##Main function
    os.mkdir(interm)
    os.mkdir(plotdir)

    fraglen = findFraglen(read_files)

    wmlen = len(open(WM).readlines())-4

    #get parameter file and predicted sites for best WM
    (params, sites) = giveMotevoParamFile(genome, wmlen, interm)
    runMotevo(motevo_path, regions, params, WM, interm)

    T2 = datetime.datetime.now()

    #get true sites and a post dict that contains non-refined coordinates as keys and lists of all posteriors as values
    sitesDict, IDstats, IDcoords = getDicts(sites, statsfile, regcov_dir)

    T3 = datetime.datetime.now()

    #create Plots and build peakstats and TFBSstats files
    getPeakPlotsTFBS(regcov_dir, sitesDict, IDstats, IDcoords, plotdir, fraglen, minpost, peakstats, TFBSstats)

    T4 = datetime.datetime.now()

    plotTFBSdist(peakstats, TFBSstats, minpost, TFBS_peakcenter_dist_hist, TFBS_post_peakcenter_dist_scatter)
    fov1 = plotStats(peakstats, plotlist, minpost)

    expcov = computeExpectedCoverage(TFBSstats, TFBSheight_TFBSpost_scatter)

    T5 = datetime.datetime.now()


    #count how many peaks have peak posterior above minpost
    posts = loadtxt(peakstats, usecols=[6], skiprows=1) #load summed posts
    totalnum = len(posts)
    truenum = len(where(posts>= minpost)[0])
    falsenum = totalnum - truenum


    text = '\n'.join(['Overall statistics:',
                      '\t-%i true, %i false out of %i peaks.' %(truenum, falsenum, totalnum),
                      '\t-%.2f percent are true.' %(100*float(truenum)/totalnum),
                      '\t-Cut-off: minimum summed posterior of %.2f' %minpost,
                      '\t-Peak plots contain TFBS of posterior >= %.2f' %minpost,
                      'Statistic for centering of TFBSs at peak centers:',
                      '\t-Expected normalized coverage (sum(post(i)*coverage(i))/(sum(post(i)) * mean(coverage))): %s' %expcov,
                      'Correlation between peak height and peak posterior:',
                      '\t-Fraction of explained variance by first principal component: %s' %fov1
                      ])


    timetext = '\n'.join(['Running time:',
                          '\t-Predicting sites on given regions: %s' %(T2-T1),
                          '\t-Loading data into dictionaries: %s' %(T3-T2),
                          '\t-Plotting coverage profiles with TFBSs and computing peak posteriors: %s' %(T4-T3),
                          '\t-Computing normalized expected coverage: %s' %(T5-T4),
                          '\t-Overall: %s' %(T5-T1)
                          ])

    lf = open(log_file, 'w')
    lf.write(text + '\n')
    lf.write(timetext)
    lf.close()


    return 0


component_skeleton.main.main(execute)
                                                                 
