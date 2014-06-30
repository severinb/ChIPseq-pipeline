#!/usr/bin/env python

import component_skeleton.main
import subprocess
import os, re
from string import *
import datetime, time
from pylab import *
from scipy.stats import gaussian_kde
import sys


def make_boxplot(xarr, yarr, binnum, outplot, xlab, ylab, plabel=None):

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
            if len(d) <= 1 or sum(d) == 0: #gaussian_kde function doesn't accept one value lists and empty lists
                continue
            try:
                k = gaussian_kde(d) #calculates the kernel density
                m = k.dataset.min() #lower bound of violin
                M = k.dataset.max() #upper bound of violin
                x = arange(m,M,(M-m)/100.) # support for violin

                ## this doesn't work as it should. Box plots and violins get somehow shifted relatively to each other.
                # #try to not plot outliers (more than 3 standard deviations away from mode):
                # kdemode = k.evaluate(x).max()
                # kdestd = sqrt(k.dataset.var())
                # m1 = max(kdemode - 3*kdestd, m)
                # M1 = min(kdemode + 3*kdestd, M)
                # x = arange(m1,M1,(M1-m1)/100.)

                v = k.evaluate(x) #violin profile (density curve)
                v = v/v.max()*w #scaling the violin to the available space
                fill_betweenx(x,p,v+p,facecolor='y',alpha=0.3)
                fill_betweenx(x,p,-v+p,facecolor='y',alpha=0.3)
            except Exception:
                print d, p, 'Couldn\'t make a violin plot'
    if bp:
        boxplot(datamatrix, positions=map(float, [round(mean([binedges[i], binedges[i+1]]), 2) for i in arange(len(binedges) -1)]), widths = w, sym='')

    xlim([min(binedges), max(binedges)])
    xticks(rotation=30)
    xlabel(xlab)
    ylabel(ylab)
    if plabel:
        title(plabel)
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
        dists.append(sites[idx])#*len(sites))
        posts.append(postDict[peakid][idx])


    print len(peakDict)

    if len(dists) != 0:
        h = hist(array(dists), sqrt(len(dists))) #, bins=linspace(-100,100,sqrt(len(dists))))

        close()
        figure()
        hbase = [mean([h[1][i], h[1][i+1]]) for i in arange(len(h[1])-1)]
        plot(hbase, h[0]/float(len(peakDict)))

        xlabel('TFBS - peakcenter offset')
        ylabel('Probability of having a TFBS inside a peak with offset')
        xlim([-100,100])

    savefig(histplot)
    close()

    figure()
    plot(dists, posts, '.')
    savefig(scatterplot)
    close()



def plotStats(peakstats, plotlist, minpost):


    # plotlist = [height_post_scatter, quality_post_scatter, height_post_violin, quality_post_violin, post_hist, post_cumulative]

    a = loadtxt(peakstats, usecols = [4,5,6], skiprows=1) #4: height, 5: quality, 6: peak posterior
    peakheights = a.T[0]
    peakquals = a.T[1]
    peakposts_nonweight = a.T[2]

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
    # solving the characteristic polynomial yields: (1-lambda)**2 = (cov/sig1*sig2)**2 and then lambda = 1 +- sqrt(cov**2/(var1*var2))
    # Thus FOV = (1 + sqrt(cov**2/(var1*var2)) )/2

    covmat = cov(peakheights, peakposts_nonweight)

    fov1 = (1 + sqrt( (covmat[0][1]**2)/(covmat[0][0] * covmat[1][1]) ))/2

    #make box plots of posterior distributions stratified by height 
    make_boxplot(peakheights, peakposts_nonweight, 20, plotlist[2], 'peak Z-score', 'number of binding sites', plabel = 'Fraction of Explained Variance: %s' %round(fov1,4) )

    make_boxplot(log10(array(peakquals)), peakposts_nonweight, 40, plotlist[3], 'log10(peak quality (RMSD))', 'number of binding sites')


    figure()
    plot(peakheights, peakposts_nonweight, '.', rasterized=True)
    #plot([0, eigvecs[:,ind][0] * scalefac], [0,eigvecs[:,ind][1] *scalefac], label='slope = %s' %baseslope)
    xlabel('peak Z-score')
    ylabel('number of binding sites')
    savefig(plotlist[0])
    close()


    figure()
    plot(log10(array(peakquals)), peakposts_nonweight, '.', rasterized=True)
    xlabel('log10(peak quality (RMSD))')
    ylabel('peak posterior')
    savefig(plotlist[1])
    close()


    figure()
    hist(peakposts_nonweight, sqrt(len(peakposts_nonweight))) #100)
    xlabel('number of binding sites')
    ylabel('number of peaks')
    savefig(plotlist[4])
    close()


    figure()
    plot(sorted(peakposts_nonweight), arange(1,len(peakposts_nonweight)+1,1))
    plot([minpost, minpost], [1, len(peakposts_nonweight)], label= 'number of binding site cut-off')
    xscale('log')
    yscale('log')
    xlabel('number of binding sites')
    ylabel('number of peaks with up to number of binding sites')
    legend(loc='lower right')
    savefig(plotlist[5])
    close()


    return fov1


def computeExpectedCoverage(TFBSstats, plotfile):
    """
    -This function plots a violin plot of TFBS posterior versus coverage at that position
    -Secondly it computes:
     sum(c(i)*w(i)) / ( sum(w(i)) * mean(c) )
     This is a measure for how centered to the peaks TFBSs are.
    """

    a = loadtxt(TFBSstats, usecols=[5,6], skiprows=1) #5: posterior, 6: coverage
    try:
        posts = a.T[0]
        TFBScovs = a.T[1]
    except IndexError:
        posts = array([1])
        TFBScovs = array([1])

    make_boxplot(TFBScovs, posts, 20, plotfile, 'coverage/height at TFBS', 'posterior of TFBS')

    expcov = sum(posts*TFBScovs)/(sum(posts) * mean(TFBScovs))

    return expcov


def execute(cf):
    """
    This component gives true regions (determined by a posterior cut-off on TFBS).
    It produces some plots: 
        -histogram of region posteriors (one with summed posteriors and one with maximum TFBS posterior per region)
        -plots peak coverage (from RegionCoverage) plots with TFBSs (above 0.5 posterior cut-off)
    """

    ##Ports and parameters
    peakstats = cf.get_input("peakstats") 
    TFBSstats = cf.get_input("TFBSstats")

    interm = cf.get_output("intermediate")
    log_file = cf.get_output("log_file")

    #plots:
    TFBS_peakcenter_dist_hist = cf.get_output("TFBS_peakcenter_dist_hist")
    TFBS_post_peakcenter_dist_scatter = cf.get_output("TFBS_post_peakcenter_dist_scatter")
    height_post_scatter = cf.get_output("height_post_scatter")
    quality_post_scatter = cf.get_output("quality_post_scatter")
    height_post_violin = cf.get_output("height_post_violin")
    quality_post_violin = cf.get_output("quality_post_violin")
    height_posterior_matches_scatter = cf.get_output("height_posterior_matches_scatter")
    TFBSheight_TFBSpost_scatter = cf.get_output("TFBSheight_TFBSpost_scatter")
    post_hist = cf.get_output("post_hist")
    post_cumulative = cf.get_output("post_cumulative")

    plotlist = [height_post_scatter, quality_post_scatter, height_post_violin,
                quality_post_violin, post_hist, post_cumulative]


    minpost = cf.get_parameter("minposterior", "float")


    T1 = datetime.datetime.now()


    plotTFBSdist(peakstats, TFBSstats, minpost, TFBS_peakcenter_dist_hist, TFBS_post_peakcenter_dist_scatter)

    fov1 = plotStats(peakstats, plotlist, minpost)

    expcov = computeExpectedCoverage(TFBSstats, TFBSheight_TFBSpost_scatter)

    #count how many peaks have peak posterior above minpost
    posts = loadtxt(peakstats, usecols=[6], skiprows=1) #load summed posts
    totalnum = len(posts)
    truenum = len(where(posts>= minpost)[0])
    falsenum = totalnum - truenum


    T2 = datetime.datetime.now()


    text = '\n'.join(['Overall statistics:',
                      '\t- %i true, %i false out of %i peaks.' %(truenum, falsenum, totalnum),
                      '\t- %.2f percent are true.' %(100*float(truenum)/totalnum),
                      '\t- Cut-off: minimum summed posterior of %.2f' %minpost,
                      '\t- Peak plots contain TFBS of posterior >= %.2f' %minpost,
                      'Statistic for centering of TFBSs at peak centers:',
                      '\t- Enrichment at binding sites: %s' %expcov,
                      'Correlation between peak Z-score and peak posterior:',
                      '\t- Fraction of explained variance by first principal component: %s' %fov1
                      ])


    timetext = '\n'.join(['Running time:',
                          '\t- Overall: %s' %(T2-T1)
                          ])

    lf = open(log_file, 'w')
    lf.write(text + '\n')
    #lf.write(timetext)
    lf.close()


    return 0


component_skeleton.main.main(execute)
                                                                 
