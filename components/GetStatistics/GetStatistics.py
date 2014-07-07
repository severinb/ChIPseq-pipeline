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
    savefig(outplot.rstrip('.pdf'))
    close()



def plotTFBSdist(peakstats, TFBSstats, region_dict, minpost, histplot, scatterplot):
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

    #TFBSstats:
    ##chrom  start   end     peakID  distance        posterior       TFBS_coverage
    TFBSdict = {}
    postDict = {}
    for line in open(TFBSstats):
        if line.startswith('#'):
            continue
        else:
            t = line.split()
            center = mean([int(t[2]),int(t[1])])
            regionid = t[3].split('.')[0]

            if float(t[5]) >= minpost: #only take TFBSs into account that are above minpost and thus assumed to be real
                try:
                    TFBSdict[regionid].append(center)
                    postDict[regionid].append(float(t[5]))
                except KeyError:
                    TFBSdict[regionid] = [center]
                    postDict[regionid] = [float(t[5])]


    dists = []
    posts = []

    # shuffle site positions and then take the nearest. I am doing this to get a background distribution
    shuffled_dists = []
    shuffled_posts = []

    noTFBS = 0
    sample_size = 1000

    for peakid in peakDict:
        regionid = peakid.split('.')[0]
        try:
            sites = array(TFBSdict[regionid])
            site_posts = array(postDict[regionid])
        except KeyError:
            noTFBS += 1
            continue

        sites -= peakDict[peakid] #subtract peak center from site position


        idx = argmin(abs(sites))
        dists.append(sites[idx])
        posts.append(site_posts[idx])

        # do shuffling sample_size times
        region_coords = region_dict[regionid]
        n = len(sites)
        shuffled_sites = randint(region_coords[0], region_coords[1], (sample_size, n))

        shuffled_sites -= peakDict[peakid]
        idxs = abs(shuffled_sites).argmin(axis=1)

        for i, idx in enumerate(idxs):
            shuffled_dists.append(shuffled_sites[i, idx])
            shuffled_posts.append(site_posts[idx])


    if len(dists) != 0:
        #h = hist(array(dists), sqrt(len(dists)), weights=posts) #, bins=linspace(-100,100,sqrt(len(dists)))
        binnum = int(sqrt(min(sum(posts),len(posts))))
        binrange = linspace(-100,100,max(2,binnum))
        h = hist(array(dists), bins=binrange, weights=posts) #, bins=linspace(-100,100,sqrt(len(dists)))
        close()
        #shuffled_h = hist(array(shuffled_dists), sqrt(len(shuffled_dists)), weights=shuffled_posts)
        shuffled_h = hist(array(shuffled_dists), bins=binrange, weights=shuffled_posts)
        close()
        figure()
        hbase = [mean([h[1][i], h[1][i+1]]) for i in arange(len(h[1])-1)]
        plot(hbase, h[0]/float(len(peakDict)), label='real sites')

        shuffled_hbase = [mean([shuffled_h[1][i], shuffled_h[1][i+1]]) for i in arange(len(shuffled_h[1])-1)]
        plot(shuffled_hbase, shuffled_h[0]/float(len(peakDict)*sample_size), label='shuffled sites') #times sample_size because I shuffled that many times per peak

        xlabel('TFBS - peakcenter offset')
        ylabel('Probability of having a TFBS at offset')
        xlim([-100,100])
        legend()


    savefig(histplot)
    savefig(histplot.rstrip('.pdf'))
    close()

    figure()
    plot(dists, posts, '.')
    savefig(scatterplot)
    close()



def plotStats(peakstats, plotlist, minpost):


    # plotlist = [height_post_scatter, quality_post_scatter, height_post_violin, quality_post_violin, post_hist, post_cumulative]

    a = loadtxt(peakstats, usecols = [4,5,6], skiprows=1) #4: zscore, 5: quality, 6: peak posterior
    peakzscores = a.T[0]
    peakquals = a.T[1]
    peakposts = a.T[2]

    ##plot height post scatter and also compute PCA
    #compute PCA to get correlation between peak posterior and peak height
    # eigvals, eigvecs = eig(cov(peakzscores, peakposts))
    # fov1 = max(eigvals)/sum(eigvals)
    # ind = argmax(eigvals)
    # baseslope = eigvecs[:,ind][1] / eigvecs[:,ind][0]

    # #also compute FOV for shuffled vectors to get a feeling for what this number means
    # aheights = array(peakzscores)
    # aposts = array(peakposts)
    # shuffle(aheights)
    # shuffle(aposts)
    # aeigvals, aeigvecs = eig(cov(aheights, aposts))
    # afov1 = max(aeigvals)/sum(aeigvals)

    # scalefac = mean([max(peakzscores), max(peakposts)])

    # Compute FOV of first principal component of the variance normalized covariance matrix. 
    # normalized cov-matrix: [[var1/(sig1*sig1), cov/(sig1*sig2)],[cov/(sig1*sig2), var2/(sig2*sig2)]], which results in: [[1, cov/(sig1*sig2)], [cov/(sig1*sig2), 1]]
    # solving the characteristic polynomial yields: (1-lambda)**2 = (cov/sig1*sig2)**2 and then lambda = 1 +- sqrt(cov**2/(var1*var2))
    # Thus FOV = (1 + sqrt(cov**2/(var1*var2)) )/2

    covmat = cov(peakzscores, peakposts)
    r_mat = corrcoef(peakzscores, peakposts)

    #fov1 = (1 + sqrt( (covmat[0][1]**2)/(covmat[0][0] * covmat[1][1]) ))/2
    print r_mat, covmat

    pearson_r = r_mat[0][1]
    fov1 = pearson_r

    global shit
    shit = False
    if isnan(pearson_r):
        shit = True

    #make box plots of posterior distributions stratified by height 
    make_boxplot(peakzscores, peakposts, 20, plotlist[2], 'ChIP-Signal (Z-score)', 'Number of Binding Sites', plabel = 'Correlation: %s' %round(fov1,4) )

    make_boxplot(log10(array(peakquals)), peakposts, 40, plotlist[3], 'log10(peak quality (RMSD))', 'number of binding sites')


    figure()
    plot(peakzscores, peakposts, '.', rasterized=True)
    #plot([0, eigvecs[:,ind][0] * scalefac], [0,eigvecs[:,ind][1] *scalefac], label='slope = %s' %baseslope)
    xlabel('peak Z-score')
    ylabel('number of binding sites')
    savefig(plotlist[0])
    close()


    figure()
    plot(log10(array(peakquals)), peakposts, '.', rasterized=True)
    xlabel('log10(peak quality (RMSD))')
    ylabel('peak posterior')
    savefig(plotlist[1])
    close()


    figure()
    hist(peakposts, sqrt(len(peakposts))) #100)
    xlabel('number of binding sites')
    ylabel('number of peaks')
    savefig(plotlist[4])
    close()


    figure()
    plot(sorted(peakposts), arange(1,len(peakposts)+1,1))
    plot([minpost, minpost], [1, len(peakposts)], label= 'number of binding site cut-off')
    xscale('log')
    yscale('log')
    xlabel('number of binding sites')
    ylabel('number of peaks with up to number of binding sites')
    legend(loc='lower right')
    savefig(plotlist[5])
    close()


    return fov1


def computeExpectedCoverage(TFBSstats, plotfile, covs_list):
    """
    -This function plots a violin plot of TFBS posterior versus coverage at that position
    -Secondly it computes:
     sum(c(i)*w(i)) / ( sum(w(i)) * mean(c) )
     This is a measure for how centered to the peaks TFBSs are.
    """

    # make histograms of coverage frequencies at sites and in total.

    a = loadtxt(TFBSstats, usecols=[5,6], skiprows=1) #5: posterior, 6: coverage
    try:
        posts = a.T[0]
        TFBScovs = a.T[1]
        make_boxplot(TFBScovs, posts, 20, plotfile, 'coverage/height at TFBS', 'posterior of TFBS')
    except IndexError:
        posts = array([1])
        TFBScovs = array([1])
        figure()
        savefig(plotfile)
        savefig(plotfile.rstrip('.pdf'))

    expcov = sum(posts*TFBScovs)/(sum(posts) * mean(covs_list)) #mean(TFBScovs))

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
    regcov_dir = cf.get_input("RegCov_dir")

    log_file = cf.get_output("log_file")

    #plots:
    TFBS_peakcenter_dist_hist = cf.get_output("TFBS_peakcenter_dist_hist")
    TFBS_post_peakcenter_dist_scatter = cf.get_output("TFBS_post_peakcenter_dist_scatter")
    zscore_post_scatter = cf.get_output("zscore_post_scatter")
    quality_post_scatter = cf.get_output("quality_post_scatter")
    zscore_post_violin = cf.get_output("zscore_post_violin")
    quality_post_violin = cf.get_output("quality_post_violin")
    TFBSheight_TFBSpost_scatter = cf.get_output("TFBSheight_TFBSpost_scatter")
    post_hist = cf.get_output("post_hist")
    post_cumulative = cf.get_output("post_cumulative")
    cov_hists = cf.get_output("coverage_histograms")

    plotlist = [zscore_post_scatter, quality_post_scatter, zscore_post_violin,
                quality_post_violin, post_hist, post_cumulative]


    minpost = cf.get_parameter("minposterior", "float")


    T1 = datetime.datetime.now()

    # read in region coverage: one dictionary regionid: start-stop and one list with all coverages
    covs_list = []
    region_dict = {}
    for regcov in os.listdir(regcov_dir):
        #chr1    120313250       120313750       reg1013598      1       3.5
        regcov_file = os.path.join(regcov_dir, regcov)
        a = loadtxt(regcov_file, usecols=[5])
        covs_list += list(a)

        with open(regcov_file) as f:
            l = f.readline().strip().split()
            region_dict[l[3]] = [int(l[1]), int(l[2])]


    plotTFBSdist(peakstats, TFBSstats, region_dict, minpost, TFBS_peakcenter_dist_hist, TFBS_post_peakcenter_dist_scatter)

    fov1 = plotStats(peakstats, plotlist, minpost)

    expcov = computeExpectedCoverage(TFBSstats, TFBSheight_TFBSpost_scatter, covs_list)

    # plot coverage histograms (coverage at sites and coverage everywhere)
    #hr = hist(covs_list, sqrt(len(covs_list)), histtype='step', normed=True, label='region coverage')
    site_covs = loadtxt(TFBSstats, usecols=[6], skiprows=1)
    site_posts = loadtxt(TFBSstats, usecols=[5], skiprows=1)
    binnum = int(sqrt(min(sum(site_posts), len(site_posts))))
    bin_range = linspace(min(covs_list), max(covs_list), max(binnum,2))

    if not shit:
        hr = hist(covs_list, bin_range)
        hs = hist(site_covs, bin_range, weights=site_posts)
        close()

        figure()
        hrbase = [mean([hr[1][i], hr[1][i+1]]) for i in arange(len(hr[1])-1)]
        plot(log(hrbase), log(hr[0]/float(len(covs_list))), label='region coverage')
        hsbase = [mean([hs[1][i], hs[1][i+1]]) for i in arange(len(hs[1])-1)]
        #plot(log(hsbase), log(hs[0]/float(len(site_covs))), label='site coverage')
        plot(log(hsbase), log(hs[0]/sum(site_posts)), label='site coverage')

    else:
        figure()

    title('Enrichment at Binding Sites: %s' %round(expcov,4))
    xlabel('log(Coverage)')
    ylabel('log(Frequency)')
    legend()
    savefig(cov_hists)
    savefig(cov_hists.rstrip('.pdf'))
    close()

    #count how many peaks have peak posterior above minpost
    posts = loadtxt(peakstats, usecols=[6], skiprows=1) #load summed posts
    totalnum = len(posts)
    truenum = len(where(posts>= minpost)[0])
    falsenum = totalnum - truenum


    T2 = datetime.datetime.now()


    text = '\n'.join(['Overall statistics:',
                      '\t- Number of true peaks out of total number of peaks: %i/%i' %(truenum, totalnum),
                      '\t- %.2f percent are true.' %(100*float(truenum)/totalnum),
                      '\t- Cut-off: minimum summed posterior of %.2f' %minpost,
                      '\t- Peak plots contain TFBS of posterior >= %.2f' %minpost,
                      'Statistic for centering of TFBSs at peak centers:',
                      '\t- Enrichment at binding sites: %s' %round(expcov,4),
                      'Correlation between peak Z-score and number of binding sites at peak: %s' %round(fov1,4)
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
                                                                 
