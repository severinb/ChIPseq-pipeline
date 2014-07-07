#!/usr/bin/env python

import os, re
from string import *
from pylab import *
import sys
import drmaa
import pickle

def getPeakPlotsTFBS(regcov_file, sitesDict, IDstats, IDcoords, plotdir, fraglen, minpost, peakstats, TFBSstats):
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

    po = open(peakstats, 'w')
    #po.write('#chrom\tstart\tend\tpeakID\theight\tquality\tsummed_posterior\n')
    to = open(TFBSstats, 'w')
    #to.write('#chrom\tstart\tend\tnearest_peakID\tnearest_peak_distance_posterior\tTFBS_coverage\n')

    # loop over all significant regions, but throw them out when there was no significant peak inside.
    # If there was no site predicted inside a peak, its posterior is 0.
    for fname in open(regcov_file):

        covs = loadtxt(fname.strip(), usecols=[4,5])

        coords = os.path.split(fname.strip())[1]

        regionstart = int(coords.split('_')[-2])
        regionend = int(coords.split('_')[-1])
        chrom = '_'.join(coords.split('_')[:-2])

        regID = IDcoords[coords] #e.g. reg1000012

        #IDstats is a dictionary containing all Gaussian peaks from mixture modeling in one region (but just the regions that had significant peaks in them!).
        #if no entry is found there was no significant peak inside the region. Thus continue and skip this region
        #regID: height, rmsd, mu, sigma, backheight, Z, peakID
        try:
            IDstats_regID = IDstats[regID]
        except KeyError:
            continue

        mus = [i[2] for i in IDstats_regID] #all mus of all Gaussians in that region
        sigmas = [i[3] for i in IDstats_regID] #all sigmas of all Gaussians in that region
        heights = array([i[0] for i in IDstats_regID]) #all heights of all Gaussians in that region
        quals = array([i[1] for i in IDstats_regID]) #all qualities of all Gaussians in that region
        zscores = array([i[5] for i in IDstats_regID])
        peakids = array([i[6] for i in IDstats_regID]) #all peak IDs of all Gaussians in that region


        #Get sites for this region. If there are none, do not continue 
        #because I want the peaks inside that region to get written to peakstats with posterior 0
        try:
            sites = sitesDict[regID]
        except KeyError:
            sites = []

        height = int(max(covs.T[1])) + 10 #just used for plotting

        #only plots profiles when there are also sites
        figure()
        plot(covs.T[0], covs.T[1], 'r', label='Coverage')

        colourlist = ['b', 'g', 'r', 'c', 'm', 'y'] #used for plotting
        ci = 0 # color index

        # To compare the histogram of distances from TFBSs to peak centers to a background,
        # we want to randomly place sites (as many as we have predicted with the WM) into the region and then take the nearest one to each peak.
        # do this 100 times or so per region. So for 3 peaks inside the region we would have 300 distances at the end
        # Be careful: now I am only computing the distance for a site that lies within 1 fragment length of the peak center. TFBSs outside peaks are not taken into account! I think this is not good!

        # loop over peaks and not over TFBS. Could be faster and better, because then a TFBS can be assigned to several peaks if these peaks are overlapping by fragment length...
        # sites: list of lists of all sites inside the region. [[start-stop (both relative to region start), strand, posterior]]
        # First reformat sites to two lists. One contains middle coordinates and the other posteriors.
        # Also plot TFBSs that are above minpost
        site_midcoords = []
        site_posts = []
        for TFBS in sites:
            if not TFBS[-1] >= 0.0001:
                continue
            else:
                s = int(TFBS[0].split('-')[0]) #start of TFBS
                e = int(TFBS[0].split('-')[1])

                midcoord = mean([s,e])

                site_midcoords.append(midcoord)
                site_posts.append(TFBS[-1])

                possibledists = array(mus) - midcoord #see which mu - TFBS distance is smallest.
                possibledists *= -1 #If TFBS is upstream, distance should be negative. Thus multiply by -1
                whichmu = argmin(abs(possibledists))

                to.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(chrom, regionstart + s, regionstart + e, peakids[whichmu], possibledists[whichmu], TFBS[-1], covs.T[1][int(midcoord)]))

                if TFBS[-1] >= minpost: #take everything for plots and statistics, but just plot TFBSs that are above minpost
                    try:
                        colr = colourlist[ci]
                    except IndexError:
                        colr = colourlist[randint(0,6)]

                    bar(s, height*TFBS[-1], e-s, color=colr, linewidth=0, alpha=0.2, label="%i-%i, post %.2f" %(s, e, TFBS[-1])) # bars with height proportional to posterior
                    ci += 1

        site_midcoords = array(site_midcoords)
        site_posts = array(site_posts)

        ### finish plotting
        xlabel('Position')
        ylabel('Coverage')
        title(os.path.split(fname)[1])
        legend(prop={'size':8})
        savefig(os.path.join(plotdir, coords))
        close()
        ###

        for p_i in arange(len(mus)):
            mu = mus[p_i]
            sigma = sigmas[p_i]
            height = heights[p_i]
            qual = quals[p_i]
            zscore = zscores[p_i]
            peakid = peakids[p_i]

            site_dists = abs(site_midcoords - mu)

            # take all sites within fragment length of peak
            peak_sites_idxs = where(site_dists <= fraglen)

            # if there is no site sum will be 0.0
            peak_post = sum(site_posts[peak_sites_idxs])
            po.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(chrom, int(regionstart + mu - sigma), int(regionstart + mu + sigma), peakid, zscore, qual, peak_post))
            

    to.close()
    po.close()


def main():

    """
    This component gives true regions (determined by a posterior cut-off on TFBS).
    It produces some plots: 
        -histogram of region posteriors (one with summed posteriors and one with maximum TFBS posterior per region)
        -plots peak coverage (from RegionCoverage) plots with TFBSs (above 0.5 posterior cut-off)
    """

    # pickled dictionaries
    IDstats = pickle.load(open(sys.argv[1], 'r'))
    IDcoords = pickle.load(open(sys.argv[2], 'r'))
    sitesDict = pickle.load(open(sys.argv[3], 'r'))

    # file containing file paths to process
    regcovfile_root = sys.argv[4]

    tfbs_outfile_root = sys.argv[5]
    peak_outfile_root = sys.argv[6]
    plotdir = sys.argv[7]

    fraglen = float(sys.argv[8])
    minpost = float(sys.argv[9])

    taskid = os.environ['SGE_TASK_ID'] 

    #create Plots and build peakstats and TFBSstats files
    regcov_file = regcovfile_root + '.' + taskid
    peak_outfile = peak_outfile_root + '.' + taskid
    tfbs_outfile = tfbs_outfile_root + '.' + taskid

    getPeakPlotsTFBS(regcov_file, sitesDict, IDstats, IDcoords, plotdir, fraglen, minpost, peak_outfile, tfbs_outfile)


if __name__ == '__main__':
    main()
                                                                 
