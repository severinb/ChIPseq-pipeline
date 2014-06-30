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

    peakheightPost = {} #dictionary that sums up for each peak (key) the posteriors of nearest TFBSs
    peakheightPostNonWeight = {}   #sums up posteriors and has rmsd as key (as peakheightPost dict)

    to = open(TFBSstats, 'w')
    to.write('#chrom\tstart\tend\tID\tdistance\tposterior\tweighted_posterior\tTFBS_coverage\n')


    for fname in open(regcov_file):

        a = loadtxt(fname, usecols=[4,5])

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


def main():

    """
    This component gives true regions (determined by a posterior cut-off on TFBS).
    It produces some plots: 
        -histogram of region posteriors (one with summed posteriors and one with maximum TFBS posterior per region)
        -plots peak coverage (from RegionCoverage) plots with TFBSs (above 0.5 posterior cut-off)
    """

    # pickled dictionaries
    IDstats = pickle.load(sys.argv[1])
    IDcoords = pickle.load(sys.argv[2])
    sitesDict = pickle.load(sys.argv[3])

    # file containing file paths to process
    regcovfile_root = sys.argv[4]

    tfbs_outfile_root = sys.argv[5]
    peak_outfile_root = sys.argv[6]
    plotdir = sys.argv[7]

    fraglen = int(sys.argv[8])
    minpost = float(sys.argv[9])

    taskid = os.environ['SGE_TASK_ID'] 

    #create Plots and build peakstats and TFBSstats files
    regcov_file = regcovfile_root + '.' + taskid
    peak_outfile = peak_outfile_root + '.' + taskid
    tfbs_outfile = tfbs_outfile_root + '.' + taskid

    getPeakPlotsTFBS(regcov_file, sitesDict, IDstats, IDcoords, plotdir, fraglen, minpost, peak_outfile, tfbs_outfile)


    return 0


component_skeleton.main.main(execute)
                                                                 
