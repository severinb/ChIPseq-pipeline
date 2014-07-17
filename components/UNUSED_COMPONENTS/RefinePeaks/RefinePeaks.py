#!/usr/bin/env python
import sys
import os
#import Gnuplot
from string import *
from pylab import *
import component_skeleton.main


def refinePeak(peak_file, ratio, plots_dir, outfile, statf):

    peak_dict = {}           #key: bp number, value: coverage
    f = open(peak_file,'r')
    total = 0                #total coverage of input peak
    bps = 0                  #run over all bp from 1 to region length
    pos = 0                  #position of highest covered bp
    coverage = 0             #variable for highest covered bp
    peak = []

    while True:
        line = f.readline()
        if line:
            bps += 1
            peak = line.split()
            peak_pos = int(peak[-2])
            peak_coverage = float(peak[-1])
            peak_dict[peak_pos] = peak_coverage
            total += float(peak[-1])
            if peak_coverage > coverage:    #To find highest covered bp
                coverage = peak_coverage    
                pos = peak_pos              #pos contains position of highest covered bp
        else:
            break

    if not peak:
        return

    chrom = peak[0]
    peak_start = int(peak[1])  #start position of non refined peak
    bp_coverage_mean = int(total/bps)
    try:
        height = (peak_dict[pos] + peak_dict[pos-1] + peak_dict[pos+1])/3. # take an average fot the height to be more robust
    except KeyError:
        height = peak_dict[pos]
    #cutoff = bp_coverage_mean * ratio
    cutoff = height * ratio  #This would be for taking the cut off from the highest covered bp
    
    f.close()

    i = 0
    stop = False
    window = [max(pos-1, 1),min(pos+1, bps)]  #start and end position of window/refined peak. Make sure that it isn't longer than the non refined region 

    while not stop:  #extend peak upstream
        i += 1
        try:
            if (peak_dict[pos-i] + peak_dict[pos-i-1] + peak_dict[pos-i+1])/3. > cutoff:  #use an average to be more robust about local variations
                window[0] = pos-i
            else:
                stop = True
                break
        except KeyError:
            stop = True
            break

    i = 0
    stop = False

    while not stop:
        i += 1
        try:
            if (peak_dict[pos+i] + peak_dict[pos+i-1] + peak_dict[pos+i+1])/3. > cutoff:
                window[1] = pos+i
            else:
                stop = True
                break
        except KeyError:
            stop = True
            break

    out = open(outfile, 'a')
    out.write('\t'.join([chrom, str(peak_start + window[0] -1), str(peak_start + window[1]) , '+']) + '\n')
    out.close()

    stat = '\t'.join([os.path.split(peak_file)[1], str(window[1] - window[0]), str(height)])
    statf.write(stat + '\n')

    peakfilename = os.path.split(peak_file)[1]
    peakplotpath = os.path.join(plots_dir, peakfilename)

    a = loadtxt(peak_file, usecols=[4,5])
    figure()
    plot(a.T[0], a.T[1])
    plot([window[0]]*int(height), range(int(height)), label=str(window[0]))
    plot([window[1]]*int(height), range(int(height)), label=str(window[1]))
    xlabel('Position')
    ylabel('Coverage')
    title(peakfilename)
    legend()
    savefig(peakplotpath + '.png')
    close()

    #g = Gnuplot.Gnuplot()
    #g('set parametric')
    #g('set trange [0:' + str(height) + ']') 
    #g('set term png')
    #g('set output \'' + peakplotpath + '.png\'')
    #g('plot \'' + peak_file + '\' u 5:6 w l, ' + str(window[0]) + ',t ,' + str(window[1]) + ',t')
    #g('plot \'' + peak_file + '\' u 5:6 w l smooth bezier, ' + str(window[0]) + ',t ,' + str(window[1]) + ',t')

def execute(cf):
    """Runs a peakrefiner that looks for base pair with highest coverage and adds neighbouring bps when they are above a cutoff.
       Outputs:
       filename_windows : file containing all windows with summed coverage (coverage=how many reads cover each basepair).
       filename.png : image of peak region with found peak
       refined_peaks : a file with the new peak regions to input for motevo or alignment... chr start end strand(+)
    """


    indir = cf.get_input("in_dir")                #directory from RegionCoverage with 500 files of regions with their coverage
    peak_plots_dir = cf.get_output("peak_plots")
    outfile = cf.get_output("refined_peaks")
    stats = cf.get_output("statistics")
    plotfile = cf.get_output("length_height")
    cutoff = cf.get_parameter("cutoff", "float")  #1.5 is default


    os.mkdir(peak_plots_dir)
    out = open(outfile, 'w')
    out.close()
    #stat file contains ID of no refined region peaklength peakheight
    statf = open(stats, 'w')
    
    #peakname = os.path.split(peak_file)[-1]
    
    regionlist = os.listdir(indir)
    regionlist_sorted = sorted(regionlist)
    for region in regionlist_sorted:
        region_file = os.path.join(indir, region)
        #print 'Processing ' + region
        refinePeak(region_file, cutoff, peak_plots_dir, outfile, statf)
    statf.close()

    #create plot
    a = loadtxt(stats, usecols=[1,2])
    figure()
    plot(a.T[0],a.T[1],'.')
    xlabel('peak length')
    ylabel('peak height')
    savefig(plotfile)
    
    figure()
    h = hist(a.T[0], 300)
    xlabel('peak length')
    ylabel('counts')
    savefig(os.path.join(os.path.split(plotfile)[0], 'lengthHist.pdf'))

    return 0

component_skeleton.main.main(execute)
