#!/usr/bin/env python

import os
from string import *
from pylab import *
from scipy import interpolate
import component_skeleton.main
import re
import datetime

def plot_peak_overlay(indir, goodlist, badlist, lengthOut, maxheight, overlayplot):

    def center_peaks(fname, d):
        a1 = loadtxt(fname, usecols=[5])
        #m = argmax(a1)
        m = argmax([mean(a1[i:i+3]) for i in arange(len(a1)-2)]) #find top by averaging coverage around i

        ##Center all peaks to 500th bp of a 1000bp region
        diff = 500 - m
        if diff < 0:
            a2 = a1[-diff:]
        else:
            a2 = list(zeros(diff)) + list(a1)
        if len(a2) >= 1000:
            a3 = a2[:1000]
        else:
            a3 = list(a2) + list(zeros(1000-len(a2)))

        x = arange(1,1001,1)
        tck = interpolate.splrep(x, a3, s=500)
        y = interpolate.splev(x, tck, der=0)

        for i in arange(len(x)):
            d[x[i]] += y[i]

        plot(x,y)
 
        return d
    

    ##for the in-frame peaks
    subplot(131)
    title('Top Z, selected')
    ylim([0,maxheight])

    goodSummedDict = {}
    for i in arange(1,1001,1):
        goodSummedDict[i] = 0.0

    for fi in goodlist:
        fname = os.path.join(indir, fi)
        goodSummedDict = center_peaks(fname, goodSummedDict)


    ##for the out of frame peaks
    subplot(132)
    title('Top Z, too ugly')
    ylim([0,maxheight])


    badSummedDict = {}
    for i in arange(1,1001,1):
        badSummedDict[i] = 0.0

    for fi in badlist:
        fname = os.path.join(indir, fi)
        badSummedDict = center_peaks(fname, badSummedDict)


    ##plot top Z-scoring peaks that didn't make the cut-offs
    subplot(133)
    title('Top Z, too long')
    ylim([0,maxheight])

    lengthSummedDict = {}
    for i in arange(1,1001,1):
        lengthSummedDict[i] = 0.0


    for fi in lengthOut:
        fname = os.path.join(indir, fi)
        lengthSummedDict = center_peaks(fname, lengthSummedDict)


    savefig(overlayplot)
    close()

    figure()
    subplot(131)
    title('Top Z, selected')
    #ylim([0,maxheight])
    plot(goodSummedDict.keys(), goodSummedDict.values())

    subplot(132)
    title('Top Z, too ugly')
    #ylim([0,maxheight])
    plot(badSummedDict.keys(), badSummedDict.values())

    subplot(133)
    title('Top Z, too long')
    #ylim([0,maxheight])
    plot(lengthSummedDict.keys(), lengthSummedDict.values())

    savefig(os.path.split(overlayplot)[0]+'/summedOverlay.pdf')
    close()


def refine(x, y, fname, interm):

    tck = interpolate.splrep(x, y, s=5000)
    ySmooth = interpolate.splev(x, tck, der=0)

    top = argmax(ySmooth)
    r = [ySmooth[top]]
    rr = [y[top]]

    for i in arange(top+1,len(ySmooth),1):
        if ySmooth[i] > ySmooth[i-1]:
            break
        else:
            r.append(ySmooth[i])
            rr.append(y[i])
    for i in arange(top-1,-1,-1):
        if ySmooth[i]>ySmooth[i+1]:
            break
        else:
            r = [ySmooth[i]] + r
            rr = [y[i]] + rr

    yref = array(r)
    yrr = array(rr)

    plotfolder = os.path.join(interm, 'smoothPlots')
    plot(arange(len(yref)),yref)
    plot(arange(len(yrr)),yrr)
    savefig('%s.png' %os.path.join(plotfolder, fname))
    close()

    #return the initial curve but refined
    return yrr



def stats(y, fname, interm):

    #find width of curve at half height
    h = max(y)/2
    mu = argmax(y)

    left = False
    try:
        left =  where(y[:mu] < h)[0][-1]
    except IndexError:
        print 'left is out'

    right = False
    try:
        right = where(y[mu:] < h)[0][0] + mu        
    except IndexError:
        print 'right is out'

    if not left and not right:
        print 'peak doesn\'t look good', fname


    if not left:
        left = mu - (right - mu) 
    if not right:
        right = mu + (mu - left)

    width = (right-left)/2.0 #half width

    sd = width/sqrt(2*log(2))

    #try to put mu more stabley to the middle of the real curve
    mum = mean([left,right])
    mu = mean([mum, mum, mu])

    plotfolder = os.path.join(interm, 'widthPlots')
    plot(arange(len(y)),y)
    plot([left,left], [0, 2*h])
    plot([right,right], [0, 2*h])
    savefig('%s.png' %os.path.join(plotfolder, fname))
    close()


    return mu, sd



def getPeakQuality(infile, fraglen, interm):
    """
    standalone main function to call within python
    """

    fname = os.path.split(infile)[1]

    a = loadtxt(infile, usecols=[4,5])
    x = a.T[0]
    y = a.T[1]

    #yref = refine(x, y, fname, interm)
    yref = y
    mu, sd = stats(yref, fname, interm)
    h = max(yref)

    #These two lines could replace the SelectPeaks component. Or the first RefinePeaks step (could give RegionCoverage directly to this component and refine the selected peaks afterwards)
    if sd * sqrt(2*log(2)) > 2*fraglen:
        return 'toolong', h, sd * sqrt(2*log(2))

    ##Now assess quality
    ##Shift gaussian by -/+ 50 around given mu to find the best fit
    best_qual = 1000.0
    best_yR = []
    best_yG = []
    best_mu = 0
    best_l = 0
    best_r = 0
    best_sd = 0

    shift=15
    step=1
    size = fraglen  #half fragment length is the hypothetical half width of a triangle at half height. Add 25 to get a little bit more

    for mui in arange(mu-shift, mu+shift+1 , step):
        for sdi in arange(sd-shift, sd+shift+1, step):
            xf = arange(len(yref))
            gaussi = h*exp(-0.5*(((xf-mui)/sdi)**2))

            #make normed curves:
            yR = yref/max(yref)
            yG = gaussi/max(gaussi)

            top = argmax(yR)

            l = max(mui - size, 0)
            r = min(mui + size, len(yR) - 1)

            qual = mean(abs(yR[l:r] - yG[l:r]))
            if qual < best_qual:
                best_qual = qual
                best_yR = array(yR)
                best_yG = array(yG)
                best_mu = mui
                best_l = l
                best_r = r
                best_sd = sdi

    plotfolder = os.path.join(interm, 'fitPlots')
    b = arange(len(best_yR))
    plot(b, best_yR, label='RMSD=%s\nh=%s' %(best_qual,h))
    plot(b, best_yG)
    plot([best_l, best_l], [0, 1])
    plot([best_r, best_r], [0, 1])
    legend()
    savefig('%s.png' %os.path.join(plotfolder, fname))
    close()

    l = best_sd*sqrt(2*log(2))
    return best_qual, h, l
    

def execute(cf):

    indir = cf.get_input("in_dir") #covfiles directory
    Z_peaks = cf.get_input("Z_scores") #from PeakMerger, used to get Z-values. Use the sorted list that contains Z-values
    outdir = cf.get_output("out_dir") #selected covfiles
    baddir = cf.get_output("bad_dir") #filtered out covfiles
    logfile = cf.get_output("log_file")
    peakscores = cf.get_output("peakscores")
    stats = cf.get_output("statistics") #this replaces the RefinePeak 0.5 step. To get lengths and heights of peaks...
    interm = cf.get_output("intermediate")
    overlayplot = cf.get_output("peakOverlay")
    readfiles = cf.get_parameter("FGfiles_string", "string")
    fraglen = cf.get_parameter("FragmentLength", "int")
    topPeaks = cf.get_parameter("topPeaks", "int") #number of peaks that should be given out based on Z and RMSD combined scores.
    giveout = cf.get_parameter("giveout", "string") #this is not really used any more, because I'm using a cutoff now.
    Q_co = cf.get_parameter("Q_cutoff", "float")   #a cut-off on Q/RMSD value: only peaks with lower score are selectable.


    os.mkdir(outdir)
    os.mkdir(baddir)
    os.mkdir(interm)

    #make folders for all three kinds of plots
    os.mkdir(os.path.join(interm, 'fitPlots'))    #Plot of each peak with fitted Gaussian (made by getPeakQuality)
    os.mkdir(os.path.join(interm, 'widthPlots'))  #Plot of each peak with computed width (made by stats function)
    os.mkdir(os.path.join(interm, 'smoothPlots')) #Plot of each peak with smoothed curve (made by refine function)

    T1 = datetime.datetime.now()

    ##try to find fragment lengths with readfiles
    if fraglen < 0 :
        if readfiles == '':
            print '\nError: No length cut-offs or read files to get maximum length cut-off are given.\n'
            return 1

        else:
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
            print '\nFragment Length: %s' %meanFraglen

        fraglen = meanFraglen


    T2 = datetime.datetime.now()

    ##Make dictionaries
    #Get quality scores (RMSD to Gaussian fit) and sort out too long peaks 
    Qdict = {}
    lengthOut = []  #list for peaknames that didn't pass the length cut-off
    qOutDict = {}  #list for peaknames that didn't pass RMSD cut-off'
    statsDict = {}

    maxheight = 0.0 #variable for highest peak


    for peak in os.listdir(indir):
        ppath = os.path.join(indir, peak)
        q, h, l = getPeakQuality(ppath, fraglen, interm) #ppath is infile, fraglen is for length cut-off and interm for different plots
        statsDict[peak] = (l, h)
        if q == 'toolong':                   #filter out too long peaks
            lengthOut.append(peak)
        elif q > 0.0 and q <= float(Q_co):   #filter out low peak shape quality peaks
            Qdict[peak] = q
        else:
            qOutDict[peak] = q
        if h > maxheight:
            maxheight = h


    T3 = datetime.datetime.now()

    #Write stats to statistics file (Note: this replaces the stats file of the RefinePeaks (0.5) comonent ran before)
    s = open(stats, 'w')
    for peak in statsDict:
        s.write('%s\t%i\t%.2f\n' %(peak, statsDict[peak][0], statsDict[peak][1]))
    s.close()



    #Get dict of Z-scores:
    fullZdict = {}  #all peaks with Z-scores are here
    for line in open(Z_peaks):
        t = line.split()
        peakpos = '_'.join(t[:3])
        fullZdict[peakpos] = float(t[3]) #here put in all peaks



    #take topPeaks from fullZdict (for comparison with actual selection)
    #topZlist = sorted(fullZdict.keys(), key= lambda k : fullZdict[k], reverse=1)[:topPeaks]
    #zqOut = [i for i in topZlist if i in qOut] #put those peaks into a list that were sorted out by Q_cutoff
    #lqOut = [i for i in topZlist if i in lengthOut] #top Z peaks that were sorted out due to length cut-off
    zqOut = []  #contains peaks that would have been in the topPeaks peaks without filtering RMSD
    zlOut = []  #contains peaks that would have been in the topPeaks peaks without filtering length
    goodpeaks = [] #contains the finally given out peaks 
    badpeaks = []  #contains peaks that made RMSD and length cut-offs but have too low Z-score, so that they didn't make it into the goodpeaks list. Should rather be called notgoodenoughpeaks
    i = 0 #assures that topPeaks peaks would be given out
    for peak in sorted(fullZdict.keys(), key= lambda k : fullZdict[k], reverse=1):
        if i >= topPeaks:
            if peak in Qdict.keys():
                badpeaks.append(peak)
                i += 1 #this is actually not necessary, I do not use i anymore 
        else:
            if peak in qOutDict.keys():
                zqOut.append(peak)
            elif peak in lengthOut:
                zlOut.append(peak)
            elif peak in Qdict.keys():
                goodpeaks.append(peak)
                i += 1
            else:  #skip peaks that are not Qdict, qOut and lengthOut. It means that they are not in the input data set (this can occur when I do not give all top Z-scoring peaks as input)
                continue  


    #write Z and Q scores of peaks (sorted by Z scores, thus top topPeaks peaks of this list is given out)
    ps = open(peakscores, 'w')
    for key in sorted(Qdict, key= lambda k: Qdict[k], reverse=0):
        ps.write( key + '\t' + str(fullZdict[key]) + '\t' + str(Qdict[key]) + '\n')

    ps.close()


    ##Create overlay plot
    plot_peak_overlay(indir, goodpeaks, zqOut, zlOut, maxheight, overlayplot)


    ##copy good covfiles to out_dir
    for peak in goodpeaks:
        os.system('cp %s %s' %(os.path.join(indir, peak), os.path.join(outdir, peak)))

    for peak in zqOut+zlOut: #copy sorted out peaks to a directory
        os.system('cp %s %s' %(os.path.join(indir, peak), os.path.join(baddir, peak)))



    #make scatter of Z RMSD score name. This includes Z and Q values for all peaks (except for the ones that didn't make the length cut-off)
    s = open(logfile + '_Z_RMSD', 'w')
    for k in Qdict.keys():
        s.write('\t'.join([k, str(fullZdict[k]), str(Qdict[k])]) + '\n')
    for k in qOutDict.keys():
        s.write('\t'.join([k, str(fullZdict[k]), str(qOutDict[k])]) + '\n')

    s.close()



    #make histogram of Q-values of given out peaks
    qvals = [Qdict[peak] for peak in goodpeaks]
    figure()
    h = hist(qvals, 60)
    savefig(os.path.join(os.path.split(logfile)[0], 'Qhist.pdf'))
    close()


    T4 = datetime.datetime.now()

    ##create log
    l = open(logfile, 'w')
    text = '\n'.join(['%i of %i peaks were sorted out due to Q cut-off of %s' %(len(qOutDict), len(qOutDict) + len(Qdict) + len(lengthOut) , Q_co),
                      '%i of %i peaks were sorted out due to length cut-off of %s (2*fragment length).' %(len(lengthOut), len(qOutDict) + len(Qdict) + len(lengthOut), (2*fraglen)),
                      '%i peaks are given out based on %s-score.' %(len(goodpeaks), giveout.upper()),
                      '%i peaks were not selected.' %(len(badpeaks)),
                      'Z-score of last given out peak: %s' %(fullZdict[goodpeaks[-1]]),
                      'Q-score of last given out peak: %s' %(Qdict[goodpeaks[-1]]),
                      '%s (Q cut-off) and %s (length cut-off) peaks of %s top Z-scoring peaks were not selected.' %(len(zqOut), len(zlOut), len(goodpeaks) + len(zqOut) + len(zlOut)),
                      '\n'
                      ])

    time = '\n'.join(['Running Time:',
                      '\t-Getting peak shape quality: %s' %(T3-T2),
                      '\t-Writing output: %s' %(T4-T3),
                      '\t-Overall: %s' %(T4-T1),
                      '\n'])

    l.write(text)
    l.write(time)
    l.close()


    return 0

component_skeleton.main.main(execute)
