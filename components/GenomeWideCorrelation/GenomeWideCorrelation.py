#!/usr/bin/env python
import component_skeleton.main
import os, re
from string import *
import datetime
from pylab import *
from scipy.stats import gaussian_kde
        
def make_boxplot(xarr, yarr, binnum, outplot, xlab, ylab, slope, intersection, r2):

    #make box plots of posterior distributions stratified by height 
    binedges = linspace(min(array(xarr)), max(array(xarr)), binnum)

    datamatrix = [[] for i in binedges]

    for i in arange(len(yarr)):
        yval = yarr[i]
        xval = xarr[i]

        index = where(binedges <= xval)[0][-1] #include smaller edge into bin
        datamatrix[index].append(yval)

    #print '-------------------'
    #print binedges
    #print datamatrix


    figure()

    #plot linear fit
    x = arange(min(binedges), max(binedges), 0.01)
    #plot(x, slope*x + intersection, label='y = %s * x + %s\nR2: %s' %(slope, intersection, r2))

    #create violin plots on an axis
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
            except Exception:
                print d, p
            m = k.dataset.min() #lower bound of violin
            M = k.dataset.max() #upper bound of violin
            x = arange(m,M,(M-m)/100.) # support for violin

            #try to not plot outliers (more than 3 standard deviations away from mode):
            kdemode = k.evaluate(x).max()
            kdestd = sqrt(k.dataset.var())
            m1 = max(kdemode - 3*kdestd, m)
            M1 = min(kdemode + 3*kdestd, M)
            x = arange(m1,M1,(M1-m1)/100.)

            v = k.evaluate(x) #violin profile (density curve)
            try:
                v = v/v.max()*w #scaling the violin to the available space
            except Exception:
                print v
            fill_betweenx(x,p,v+p,facecolor='y',alpha=0.3)
            fill_betweenx(x,p,-v+p,facecolor='y',alpha=0.3)
    if bp:
        boxplot(datamatrix, positions=map(float, [round(mean([binedges[i], binedges[i+1]]), 2) for i in arange(len(binedges) -1)]), widths = w, sym='') #ones(len(binedges)) + 4)

    xlim([min(binedges), max(binedges)])
    xticks(rotation=30)
    xlabel(xlab)
    ylabel(ylab)
    legend()
    savefig(outplot)
    close()


def execute(cf):
    """
    This component calls the program that computes WM scores over whole genome
    """

    ##Ports and parameters
    indir = cf.get_input("indir") #outdir of PeakCallNewNoiseModel which contains outzvals file (Z values of windows)
    infile = cf.get_input("binnedWMscores") #output file of BinReads that was ran with computeWMScores outfile 

    datafile = cf.get_output("datafile")
    log_file = cf.get_output("log_file")
    correlation_scatter = cf.get_output("correlation_scatter")
    correlation_violin = cf.get_output("correlation_violin")

    T1 = datetime.datetime.now()

    #load binned WM scores into dictionary:
    #chr1    10000   10500   10250   7.4158
    scoreDict = {}

    for line in open(infile):
        if not line.startswith("#"):
            t = line.strip().split()
            ID = '_'.join(t[:3])
            scoreDict[ID] = t[-1]

    T2 = datetime.datetime.now()

    o = open(datafile, 'w')

    for line in open(os.path.join(indir, 'outzvals')):
        if not line.startswith("#"):
            t = line.strip().split()
            ID = '_'.join(t[:3])

            try:
                wmscore = scoreDict[ID]
                o.write('\t'.join([ID, t[-2], wmscore]) + '\n')
            except KeyError:
                continue

    o.close()

    T3 = datetime.datetime.now()

    #sort datafile by Z-score to be able to make linear regression with equal number of regions above and below Z-score cut-off (default=3)
    sortedfile = os.path.join( os.path.split(datafile)[0], 'sortedfile')
    os.system('sort -k2gnr %s > %s' %(datafile, sortedfile))

    a = loadtxt(sortedfile, usecols=[1,2])


    #compute linear regression
    try:
        abovenum = where(a.T[0] >= 3.0)[0][-1] + 1
    except IndexError:
        # if there are no WM-scores above the cut-off (if for example a WM is given to ComputeWMScores consisting just of ones),
        # then the file datafile and sortedfile will be empty. Then there is no sense in computing correlation.

        # produce dummy plots so that no one can complain
        figure()
        xlabel('Z value')
        ylabel('Summed WM score')
        savefig(correlation_scatter)
        close()

        figure()
        xlabel('Z value bins')
        ylabel('Summed WM scores')
        savefig(correlation_violin)
        close()

        os.system('rm %s' %sortedfile)

        lf = open(log_file, 'w')
        lf.write('\n'.join(['No linear fit made, because most likely there were no WM scores above cut-off. Check the output file \'datafile\'',
                            'Running times:',
                            '\t-Loading WM scores: %s' %str(T2-T1),
                            '\t-Writing file with Z values and WM scores: %s' %str(T3-T2),
                            '\t-Overall: %s' %str(T3-T1)
                            ]))
        lf.close()


        return 0


    belownum = len(a.T[0]) - abovenum

    #take an equal number of above and below Z-cutoff windows
    if abovenum <= belownum:
        print 'More windows are below than above Z-score cut-off of 3'

        indxs = arange(belownum) + abovenum
        shuffle(indxs)

        zs = a.T[0][:abovenum] + a.T[0][indxs[:abovenum]]
        ps = a.T[1][:abovenum] + a.T[0][indxs[:abovenum]]

    else:
        print 'More windows are above than below Z-score cut-off of 3'

        indxs = arange(abovenum)
        shuffle(indxs)

        zs = a.T[0][abovenum:] + a.T[0][indxs[:belownum]]
        ps = a.T[1][abovenum:] + a.T[0][indxs[:belownum]]


    A = vstack([ps, np.ones(len(zs))]).T
    b = lstsq(A, zs)
    slope = b[0][0]
    intersection = b[0][1]

    #compute R**2
    sstot = sum((zs - mean(zs))**2)
    sserr = sum((zs - (slope*ps + intersection))**2)

    #coefficient of determination
    r2 = 1 - (sserr/sstot)
    #Pearson correlation coefficient
    r = cov(zs,ps)/( sqrt(var(zs)) * sqrt(var(ps)))

    print 'slope: %s\nintersection: %s\nR2: %s\n' %(slope, intersection, r2)
    print 'correlation coefficient: %s\n' %(r)

    figure()
    plot(a.T[1], a.T[0], '.', rasterized=True)
    xlabel('Summed WM score')
    ylabel('Z score')
    savefig(correlation_scatter)
    close()

    make_boxplot(a.T[1], a.T[0], 30, correlation_violin, 'WM score bins', 'Z scores', slope, intersection, r2)

    T4 = datetime.datetime.now()

    os.system('rm %s' %sortedfile)

    lf = open(log_file, 'w')
    lf.write('\n'.join(['Linear Fit:',
                        '\t-Z-score = %s * WM-score + %s' %(slope, intersection),
                        '\t-R2 = %s' %r2,
                        'Correlation coefficient: %s' %r[0][1],
                        'Running times:',
                        '\t-Loading WM scores: %s' %str(T2-T1),
                        '\t-Writing file with Z values and WM scores: %s' %str(T3-T2),
                        '\t-Plotting correlation plot: %s' %str(T4-T3),
                        '\t-Overall: %s' %str(T4-T1)
                        ]))
    lf.close()

    return 0


component_skeleton.main.main(execute)
                                                                 
