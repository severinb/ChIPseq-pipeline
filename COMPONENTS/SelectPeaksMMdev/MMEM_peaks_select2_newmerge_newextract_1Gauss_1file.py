#!/import/bc2/home/nimwegen/GROUP/local/bin/python

import sys
import os
from pylab import *
from scipy.stats import *
from string import *
from scipy import interpolate

def EM2(x,rho,mu,sigma,maxDiff=0.01):

    logLik = []

    if len(rho) != len(mu) or len(mu) != len(sigma):
        print 'rho, sigma or mu not right\n'
        sys.exit(0)

    mix = len(rho)

    print '------------------------'

    while True:

        # calculate likelihood of x
        pG = [norm.pdf(x, mu[i], sigma[i]) for i in arange(mix)] #list of arrays
        pE = uniform.pdf(x,loc=min(x), scale=max(x))

        logLik.append(sum(log( sum([pG[i]*rho[i] for i in arange(mix)], axis=0) + (1-sum(rho))*pE)))
                
        # check for convergence
        if len(logLik)>1:
            print logLik[-1] - logLik[-2]
        if len(logLik)>1 and (logLik[-1] - logLik[-2]) < maxDiff:
            break
        # calculate probabilities
        ppG = [rho[i] * pG[i] /( sum([pG[j]*rho[j] for j in arange(mix)], axis=0) + (1-sum(rho)) * pE ) for i in arange(mix)] #list of arrays of probabilities
        ppE = (1-sum(rho))*pE/( sum([pG[j]*rho[j] for j in arange(mix)], axis=0) + (1-sum(rho)) * pE )
        # update parameter
        mu = [sum(x*ppG[i])/sum(ppG[i]) for i in arange(mix)]
        sigma = [sqrt(sum(ppG[i]*(x-mu[i])**2)/sum(ppG[i])) for i in arange(mix)]
        #sigma = [min( max(10, sqrt(sum(ppG[i]*(x-mu[i])**2)/sum(ppG[i]))) , 80) for i in arange(mix)]
        rho = [sum(ppG[i])/len(x) for i in arange(mix)]


    print logLik[-1]
    print mu
    print sigma
    print rho

    return (rho, mu, sigma)


def extractPeaks(rho, mu, sig, unif, x):
    """
    This function returns Gaussian looking regions from the input region. 
    It also returns the likelihood of the data in the extracted region assuming this Gaussian.
    """

    ##first try to find overlapping Gaussians that actually represent just one peak
    #For this make windows that are 2sigmas around mu
    windows = []
    width = 1
    for i in arange(len(rho)):
        windows.append((rho[i], mu[i]-width*sig[i], mu[i]+width*sig[i]))

    #try to find overlapping windows.
    sortwins = sorted(windows, key=lambda k: k[1])
    overlapping_wins = []

    cstart = None #current start
    cstop = None

    i = 0
    for win in sortwins:

        if not cstart and not cstop:
            cstart = win[1]
            cstop = win[2]
            overlapping_wins.append([win])

        elif win[1] < cstart and win[2] > cstart:
            overlapping_wins[i].append(win)
            cstart = win[1]
            if cstop < win[2]:
                cstop = win[2]

        elif win[1] > cstart and win[1] < cstop:
            overlapping_wins[i].append(win)
            if cstop < win[2]:
                cstop = win[2]

        else:
            cstart = win[1]
            cstop = win[2]
            overlapping_wins.append([win])

            i += 1


    rho_co = 1.0/(len(rho)+4)
    print 'rho_co: ', rho_co
    #merge overlapping windows (just to get a good sigma range for quality assessment):
    sigma_ranges = []
    for wins in overlapping_wins:
        start = 0
        stop = 0
        totrho = 0.0
        #get total rho for normalization
        for win in wins:
            totrho += win[0]

        if totrho < rho_co:
            continue
        #now merge
        for win in wins:
            start += (win[0]/totrho)*win[1]
            stop += (win[0]/totrho)*win[2]
        if not [int(start),int(stop)] in sigma_ranges:
            sigma_ranges.append([int(start), int(stop)])


    print ''
    print sigma_ranges


    #extend to ~two sigma
    for r in sigma_ranges:
        #half = (r[1]-r[0]) * ((width -1)/2)
        half = (r[1]-r[0]) * ((2-width)/3.0)
        r[0] -= half
        r[1] += half

    print sigma_ranges


    #compute likelihoods and heights
    peaks = []
    x = array(x)
    for r in sigma_ranges:
        d1 = x[where(x >=r[0])]
        datapoints = d1[where(d1 <= r[1])]
        print min(datapoints), max(datapoints)
        likeli = sum((1-sum(rho))*unif.pdf(datapoints))
        for i in arange(len(rho)):
            likeli += sum(rho[i] * norm.pdf(datapoints, mu[i], sig[i]))
        likeli /= r[1]-r[0]
        peaks.append([likeli, r[0], r[1]])


    print ''
    print peaks
    return peaks


def DoTheyOverlap(rho, mu, sig):

    ##first try to find overlapping Gaussians that actually represent just one peak
    #For this make windows that are 2sigmas around mu
    windows = []
    width = 1
    for i in arange(len(rho)):
        windows.append((rho[i], mu[i]-width*sig[i], mu[i]+width*sig[i]))

    #try to find overlapping windows.
    sortwins = sorted(windows, key=lambda k: k[1])
    overlapping_wins = []

    cstart = None #current start
    cstop = None

    i = 0
    for win in sortwins:

        if not cstart and not cstop:
            cstart = win[1]
            cstop = win[2]
            overlapping_wins.append([win])

        elif win[1] < cstart and win[2] > cstart:
            overlapping_wins[i].append(win)
            cstart = win[1]
            if cstop < win[2]:
                cstop = win[2]
            return True

        elif win[1] > cstart and win[1] < cstop:
            overlapping_wins[i].append(win)
            if cstop < win[2]:
                cstop = win[2]
            return True

        else:
            cstart = win[1]
            cstop = win[2]
            overlapping_wins.append([win])

            i += 1


    return False

def getHeight(mu, sig, rho, y):
    """
    y is the counts variable, i.e. sum(y) gives the integral under the observed profile.
    The height here is not the height of the mixture model. It's the height of the collapsed Gaussians.    
    This is not the optimal way, but easier and good enough, I guess.
    """
    #integral of observed profile
    N = sum(y)

    heights = []

    for i in arange(len(rho)):
        a = norm.pdf(mu[i], loc=mu[i], scale=sig[i]) #height of Gaussian

        heights.append(rho[i] * a * N)

    return heights



def getRMSD(mu, sig, rho, x, y):
    """
    x is the x-axis, e.g. 1-1000 or so
    y is the normed profile, i.e. counts/sum(counts)
    This function calculates for each mu-sigma pair the deviation from the Gaussian to the observed profile (x).
    """

    #build mixture profile
    gy = zeros(len(x))
    for i in arange(len(rho)):
        gy += rho[i] * norm.pdf(x, mu[i], sig[i])


    #get RMSD for each sigma-mu pair
    
    qs = []
    for i in arange(len(mu)):
        l = int(mu[i] - sig[i])
        r = int(mu[i] + sig[i])

        rmsd = mean(abs(y[l:r]-gy[l:r]))

        qs.append(rmsd)

    return qs



def main(covfile, order, outfile, outplot):

    name = os.path.split(covfile)[1]
    print name

    a = loadtxt(covfile, usecols=[4,5])

    xs = a.T[0]
    counts = a.T[1]

    global minx, maxx 
    minx = min(xs)
    maxx = max(xs)

    #global width
    #width = 1.3

    #reduce histogram (coverage) if total coverage is too high (EM gets slow)                                                                                                                                                                
    if sum(counts) >= 100000:
        ys = counts/sum(counts)
        ys *= 100000
    else:
        ys = counts

    ##produce data counts (like histogram)
    N = len(xs) #data range 1-1000
    x = []

    for i in arange(N):
        for j in arange(ys[i]):
            x.append(xs[i])

    x = array(x)


    ##start EM
    ##First start with order given as argument and then decrement by one.
    ##Decrement as long as there are no longer any overlapping Gaussian. For every peak I want one Gaussian!
    while 1:
        ##produce initial values for EM
        mu = []
        sig = []
        rho = []

        start = (maxx-minx)/(order+2)
        for i in arange(order):
            #mu.append(randint.rvs(minx, maxx))
            print start, maxx
            #mu = linspace(minx, maxx, order)
            #mu.append(randint.rvs(minx, maxx))
            #mu.append(randint.rvs(start, maxx-start))
            mu.append(randint.rvs(minx+start-1, maxx-start+1))
            print mu
            sig.append(randint.rvs(10,100))
            rho.append((1.0/order)-0.01)

        newrho, newmu, newsig = EM2(x, rho, mu, sig)

        overlap = DoTheyOverlap(newrho, newmu, newsig)

        rho = newrho
        mu = newmu
        sig = newsig

        if not overlap or order == 1:
            break
        else:
            order -= 1
            continue


    RMSD = getRMSD(mu, sig, rho, xs, counts/sum(counts))

    height = getHeight(mu, sig, rho, counts)

    print 'RMSD: ', RMSD
    print 'height: ', height

    ##produce plots
    unif = (1-sum(rho)) * uniform.pdf(xs , loc=min(xs), scale = max(xs))
    y = unif

    for i in arange(len(mu)):
        y += rho[i] * norm.pdf(xs, mu[i], sig[i])
        plot(xs, rho[i] * norm.pdf(xs, mu[i], sig[i]), '--')

    
    plot(xs, counts/sum(counts), 'g', linewidth=2)
    plot(xs,y, 'k', linewidth=1.5)

    ##extract peaks
    #peaks = extractPeaks(rho, mu, sig, uniform(loc=min(xs), scale = max(xs)), x)

    #write data
    n = name.split('_')
    chrom = '_'.join(n[:-2])
    start = int(n[-2])

    o = open(outfile, 'a')

    for i in arange(len(mu)):
        start1 = mu[i] - sig[i]
        stop1 = mu[i] + sig[i]
        q = RMSD[i]
        h = height[i]
        print start1, stop1, q, h

        o.write('%s\t%i\t%i\tid\t+\t%.2f\t%s\n' %(chrom, start+start1, start+stop1, h, str(q)))
        plot([start1, start1],[0,max(y)], label = '%i %i %.4f' %(start1, stop1, h))
        plot([stop1,stop1],[0,max(y)])
        legend()

    o.close()

    savefig(outplot)
    close()


if __name__ == "__main__":


    if len(sys.argv) != 6:
        print '\nUsage: python prog.py filesfile outfile plotdir fraglen order\n'
        sys.exit(0)

    filesfileroot = sys.argv[1] #a file containing paths to files. filesfileroot.SGE_TASK_ID is actual file
    outfileroot = sys.argv[2] #files that contains final windows with scores and statistics. outfileroot.SGE_TASK_ID is used.
    plotdir = sys.argv[3] #where to store plots to. Folder is not created, it's just written to
    fraglen = float(sys.argv[4]) #fragment length to constrain sigmas
    order = int(sys.argv[5])

    taskid = os.environ['SGE_TASK_ID']
    print taskid
    filesfile = filesfileroot + '.%s' %taskid
    outfile = outfileroot + '.%s' %taskid

    print filesfile, outfile

    for peak in open(filesfile):
        ppath = peak.strip()
        name = os.path.split(ppath)[1]
        print name

        outplot = os.path.join(plotdir, name + '.png')
        main(ppath, order, outfile, outplot)
