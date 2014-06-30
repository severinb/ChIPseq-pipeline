#!/import/bc2/home/nimwegen/GROUP/local/bin/python

import os, sys
from string import *
from pylab import *
from scipy.stats import *

def EM2(x, rho, mu, sigma, fraglen, maxDiff=0.01):

    logLik = []

    if len(rho) != len(mu) or len(mu) != len(sigma):
        print 'rho, sigma or mu not right\n'
        sys.exit(0)

    mix = len(rho)

    #the linear relation of maximum fragment length and mode of the sigma distribution of high peaks was 0.3*fraglen+28. So I took a range of +- 30
    #If a negative fragment length is given, sigma is not constrained
    if fraglen < 0:
        minsig = 0
        maxsig = 10000
    else:
        minsig = 0.3*fraglen #- 10
        maxsig = 0.3*fraglen + 60

    print "Allow sigma range of %s-%s" %(minsig, maxsig)

    while True:
        print '------------------------'

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
        #sigma = [sqrt(sum(ppG[i]*(x-mu[i])**2)/sum(ppG[i])) for i in arange(mix)]
        sigma = [min( max(minsig, sqrt(sum(ppG[i]*(x-mu[i])**2)/sum(ppG[i]))) , maxsig) for i in arange(mix)]
        rho = [sum(ppG[i])/len(x) for i in arange(mix)]


    print logLik[-1]
    print mu
    print sigma
    print rho

    return (rho, mu, sigma)


def collapseOverlap(mu, sig, rho, width):
    """
    Overlapping windows are bad for phylogibbs (results in exact window matches).
    One nice peak can also be constructed of several Gaussians with similar mu.
    Therefore combine overlapping Gaussians, add up rhos and extend sigmas.
    """

    windows = []
    for i in arange(len(rho)):
        windows.append((rho[i], mu[i]-width*sig[i], mu[i]+width*sig[i], mu[i], sig[i]))

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

    #now combine overlapping ones. mu and sigmas get added up in a rho-weighted manner
    newrho = []
    newmu = []
    newsig = []

    for wins in overlapping_wins:
        totmu = 0.0
        totsig = 0.0
        totrho = 0.0

        #just take the whole sigma range of all overlapping windows. I want to capture all binding sites.
        #most left resp. right position of overlapping windows
        #l = 100000.0 
        #r = 0.0
        
        #get total rho for normalization
        for win in wins:
            totrho += win[0]

        #now merge
        for win in wins:
            totmu += (win[0]/totrho) * win[-2]
            totsig += (win[0]/totrho) * win[-1]
            #if (win[-2]-win[-1]) < l:
            #    l = (win[-2]-win[-1])
            #if (win[-2]+win[-1]) > r:
            #    r = (win[-2]+win[-1])

        #totsig = (l+r)/2.0

        newrho.append(totrho)
        newmu.append(totmu)
        newsig.append(totsig)


    print newrho
    print newmu
    print newsig

    return newrho, newmu, newsig


def getRMSD(gy, y, Cmu, Csig):
    """
    gy is the mixture model fitted to the profile (y)
    y is the normed profile, i.e. counts/sum(counts)
    This function calculates for each mu-sigma pair the deviation from the Gaussian to the observed profile (x).
    """

    #get RMSD for each sigma-mu pair
    qs = []
    for i in arange(len(Cmu)):
        l = int(Cmu[i] - Csig[i])
        r = int(Cmu[i] + Csig[i])

        rmsd = mean(abs(y[l:r]-gy[l:r]))
        if isnan(rmsd):
            rmsd = 0

        qs.append(rmsd)

    return qs


def getHeight(gy, mu, sig, rho, y):
    """
    y is the counts variable, i.e. sum(y) gives the integral under the observed profile.
    The height here is not the height of the mixture model. It's the height of the collapsed Gaussians.
    This is not the optimal way, but easier and good enough, I guess.
    """

    #integral of observed profile
    N = sum(y)
    x = arange(1, len(y)+1, 1)

    heights = []

    for i in arange(len(rho)):
        #a = norm.pdf(mu[i], loc=mu[i], scale=sig[i]) #height of Gaussian
        #heights.append(rho[i] * a * N)
        a = gy[int(mu[i])]
        heights.append(a * N)
        
        #plot(x, rho[i]*norm.pdf(x, mu[i], sig[i]), ':')
        plot([x[0],x[-1]],[a, a])

    return heights


def main(covfile, order, outfile, outplot, fraglen, width, rmsd_co):

    f = open(covfile)
    ID = f.readline().split()[3]
    f.close()

    a = loadtxt(covfile, usecols=[4,5])

    xs = a.T[0]
    counts = a.T[1]

    #reduce histogram (coverage) if total coverage is too high (EM gets slow)
    if sum(counts) >= 100000:
        ys = counts/sum(counts)
        ys *= 100000
    else:
        ys = counts

    minx = min(xs)
    maxx = max(xs)

    ##produce data counts (like histogram)
    N = len(xs) #data range 1-1000
    x = []
    for i in arange(N):
        for j in arange(ys[i]):
            x.append(xs[i])
    x = array(x)

    ##start EM
    ##produce initial values for EM
    start = (maxx-minx)/(order+2)

    mu = linspace(minx+start-1, maxx-start+1, order) #use some equally spaced initial mu.
    sig = []
    rho = []

    for i in arange(order):
        sig.append(randint.rvs(10,100))
        rho.append((1.0/order)-0.01)


    for ijk in [1,2,3,4,5]: #run MMEM at most 5 times. If then cut-off still isn't made, then it's a bad peak...

        rho, mu, sig = EM2(x, rho, mu, sig, fraglen)

        ##build mixture model
        y = (1-sum(rho)) * uniform.pdf(xs , loc=min(xs), scale = max(xs))
        for i in arange(len(mu)):
            y += rho[i] * norm.pdf(xs, mu[i], sig[i])


        #Collapsed statistics. They are just used for the range of RMSD assessment and height computation.
        #Mixture is still built with EM sigma, mu and rho
        Crho, Cmu, Csig = collapseOverlap(mu, sig, rho, width)

        RMSD = getRMSD(y, counts/sum(counts), Cmu, Csig)
        if len(where(array(RMSD) > rmsd_co)[0]) != 0:
            print "\n%s has not made RMSD cut-off in run %s\nRMSDs: %s" %(name, ijk, RMSD)
            continue
        else:
            break

        #else do this:

    #produce plots of single Gaussians
    for i in arange(len(mu)):
        plot(xs, rho[i] * norm.pdf(xs, mu[i], sig[i]), '--')


    #subtract uniform form mixture profile again, because this is the background I do not want to use for height (I think. or do heights then become uncomparable between regions?)
    height = getHeight( y-(1-sum(rho)) * uniform.pdf(xs , loc=min(xs), scale = max(xs)) , Cmu, Csig, Crho, counts)

    
    plot(xs, counts/sum(counts), 'g', linewidth=2)
    plot(xs,y, 'k', linewidth=1.5)


    o = open(outfile, 'a')

    n = name.split('_')
    chrom = '_'.join(n[:-2])
    start = int(n[-2])

    for i in arange(len(Cmu)):
        l = Cmu[i]-(width*Csig[i])
        r = Cmu[i]+(width*Csig[i])
        plot([l,l],[0,max(y)], label = '%i %i %.1f %s' %(l, r, height[i], '{:.3e}'.format(float(RMSD[i])) ))
        plot([r,r],[0,max(y)])

        o.write('%s\t%i\t%i\t%s.p%s\t%.3f\t%s\t%.3f\t%.3f\n' %(chrom, int(start+l), int(start+r), ID, str(i+1), height[i], RMSD[i], Cmu[i], Csig[i]))

    legend(bbox_to_anchor=(1.1, 1.1))
    savefig(outplot)
    close()
    o.close()


if __name__ == '__main__':

    if len(sys.argv) != 8:
        print '\nUsage: python prog.py filesfile outfile plotdir fraglen order width RMSD_cutoff\n'
        sys.exit(0)

    filesfileroot = sys.argv[1] #a file containing paths to files. filesfileroot.SGE_TASK_ID is actual file
    outfileroot = sys.argv[2] #files that contains final windows with scores and statistics. outfileroot.SGE_TASK_ID is used.
    plotdir = sys.argv[3] #where to store plots to. Folder is not created, it's just written to
    fraglen = float(sys.argv[4]) #fragment length to constrain sigmas
    order = int(sys.argv[5])
    width = float(sys.argv[6])
    rmsd_co = float(sys.argv[7]) #this is not used to sort out peaks. it is just used to rerun MMEM fitting when RMSD values are below this one. Peaks that do not make the cut-off get sorted out by SelectPeaksMM.py

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
        main(ppath, order, outfile, outplot, fraglen, width, rmsd_co)
