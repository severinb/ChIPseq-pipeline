#!/usr/bin/env python

import os, sys
from string import *
from pylab import *
from scipy.stats import *

def EM2(x, rho, mu, sigma, fraglen, maxDiff=0.1):

    logLik = []

    if len(rho) != len(mu) or len(mu) != len(sigma):
        print 'rho, sigma or mu not right\n'
        sys.exit(0)

    mix = len(rho)

    #the linear relation of maximum fragment length and mode of the sigma distribution of high peaks was 0.3*fraglen+28. So I took a range of +- 30
    #If a negative fragment length is given, sigma is not constrained
    if fraglen < 0:
        minsig = 10
        maxsig = 110
    else:
        #fit: 0.41387037 * fraglen + 17.28960321 (R2:  0.113487283812), std residuals:  19.0427328655
        #+- 1.5 std residuals
        minsig = max(1, 0.414*fraglen - 11.5)
        maxsig = 0.414*fraglen + 45.5


    minx = min(x)
    maxx = max(x)

    print "Allow sigma range of %s-%s" %(minsig, maxsig)

    while True:

        # calculate likelihood of x
        pG = [norm.pdf(x, mu[i], sigma[i]) for i in arange(mix)] #list of arrays
        pE = uniform.pdf(x,loc=min(x), scale=max(x))

        logLik.append(sum(log( sum([pG[i]*rho[i] for i in arange(mix)], axis=0) + (1-sum(rho))*pE)))
                
        # check for convergence
        if len(logLik)>1 and (logLik[-1] - logLik[-2]) < maxDiff:
            break

        # calculate probabilities
        ppG = [rho[i] * pG[i] /( sum([pG[j]*rho[j] for j in arange(mix)], axis=0) + (1-sum(rho)) * pE ) for i in arange(mix)] #list of arrays of probabilities
        ppE = (1-sum(rho))*pE/( sum([pG[j]*rho[j] for j in arange(mix)], axis=0) + (1-sum(rho)) * pE )

        # update parameter
        #do not allow mu to become smaller than min(x) or larger than max(x). Otherwise naNs can arise
        candmu = [sum(x*ppG[i])/sum(ppG[i]) for i in arange(mix)]
        mu = [min( max( minx, candmu[i]), maxx) for i in arange(mix)]

        #constrain sigma by minsig and maxsig
        #sigma = [sqrt(sum(ppG[i]*(x-mu[i])**2)/sum(ppG[i])) for i in arange(mix)]
        sigma = [min( max(minsig, sqrt(sum(ppG[i]*(x-mu[i])**2)/sum(ppG[i]))) , maxsig) for i in arange(mix)]

        #add pseudo count to rho. If pG and ppG get accidentally 0 (because mu is somehow out of range and probabilities get 0), mu gets nan (divide by zero), rho gets 0 which results in ppG = 0 again. So mu stays nan which results in problems. It also rescues this mu, i.e. the Gaussian would get lost otherwise.
        pc = 0.001
        rho = array([sum(ppG[i])/len(x) for i in arange(mix)])
        rho += pc
        rho /= (1 + mix*pc)


    print logLik[-1]
    print mu
    print sigma
    print rho

    return (rho, mu, sigma, logLik[-1])


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

        if not cstart and not cstop: #for initialization
            cstart = win[1]
            cstop = win[2]
            overlapping_wins.append([win])

        #use <= (including edges). This should deal with identical windows (otherwise identical windows are not treated as overlapping)
        elif win[1] <= cstart and win[2] > cstart: #next window is overlapping start of first window (I think this shouldn't happen since sortwins list is sorted)
            overlapping_wins[i].append(win)
            cstart = win[1]
            if cstop < win[2]:
                cstop = win[2]

        elif win[1] > cstart and win[1] < cstop: #whether start of next window is inside last (current window)
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

    # use smoothed profile, using moving average.
    win = ones(5)/5.
    y_av = convolve(y, win, 'same')
    y = y_av

    minl = 0
    maxr = len(y) - 1

    #get RMSD for each sigma-mu pair
    qs = []
    for i in arange(len(Cmu)):
        try:
            l = int(Cmu[i] - Csig[i])
        except ValueError:
            print Cmu, Csig
            sys.exit(1)

        r = int(Cmu[i] + Csig[i])

        l = max(l, minl)
        r = min(r, maxr)

        rmsd = mean(abs(y[l:r]-gy[l:r]))
        if isnan(rmsd):
            rmsd = 10.0 #give a large number (not 0, to not get problems with logarithms). I could also give a very small number..

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
        a = gy[int(mu[i]) - 1] #mus go from 1 to region length, indices from 0 to region length -1 
        heights.append(a * N)

        # plot([x[0],x[-1]],[a, a])

    return heights


def main(covfile, order, outfile, outplot, fraglen, width):

    f = open(covfile)
    ID = f.readline().split()[3]
    f.close()

    a = loadtxt(covfile, usecols=[4,5])

    xs = a.T[0]
    counts = a.T[1]

    #reduce histogram (coverage) if total coverage is too high (EM gets slow)
    if sum(counts) >= 50000:
        ys = counts/sum(counts)
        ys *= 50000
    else:
        ys = counts

    minx = min(xs)
    maxx = max(xs)

    #theoretically there is at least enough space for region length/ (2*fragment length) peaks. Thus set order to this number. Take 4 times fragment length, as this gives still enough Gaussians
    order = max(2, min(int(maxx/ (4*fraglen)), 30)) #constrain Gaussian number. Otherwise it can get very slow...
    print 'chosen number of Gaussians: ', order

    ##produce data counts (like histogram)
    N = len(xs) #data range 1-1000
    x = []
    for i in arange(N):
        for j in arange(floor(ys[i])): #take floor, so that numbers between 0 and 1 get 0 and thus no weight
            x.append(xs[i])
    x = array(x)

    #constrain sigma
    if fraglen < 0:
        minsig = 10
        maxsig = 110
    else:
        minsig = max(1, 0.414*fraglen - 11.5)
        maxsig = 0.414*fraglen + 45.5

    ##start EM
    ##produce initial values for EM
    start = (maxx-minx)/(order+2)

    #variables for best parameter fits
    Brmsd = 100
    Bmu = None
    Bsig = None
    Brho = None
    Bll = None

    for ijk in [1,2,3,4,5]: #run MMEM 5 times. Keep track of the best fit... Best fit is defined as the one that has least RMSD over the whole region.

        #initialize all parameters to fit. I initialize mu with random parameters because linspace results in the fit most of the time. Since I'm running things multiply I can allow random initiation here
        mu = array([randint.rvs(minx+start-1, maxx-start+1) for i in arange(order)])  #linspace(minx+start-1, maxx-start+1, order) #use some equally spaced initial mu.
        sig = []
        rho = []


        for i in arange(order):
            sig.append(randint.rvs(minsig, maxsig))
            rho.append((1.0/order)-0.01)

        rho, mu, sig, ll = EM2(x, rho, mu, sig, fraglen)

        ##build mixture model
        y = (1-sum(rho)) * uniform.pdf(xs , loc=min(xs), scale = max(xs))
        for i in arange(len(mu)):
            y += rho[i] * norm.pdf(xs, mu[i], sig[i])


        rmsd = mean(abs(counts/sum(counts)-y))
        print 'RMSD: ', rmsd, ijk

        #if rmsd < Brmsd:
        if ll > Bll or Bll == None:
            Brmsd = rmsd
            Bmu = mu
            Bsig = sig
            Brho = rho
            Bijk = ijk
            Bll = ll

    print 'Best Run: ', Bijk
    print 'Best whole region RMSD: ', Brmsd
    print 'Best log-likelihood: ', Bll
    print 'Best mu, sigma and rho: ', Bmu, Bsig, Brho


    #Collapsed statistics. They are just used for the range of RMSD assessment and height computation.
    #Mixture is still built with EM sigma, mu and rho
    Crho, Cmu, Csig = collapseOverlap(Bmu, Bsig, Brho, width)


    #build mixture model and produce plots of single Gaussians
    y = (1-sum(Brho)) * uniform.pdf(xs , loc=min(xs), scale = max(xs))
    for i in arange(len(Bmu)):
        y += Brho[i] * norm.pdf(xs, Bmu[i], Bsig[i])
        plot(xs, Brho[i] * norm.pdf(xs, Bmu[i], Bsig[i]), '--', color='0.6')


    RMSD = getRMSD(y/max(y), (counts/sum(counts))/max(y), Cmu, Csig) #normalize curves to 1, so that peaky peaks do not get penalized by larger deviations.

    #subtract uniform form mixture profile again, because this is the background I do not want to use for height (I think. or do heights then become uncomparable between regions?)
    height = getHeight( y-(1-sum(Brho)) * uniform.pdf(xs , loc=min(xs), scale = max(xs)) , Cmu, Csig, Crho, counts)
    uniform_height = (1-sum(Brho)) * uniform.pdf(int(mean(xs)) , loc=min(xs), scale = max(xs)) * sum(counts)  #it's not important that I look for the uniform height at the mean...

    #now normalize RMSD by height. I don't know why I did the following line. I think it doesn't make sense at all
    #RMSD = array(RMSD)/(array(height)+0.01)  #add small pseudo count to not divide by zero

    plot(xs, counts/sum(counts), 'r', linewidth=2)
    plot(xs,y, 'k', linewidth=1.5)


    o = open(outfile, 'a')

    n = name.split('_')
    chrom = '_'.join(n[:-2])
    start = int(n[-2])

    colourlist = ['b', 'g', 'r', 'c', 'm', 'y'] #used for plotting
    for i in arange(len(Cmu)):
        l = Cmu[i]-(width*Csig[i])
        r = Cmu[i]+(width*Csig[i])

        try:
            colr = colourlist[i]
        except IndexError:
            colr = colourlist[randint.rvs(0,6)]

        #plot([l,l],[0,max(y)], label = '%i %i %.1f %s' %(l, r, height[i], '{:.3e}'.format(float(RMSD[i])) ))

        # plot([l,l],[0,max(y)], colr, label = '%i %i %.1f' %(l, r, height[i]))
        # plot([r,r],[0,max(y)])

        bar(l,max(y), r-l, color=colr, linewidth=0, alpha=0.2, label="%i-%i %.1f" %(l, r, height[i]))

        o.write('%s\t%i\t%i\t%s.p%s\t%.3f\t%s\t%.3f\t%.3f\t%.3f\n' %(chrom, int(start+l), int(start+r), ID, str(i+1), height[i], RMSD[i], Cmu[i], Csig[i], uniform_height))

    legend(bbox_to_anchor=(1.1, 1.1))
    savefig(outplot)
    close()
    o.close()


if __name__ == '__main__':

    if len(sys.argv) != 7:
        print '\nUsage: python prog.py filesfile outfile plotdir fraglen order width\n'
        sys.exit(0)

    filesfileroot = sys.argv[1] #a file containing paths to files. filesfileroot.SGE_TASK_ID is actual file
    outfileroot = sys.argv[2] #files that contains final windows with scores and statistics. outfileroot.SGE_TASK_ID is used.
    plotdir = sys.argv[3] #where to store plots to. Folder is not created, it's just written to
    fraglen = float(sys.argv[4]) #fragment length to constrain sigmas
    order = int(sys.argv[5])
    width = float(sys.argv[6])

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
        main(ppath, order, outfile, outplot, fraglen, width)
