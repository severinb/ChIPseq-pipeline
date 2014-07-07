################################
## Spectral clustering example, using python and scipy
##
## Coded by Nicolau Werneck <nwerneck@gmail.com> in 2011-03-10
## Based on "A Tutorial on Spectral Clustering", by Ulrike von Luxburg
################################

from pylab import *
#import scipy
from scipy.sparse.linalg import eigen
import numpy as np

## The similarity function
def dist(x,y):
    sig = 3.0
    return exp( -sig*(x-y)**2 )


## The "main" function... notice that running
## this will save three images to the local
## directory. If you want just ot see them you
## need to uncomment the ion(), and comment
## the savefig()s

if __name__=='__main__':

    #ion()

    print eigen

    Ns = 1000
    who = array(floor(rand(Ns)*5), dtype=np.int)

    mu = array([0, 6, 12, 12, 18])
    sig = array([1,1,1,1,1]) * 1.05

    q = randn(Ns)*sig[who]+mu[who]
   
    q.sort()

    figure(1, figsize=(6,3))
    subplot(1,2,1)
    title('Input data, sorted')
    plot(q)
    subplot(1,2,2)
    title('Hitogram of input data')
    hist(q, bins=sqrt(Ns), range=(-3,21)  )

    savefig('spec1.png', dpi=100)


    ## The graph matrix with similarity weights...
    W = zeros((Ns,Ns))

    Knn = 20

    ## Put the values in the matrix according to k-nearest neighbours
    ## approach. OBS: we have are considering the inputs have been
    ## sorted...
    for n in range(Ns):
        for j in range( max(n-Knn, 0), min(Ns, n+Knn) ):
            W[n,j] = dist(q[n], q[j])
            W[j,n] = W[n,j]


    ## The degree matrix

    D = diag(np.sum(W,0))

    L = identity(Ns) - dot(inv(D), W).T


    vw=np.max(abs(W))
    vl=np.max(abs(L))


    figure(2, figsize=(6,3))
    suptitle('Weight and Laplacian matrices')
    subplot(1,2,1)
    imshow(W, interpolation='nearest', cmap='RdBu', vmin=-vw, vmax=vw)
    subplot(1,2,2)
    imshow(L, interpolation='nearest', cmap='RdBu', vmin=-vl/20, vmax=vl/20)
    savefig('spec2.png', dpi=100)


    figure(3, figsize=(6,3))
    lam,u = eigen(L , k=4, which='SR')
    lamall,ulixo = eigen(L , k=10, which='SR')

    lamall=real(lamall)

    lam = real(lam)
    u = real(u)

    print lam
    subplot(1,2,1)
    title('First eigenvalues')
    plot(sort(lamall), '-+')
    subplot(1,2,2)
    title('First eigenvectors')
    plot(u[:,:4])
    savefig('spec3.png', dpi=100)
