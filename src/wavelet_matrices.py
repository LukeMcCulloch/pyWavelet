#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 29 15:02:35 2017

@author: TLM

nullspace of a matrix:
    http://scipy-cookbook.readthedocs.io/items/RankNullspace.html
  is a better answer to:
    https://stackoverflow.com/questions/5889142/python-numpy-scipy-finding-the-null-space-of-a-matrix
"""

import numpy as np
import scipy as sp
import scipy.misc
import sympy as sy #LCD???  ...rational approximations.. why?

#factorial = sp.math.factorial
factorial = sp.misc.factorial
rank = np.linalg.matrix_rank

rat = sy.Rational
frac = sy.fraction

svd = np.linalg.svd


#
#***********************************************************
#
def Factorial(m):
    try:
        if len(m) > 1:
            tisarray = 1
        else:
            tisarray = 0   
    except:
        tisarray = 0
        
    if tisarray == 1:
        r,c = np.shape(m)
        f = np.zeros((r,c))
        for i in range(0,r):
            for j in range(0,c):
                f[i,j] = np.prod( np.arange(2,m[i,j]+1) )
        return f
    else:
        if m<0.:
            return -1.
        else:
            return factorial(m)

#
#***********************************************************
#              

def matlab_arii(i,e):
    return np.arange(i,e+1)

#
#***********************************************************
#
def Knots(d,j):
    """
        x = Knots(d, j) returns a vector 
        of knot values for B-spline scaling
        functions of degree d, level j.
    """
    aa = matlab_arii(0.,2.**j-1)/(2.**j)
    x = np.asarray([0. for el in range(d-1)] + list(aa) + [1. for el in range(d)])
    return x
    
#
#***********************************************************
#

def Greville(d,u):
    """verified
        x = Greville(d, u) 
        returns the vector of Greville abscissa values
        corresponding to degree d and knot vector u.
    """
    l = len(u)
    x = u[0:1-d]
    for k in range(2,d+1):
        x = x + u[k-1:l-d+k]
    return x / d

#
#***********************************************************
#
#def Choose(n,r):
#    return Factorial(n) / (Factorial(r) * Factorial(n-r))
def Choose(i,d):
    return np.divide( Factorial(i) , np.multiply( Factorial(d) , Factorial(i-d) ) )
#def hbj(p,j):
#    return factorial(p+1)/(factorial(j)*factorial(p+1-j))
#def Choose(p,j):
#    return factorial(p+1)/(factorial(j)*factorial(p+1-j))

def BernsteinInner(d):
    """
        I = BernsteinInner(d) returns the 
        matrix of inner products of Bernstein
        polynomials of degree d.
    """
    i = np.ones((d+1, 1),int)*np.arange(0,d+1,1,int)
    j = i.T
    I = np.divide( np.multiply( Choose(d, i) , Choose(d, j) ) , (Choose(2*d, i+j)*(2*d + 1)) )
    return I

def BernsteinWeights(d,j):
    w = np.identity(2**j + d)
    if d==0:
        return w
    u = Knots(d,j)
    g = Greville(d,u)
    
    for i in range(0,2**j-1):
        for r in range(0,d):
            u,g,w = InsertKnot(d,u,g,w,(i+1.)/(2.**j))
    return w

def Inner(d,j):
    I0 = BernsteinInner(d)
    n = 2**j + d
    I = np.zeros((n,n))
    w = BernsteinWeights(d,j)
    for k in range(0,n):
        #w1 = np.reshape(w[:,k],d+1,2**j)
        #w1 = w[:,k].reshape(d+1,2**j)
        w1 = w[:,k].reshape(2**j,-1)
        w1 = w1.reshape(2**j,d+1).T
        #setshape??
        for l in range(0,n):
            #w2 = w[:,l].reshape(d+1,2**j)
            
            w2 = w[:,l].reshape(2**j,-1)
            w2 = w2.reshape(2**j,d+1).T
            
            #I[k,l] = np.trace(w1.T*I0*w2)
            I[k,l] = np.matmul(np.matmul(w1.T,I0),w2).trace()
            I[l,k] = I[k,l]
    I = I / 2.**j
    return I

###----------------------------------------------------------------------------
### 2 Scale Relation for Basis Functions  TLMTLMTLMTLM below!
###----------------------------------------------------------------------------

"""
NOT USED!
#import scipy as sp
"""

def binomial(n_,i_):
    """
        P&T : i is scipy k
        (n,i) <=> (n,k) so i<=>k in the literature
        where
            n is `top`
            i is `bottom`
        
    """
    return sp.special.binom( n_,i_)
    
def TwoScale(i,k):
    return binomial(k,i)/(2.**(k-1))

def checkstein(i,n,u):
    return binomial(n,i)*(u**i)*((1.-u)**(n-i))

def Bernstein(i,n,u):
    """return the ith (n+1)th order, (n)th degree Bernstein polynomial
        only if nonzero at fixed u
        i : ith conrol point goes with ith Basis function (ith span index?)
        n : degree
        u : knot location
        B : Polynomial value
    Piegel and Tiller, page 20
    """
    K = n+1
    B      = np.zeros((K),float)
    B[n-i] = 1.0
    u1     = 1.0-u
    for k in range(1,K):
        for j in range(n,k-1,-1): #careful - sly index!
            B[j] = u1*B[j] + u*B[j-1]
    return B[n]

def AllBernstein(n,u):
    """return all of the ith (n+1)th degree Bernstein polynomial
        only compute if nonzero at fixed u
        n : degree
        u : knot location
        B : array of 
    Piegel and Tiller, page 21
    """
    K       = n+1
    B       = np.zeros((K),float)
    B[0]    = 1.0
    u1      = 1.0-u
    for j in range(1,K):
        saved = 0.
        for k in range(0,j):
            temp = B[k]
            B[k] = saved + u1*temp
            saved = u*temp
        B[j] = saved
    return B

#
#***********************************************************
#

def PolyEval(g, p, gnew):
    """
    % pret = PolyEval(g, p, gnew) returns the values of a control polygon
    % defined by abscissas g and ordinates p, evaluated at gnew.
    """
    m, n = np.shape(p)
    assert(np.size(g) == m),'PolyEval: Length of g and rows of p must be the same.'
    
    lgn = np.size(gnew)
    pret = np.zeros((lgn,n)) #TLM GUESS!
    #
    #*******
    #COMMON MISTAKE:
    #
    #for i in range(1,len(gnew)+1):
    #
    #****
    #Correction:
    #
    for i in range(0,len(gnew)):
      #row = max(find(g <= gnew(i)))
      row = max( (g <= gnew[i]).nonzero() )[-1]
      #row = (g <= gnew[i]).nonzero().max()
      #aaa = g <= gnew[i]
      #ara = aaa.nonzero()
      # when i=0 row = 1
      # corresponds to matlab
      # when i=1, row=1
      if row == m-1:
        pret[i,:] = p[m-1,:]
      else:
        frac = (g[row+1] - gnew[i])/(g[row+1] - g[row])
        pret[i,:] = frac*p[row,:] + (1 - frac)*p[row+1,:]
    return pret

#
#***********************************************************
#

def InsertKnot(d, u, g, p, unew):
    """
    % [uret, gret, pret] = InsertKnot(d, u, g, p, unew) inserts a new knot at
    % unew for B-spline scaling functions of degree d, thereby modifying knot
    % vector u, Greville abscissas g, and synthesis matrix p.
    """
    uret = np.sort(np.concatenate( (u , [unew]),axis=0))
    gret = Greville(d, uret)
    pret = PolyEval(g, p, gret)
    return uret, gret, pret
#
#***********************************************************
#
"""
NOT USED!
"""
#def matlab_cat1_old(Vec, Array):
#    Array[:,0] = Array[:,0]+Vec
#    return Array
#def matlab_cat2_old(Array, Vec):
#    Array[:,1] = Array[:,1]+Vec
#    return Array
#
#def matlab_cat1(Vec, Array):
#    return np.concatenate((Vec,Array), axis=1)
#def matlab_cat2(Array, Vec):
#    return np.concatenate((Array,Vec), axis=1)
#
#***********************************************************
#
def FindP(d,j):
    """
        returns the P matrix for B-spline scaling functions of 
        degree  : d 
        level   : j
    """
    d = int(np.fix(d))
    assert(d>=0),'Error, FindP:  Must have d >= 0.'
    assert(j >= 1),'Error, FindP:  Must have j >= 1.'
    if d == 0:
        #P = np.asarray([[1],[1]])
        P = np.asarray([[1.,1.]]).T
        
        #mm = np.zeros(shape = [2**j,2**(j-1)])
        for i in range(2,j+1):
            print i
            sp = np.shape(P)
            p1 = np.concatenate((P,np.zeros(sp)), axis=1)
            p2 = np.concatenate((np.zeros(sp),P), axis=1)
            P = np.array([list(p1) + list(p2)])[0]
            
            """
            P = np.array([ [ list(el) for el in p1]+
                            [list(el) for el in p2] ])
            #"""
    else:
        u = Knots(d,j-1)
        g = Greville(d,u)
        P = np.identity(2**(j-1) + d)
        for k in range(0, 2**(j-1)-1+1  ):
            u,g,P = InsertKnot(d, u, g, P, (2*k+1.)/2**j  )
    
    return P




#
#***********************************************************
#

#def null(A, eps=1e-15):
#    u, s, vh = sp.linalg.svd(A)
#    null_mask = (s <= eps)
#    null_space = sp.compress(null_mask, vh, axis=0)
#    return sp.transpose(null_space)

def nullspace(A, atol=1e-13, rtol=0):
    A = np.atleast_2d(A)
    u, s, vh = svd(A)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    return ns
null = nullspace

def gcd(*numbers):
    """Return the greatest common divisor of the given integers"""
    from fractions import gcd
    return reduce(gcd, numbers)

def lcm(*numbers):
    """Return lowest common multiple."""    
    def lcm(a, b):
        return (a * b) // gcd(a, b)
    return reduce(lcm, numbers, 1)

#
#***********************************************************
#
vfrac = np.vectorize(frac)
vrat = np.vectorize(rat)
def LCD(m):
    num, denom = vrat(m)
    return d
#LCD = lcm
#
#***********************************************************
#
normalization='L2'
def FindQ(d, j, normalization='L2'):
    P = FindP(d,j)
    I = Inner(d,j)
    #
    M = np.matmul(P.T,I)
    m1,m2 = np.shape(M)
    n = m2 - rank(M) #M.ndim #np.rank(M)
    Q = np.zeros((m2,n))
    found = 0
    start_col = 0
    while ( (found < n/2.) and (start_col < 2) ):
        #beware the matlab indices!  (used verbatum here)
        start_col =  start_col + 1 + int(found > d)
        width = 0
        rank_def = 0
        while(  not rank_def and (width < m2 - start_col +1) ):
            width = width + 1
            #submatrix = M[:,start_col:start_col+width-1]
            submatrix = M[:,start_col-1:start_col+width-1] #adjusted indices here!
            rank_def = width - rank(submatrix)
        if rank_def:
            print 'width = ', width
            q_col = null(submatrix)
            #--------------------------------------------------------
            if normalization == 'min':
                q_col = q_col/min(abs(q_col + 1e38*(abs(q_col) < 1e-10)))
            elif normalization == 'max':
                q_col = q_col/max(abs(q_col))
            elif normalization == 'lcd':
                print 'error LCD not implemented yet'
                pass
                q_col = q_col/min(abs(q_col + 1e38*(abs(q_col) < 1e-10)))
                q_col = q_col*LCD(q_col)
            #--------------------------------------------------------
            # change sign to give consistent orientation
            q_col = q_col*(-1)**(start_col + np.floor((d+1.)/2.) + (q_col[0,0] > 0))
            # correct any slight error for answers that should be integers
            #if np.all(abs(submatrix*np.round(q_col)) < 1e-10) and np.any(np.round(q_col) != 0):
            if np.all(abs(np.matmul(submatrix,np.round(q_col)) ) < 1e-10) and np.any(np.round(q_col) != 0):
                q_col = np.round(q_col)
        # put column into strcmpleft half of Q
        found = found + 1
        #Q[start_col:start_col+width-1+1,found] = q_col[:,0]
        Q[start_col-1:start_col-1+width-1+1,found-1] = q_col[:,0]
        #Q[start_col-1:start_col-1+width-1+1,found] = q_col[:,0]
        
        # use symmetry to put column into right half of Q in reverse order
        # and negated if degree is even
        Q[:,n-found] = np.flipud(Q[:,found-1])*(-1.)**(d+1.)
    if normalization=='L2':
        ip = np.matmul(Q.T,np.matmul(I,Q))
        Q = np.matmul(Q,np.diag(1./np.sqrt(np.diag(ip))))
    return Q

#
#***********************************************************
#

if __name__ == '__main__':
    u = np.asarray([1,1,1,0,0,0,1,1,1],float)
    print np.matrix(Greville(3,u)).T
    """
        matlab:
            ans =

                   1.00000
                   0.66667
                   0.33333
                   0.00000
                   0.33333
                   0.66667
                   1.00000
    """
    u = np.asarray([1,1,1,1,0,0,0,1,1,1,1],float)
    print np.matrix(Greville(4,u)).T
    """
        matlab:
            ans =

                   1.00000
                   0.75000
                   0.50000
                   0.25000
                   0.25000
                   0.50000
                   0.75000
                   1.00000
    """
    
    u = np.asarray([1.,2.,3.,4.,5.,6.])
    
    #degree:
    d = 2
    
    #level:
    j = 2
    
    t = Knots(d,j)
    """
        >> t = Knots(d,j);
        >> t
        t =
        
           0.00000   0.00000   0.25000   0.50000   0.75000   1.00000   1.00000
    #"""
    
    """
            >>
        >> Knots(d, j - 1)
        ans =
        
           0.00000   0.00000   0.50000   1.00000   1.00000
           >>> Knots(d, j - 1)
           array([ 0. ,  0. ,  0.5,  1. ,  1. ])
        
        >>
        >>
        >> Greville(d, u)
        ans =
        
           1.5000   2.5000   3.5000   4.5000   5.5000
           
          >>> Greville(d, u)
          array([ 1.5,  2.5,  3.5,  4.5,  5.5])

    #"""
    
#    u = Knots(d,j-1)
#    g = Greville(d,u)
#    P = np.identity(2**(j-1) + d)
#    k=0 #k=1
#    p = P
#    unew = (2*k+1.)/(2**j)
#    uret = np.sort(np.concatenate( (u , [unew]),axis=0))
#    gret = Greville(d, uret)
#    
#    gnew = gret
#    
#    m, n = np.shape(p)
#    
#    
#    u,g,P = InsertKnot(d, u, g, P, (2*k+1.)/(2**j)  )
#    
#    
    P1 = FindP(2,2)