#@+leo-ver=4-thin
#@+node:gcross.20100610125246.1321:@thin directly-simulate2.py
from __future__ import division
from numpy import *
from numpy.linalg import *

import sys

try:
    lam = float(sys.argv[1])
    number_of_sites = int(sys.argv[2])
except:
    print "Usage: %s <lambda> <number of sites> [number of levels]" % sys.argv[0]
    sys.exit()

try:
    number_of_levels = int(sys.argv[3])
except IndexError:
    number_of_levels = 2
except ValueError:
    print "Usage: %s <lambda> <number of sites> [number of levels]" % sys.argv[0]
    sys.exit()

if number_of_sites < 3 or number_of_sites % 2 == 0:
    print "The number of sites must be an odd integer >= 3."

n = (number_of_sites - 3) // 2

if n > 2:
    print "Program override:  %i is too big of a value for the number of sites, idiot." % number_of_sites
    sys.exit()

#@+others
#@+node:gcross.20100610125246.1322:Utility functions
#@+node:gcross.20100610125246.1323:otimes
def otimes(A,B):
    return multiply.outer(A,B).transpose(0,2,1,3).reshape(A.shape[0]*B.shape[0],A.shape[1]*B.shape[1])
#@-node:gcross.20100610125246.1323:otimes
#@+node:gcross.20100610125246.1324:vtimes
def vtimes(A,B):
    return multiply.outer(A,B).ravel()
#@-node:gcross.20100610125246.1324:vtimes
#@+node:gcross.20100610125246.1325:projector
def projector(s):
    return multiply.outer(s,s)

#@-node:gcross.20100610125246.1325:projector
#@+node:gcross.20100610125246.1326:operators
def operators(*ops):
    return reduce(otimes,ops)
#@-node:gcross.20100610125246.1326:operators
#@-node:gcross.20100610125246.1322:Utility functions
#@+node:gcross.20100610125246.1327:Single site objects
#@+node:gcross.20100610125246.1328:q, d = 5
_5_O0 = array([1,0,0,0,0])
_5_O1 = array([0,1,0,0,0])
_5_C0 = array([0,0,1,0,0])
_5_C1 = array([0,0,0,1,0])
_5_I  = array([0,0,0,0,1])

_5_id = identity(5)
_5_U0 = array([1,0,1,0,0])/sqrt(2)
_5_U1 = array([0,1,0,1,0])/sqrt(2)

_5_H_in = projector(_5_U1)
_5_H_U = (1/2) * array(
    [[ 1, 0,-1, 0,0],
     [ 0, 1, 0,-1,0],
     [-1, 0, 1, 0,0],
     [ 0,-1, 0, 1,0],
     [ 0, 0, 0, 0,0]]
    )
#@-node:gcross.20100610125246.1328:q, d = 5
#@+node:gcross.20100610125246.1329:b, d = 3
_3_0 = array([1,0,0])
_3_1 = array([0,1,0])
_3_I = array([0,0,1])

_3_id = identity(3)
#@-node:gcross.20100610125246.1329:b, d = 3
#@+node:gcross.20100610125246.1330:q, d = 2
_2_0 = array([1,0])
_2_1 = array([0,1])

_2_id = identity(2)
#@-node:gcross.20100610125246.1330:q, d = 2
#@-node:gcross.20100610125246.1327:Single site objects
#@+node:gcross.20100610125246.1331:Two-site objects
#@+node:gcross.20100610125246.1332:bq, d = 3x2
_32_H_Phi = zeros((6,6),float64)
terms = \
    [(1,1,1)
    ,(4,4,1)
    ,(2,2,1/2)
    ,(2,3,1/2)
    ,(3,2,1/2)
    ,(3,3,1/2)
    ]
for i, j, value in terms:
    _32_H_Phi[i-1,j-1] = value
#@-node:gcross.20100610125246.1332:bq, d = 3x2
#@+node:gcross.20100610125246.1333:bq, d = 3x5
_35_USing = projector(vtimes(_3_0,_5_U1)-vtimes(_3_1,_5_U0))/2
_35_H_Phi = otimes(diag([1,1,0]),projector(_5_U0)+projector(_5_U1))-_35_USing
#@-node:gcross.20100610125246.1333:bq, d = 3x5
#@+node:gcross.20100610125246.1334:qb, d = 5x3
_53_H_idle = otimes(diag([1,1,1,1,0]),diag([0,0,1])) + otimes(diag([0,0,0,0,1]),diag([1,1,0]))

_53_USing = 1/sqrt(2) * (vtimes(_5_U1,_3_0) - vtimes(_5_U0,_3_1))
_53_Lambda = 1/sqrt(1+lam**2)*(lam*_53_USing + vtimes(_5_I,_3_I))

_53_H_L = projector(_53_Lambda) + _53_H_idle
#@-node:gcross.20100610125246.1334:qb, d = 5x3
#@-node:gcross.20100610125246.1331:Two-site objects
#@+node:gcross.20100610125246.1335:Hamiltonian
H_in = operators(_5_H_in,identity(((3*5)**n)*3*2))
H_U = sum(operators(
        identity((5*3)**i),
        _5_H_U,
        identity((3*5)**(n-i)*3*2)
    )   for i in xrange(n+1)
    )
H_Phi = sum(operators(
        identity(5*(3*5)**i),
        _35_H_Phi,
        identity((3*5)**(n-i-1)*3*2)
    ) for i in xrange(n)
    ) \
    + operators(identity(5*(3*5)**n),_32_H_Phi)
H_Lambda = sum(operators(
        identity((5*3)**i),
        _53_H_L,
        identity((5*3)**(n-i)*2)
    )   for i in xrange(n+1)
    )

H = H_in+H_U+H_Phi+H_Lambda
#@-node:gcross.20100610125246.1335:Hamiltonian
#@-others

from scipy.sparse.linalg.eigen.arpack.speigs import ARPACK_eigs
def matvec(x):
    return dot(H,x)
#evals, evecs = eigh(H)
#evals = evals[:number_of_levels]
evals, evecs = ARPACK_eigs(matvec,H.shape[0],number_of_levels,which='SR')
for eval in sorted(evals):
    print eval
#@-node:gcross.20100610125246.1321:@thin directly-simulate2.py
#@-leo
