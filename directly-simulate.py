#@+leo-ver=4-thin
#@+node:gcross.20100528235053.1304:@thin directly-simulate.py
from __future__ import division
from numpy import *
from numpy.linalg import *

import sys

try:
    lam = float(sys.argv[1])
    n = int(sys.argv[2])
except ValueError:
    print "Usage: %s [lambda] [number of middle 5-3 clusters]" % sys.argv[0]
    sys.exit()

if n > 2:
    print "Program override:  %i is too big of a value for n, idiot." % n
    sys.exit()

#@+others
#@+node:gcross.20100528235053.1313:Utility functions
#@+node:gcross.20100528235053.1314:otimes
def otimes(A,B):
    return multiply.outer(A,B).transpose(0,2,1,3).reshape(A.shape[0]*B.shape[0],A.shape[1]*B.shape[1])
#@-node:gcross.20100528235053.1314:otimes
#@+node:gcross.20100528235053.1315:vtimes
def vtimes(A,B):
    return multiply.outer(A,B).ravel()
#@-node:gcross.20100528235053.1315:vtimes
#@+node:gcross.20100528235053.1316:projector
def projector(s):
    return multiply.outer(s,s)

#@-node:gcross.20100528235053.1316:projector
#@+node:gcross.20100528235053.1319:operators
def operators(*ops):
    return reduce(otimes,ops)
#@-node:gcross.20100528235053.1319:operators
#@-node:gcross.20100528235053.1313:Utility functions
#@+node:gcross.20100528235053.1305:Single site objects
#@+node:gcross.20100528235053.1306:q, d = 5
_5_O0 = array([1,0,0,0,0])
_5_O1 = array([0,1,0,0,0])
_5_C0 = array([0,0,1,0,0])
_5_C1 = array([0,0,0,1,0])
_5_I  = array([0,0,0,0,1])

_5_id = identity(5)
_5_id_proj = diag([1,1,0,0,0])
_5_U = (1/2) * array(
    [[ 1, 0,-1, 0,0],
     [ 0, 1, 0,-1,0],
     [-1, 0, 1, 0,0],
     [ 0,-1, 0, 1,0],
     [ 0, 0, 0, 0,0]]
    )
_5_in = diag([0,1,0,0,0])
#@-node:gcross.20100528235053.1306:q, d = 5
#@+node:gcross.20100528235053.1307:b, d = 3
_3_0 = array([1,0,0])
_3_1 = array([0,1,0])
_3_I = array([0,0,1])

_3_id = identity(3)
_3_id_proj = diag([1,1,0])
#@-node:gcross.20100528235053.1307:b, d = 3
#@+node:gcross.20100528235053.1308:q, d = 2
_2_0 = array([1,0])
_2_1 = array([0,1])

_2_id = identity(2)
#@-node:gcross.20100528235053.1308:q, d = 2
#@-node:gcross.20100528235053.1305:Single site objects
#@+node:gcross.20100528235053.1309:Two-site objects
#@+node:gcross.20100528235053.1310:bq, d = 3x2
_32_10 = vtimes(_3_0,_2_1)
_32_01 = vtimes(_3_1,_2_0)

s = (1/sqrt(2))*(_32_01-_32_10)
_32_id_proj = otimes(_3_id_proj,_2_id)
_32_Phi = _32_id_proj - projector(s)
#@-node:gcross.20100528235053.1310:bq, d = 3x2
#@+node:gcross.20100528235053.1311:bq, d = 3x5
_35_10 = vtimes(_3_0,_5_O1)
_35_01 = vtimes(_3_1,_5_O0)

s = (1/sqrt(2))*(_35_01-_35_10)
_35_id_proj = otimes(_3_id_proj,_5_id_proj)
_35_Phi = _35_id_proj - projector(s)
#@-node:gcross.20100528235053.1311:bq, d = 3x5
#@+node:gcross.20100528235053.1312:qb, d = 5x3
_53_I = otimes(projector(_5_I),identity(3)-projector(_3_I)) + otimes(identity(5)-projector(_5_I),projector(_3_I))

s_lam = 1/sqrt(1+lam**2)*((lam/sqrt(2))*(vtimes(_5_C1,_3_0)-vtimes(_5_C0,_3_1))+vtimes(_5_I,_3_I))
_53_Lambda = projector(s_lam) + _53_I
#@-node:gcross.20100528235053.1312:qb, d = 5x3
#@-node:gcross.20100528235053.1309:Two-site objects
#@+node:gcross.20100528235053.1317:Hamiltonian
H_in = operators(_5_in,identity(((3*5)**n)*3*2))
H_U = sum(operators(
        identity((5*3)**i),
        _5_U,
        identity((3*5)**(n-i)*3*2)
    )   for i in xrange(n+1)
    )
H_Phi = sum(operators(
        identity(5*(3*5)**i),
        _35_Phi,
        identity((3*5)**(n-i-1)*3*2)
    ) for i in xrange(n)
    ) \
    + operators(identity(5*(3*5)**n),_32_Phi)
H_Lambda = sum(operators(
        identity((5*3)**i),
        _53_Lambda,
        identity((5*3)**(n-i)*2)
    )   for i in xrange(n+1)
    )

H = H_in+H_U+H_Phi+H_Lambda
#@-node:gcross.20100528235053.1317:Hamiltonian
#@-others

#evals, evecs = eigh(H)
from scipy.sparse.linalg.eigen.arpack.speigs import ARPACK_eigs
def matvec(x):
    return dot(H,x)
evals, evecs = ARPACK_eigs(matvec,H.shape[0],2)
for eval in evals[:2]:
    print eval
#@-node:gcross.20100528235053.1304:@thin directly-simulate.py
#@-leo
