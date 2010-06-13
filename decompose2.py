#@+leo-ver=4-thin
#@+node:gcross.20100612231244.1339:@thin decompose2.py
from __future__ import division
from numpy import *
from numpy.linalg import *

set_printoptions(linewidth=132)

import sys

try:
    lam = float(sys.argv[1])
    number_of_sites = int(sys.argv[2])
except:
    print "Usage: %s <lambda> <number of sites> [number of levels]" % sys.argv[0]
    sys.exit()

#@+others
#@+node:gcross.20100612231244.1340:Utility functions
#@+node:gcross.20100612231244.1341:otimes
def otimes(A,B):
    return multiply.outer(A,B).transpose(0,2,1,3).reshape(A.shape[0]*B.shape[0],A.shape[1]*B.shape[1])
#@-node:gcross.20100612231244.1341:otimes
#@+node:gcross.20100612231244.1342:vtimes
def vtimes(A,B):
    return multiply.outer(A,B).ravel()
#@-node:gcross.20100612231244.1342:vtimes
#@+node:gcross.20100612231244.1343:projector
def projector(s):
    return multiply.outer(s,s)

#@-node:gcross.20100612231244.1343:projector
#@+node:gcross.20100612231244.1344:operators
def operators(*ops):
    return reduce(otimes,ops)
#@-node:gcross.20100612231244.1344:operators
#@+node:gcross.20100612231244.1384:svd_decompose
def svd_decompose(O,d1,d2):
    U,S,V = svd(O.reshape(d1,d2,d1,d2).transpose(0,2,1,3).reshape(d1*d1,d2*d2))
    d = 0
    while (d < len(S) and abs(S[d]) > 1e-10):
        d += 1
    print "U="
    print_matrices(U[:d]*S[:d].reshape(d,1))
    V = V.transpose()
    print "V="
    print_matrices(V[:d]*S[:d].reshape(d,1))
#@-node:gcross.20100612231244.1384:svd_decompose
#@+node:gcross.20100612231244.1385:print_matrices
def print_matrices(matrices):
    d = int(sqrt(matrices.shape[1]))
    matrices = matrices.reshape(matrices.shape[0],d,d)
    indent = "    ["
    for m in matrices:
        for i in xrange(d):
            print ("%s(" + "(%f):."*d + "()):.") % tuple([indent] + list(m[i]))
            indent = "     "
        print "     ()"
        indent = "    ,"
#@-node:gcross.20100612231244.1385:print_matrices
#@-node:gcross.20100612231244.1340:Utility functions
#@+node:gcross.20100612231244.1345:Single site objects
#@+node:gcross.20100612231244.1346:q, d = 5
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
#@-node:gcross.20100612231244.1346:q, d = 5
#@+node:gcross.20100612231244.1347:b, d = 3
_3_0 = array([1,0,0])
_3_1 = array([0,1,0])
_3_I = array([0,0,1])

_3_id = identity(3)
#@-node:gcross.20100612231244.1347:b, d = 3
#@+node:gcross.20100612231244.1348:q, d = 2
_2_0 = array([1,0])
_2_1 = array([0,1])

_2_id = identity(2)
#@-node:gcross.20100612231244.1348:q, d = 2
#@-node:gcross.20100612231244.1345:Single site objects
#@+node:gcross.20100612231244.1349:Two-site objects
#@+node:gcross.20100612231244.1350:bq, d = 3x2
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
svd_decompose(_32_H_Phi,3,2)
#@-node:gcross.20100612231244.1350:bq, d = 3x2
#@+node:gcross.20100612231244.1351:bq, d = 3x5
_35_USing = projector(vtimes(_3_0,_5_U1)-vtimes(_3_1,_5_U0))/2
_35_H_Phi = otimes(diag([1,1,0]),projector(_5_U0)+projector(_5_U1))-_35_USing
#@-node:gcross.20100612231244.1351:bq, d = 3x5
#@+node:gcross.20100612231244.1352:qb, d = 5x3
_53_H_idle = otimes(diag([1,1,1,1,0]),diag([0,0,1])) + otimes(diag([0,0,0,0,1]),diag([1,1,0]))

_53_USing = 1/sqrt(2) * (vtimes(_5_U1,_3_0) - vtimes(_5_U0,_3_1))
_53_Lambda = 1/sqrt(1+lam**2)*(lam*_53_USing + vtimes(_5_I,_3_I))

_53_H_L = projector(_53_Lambda) + _53_H_idle
#@-node:gcross.20100612231244.1352:qb, d = 5x3
#@-node:gcross.20100612231244.1349:Two-site objects
#@-others

#@-node:gcross.20100612231244.1339:@thin decompose2.py
#@-leo
