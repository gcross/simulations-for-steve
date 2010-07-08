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
    U = U.transpose().reshape(d1*d1,d1,d1)
    V = V.reshape(d2*d2,d2,d2)
    d = 0
    while (d < len(S) and abs(S[d]) > 1e-10):
        d += 1
    U = U[:d]
    V = V[:d]
    S = S[:d]
    VS = V*S.reshape(d,1,1)
    print "U="
    print_matrices(U)
    print "V="
    print_matrices(VS)
    M = recombine(U,VS,d1*d2)
    assert allclose(M,O)
#@-node:gcross.20100612231244.1384:svd_decompose
#@+node:gcross.20100707155745.1405:recombine
def recombine(U,V,d):
    M = zeros((d,d),complex128)
    for u, v in zip(U,V):
        M += otimes(u,v)
    return M
#@-node:gcross.20100707155745.1405:recombine
#@+node:gcross.20100612231244.1385:print_matrices
def print_matrices(matrices):
    d = matrices.shape[1]
    indent = "    ["
    for m in matrices:
        for i in xrange(d):
            print ("%s(" + "(%f):."*d + "()):.") % tuple([indent] + list(m[i]))
            indent = "     "
        print "     ()"
        indent = "    ,"
    print "    ]"
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
#print "Decomposing 3 x 2 Phi..."
#svd_decompose(_32_H_Phi,3,2)
#@-node:gcross.20100612231244.1350:bq, d = 3x2
#@+node:gcross.20100612231244.1351:bq, d = 3x5
_35_USing = projector(vtimes(_3_0,_5_U1)-vtimes(_3_1,_5_U0))/2
_35_H_Phi = otimes(diag([1,1,0]),projector(_5_U0)+projector(_5_U1))-_35_USing
#print "Decomposing 3 x 5 Phi..."
#svd_decompose(_35_H_Phi,3,5)
#@-node:gcross.20100612231244.1351:bq, d = 3x5
#@+node:gcross.20100612231244.1352:qb, d = 5x3
_53_H_idle = otimes(diag([1,1,1,1,0]),diag([0,0,1])) + otimes(diag([0,0,0,0,1]),diag([1,1,0]))

_53_USing = 1/sqrt(2) * (vtimes(_5_U1,_3_0) - vtimes(_5_U0,_3_1))
_53_Lambda = 1/sqrt(1+lam**2)*(lam*_53_USing + vtimes(_5_I,_3_I))
_53_H_L = projector(_53_Lambda) + _53_H_idle

A1U = 1/sqrt(2) * _5_U1
B1U = -1/sqrt(2) * _5_U0
C0U = _5_I

A1V = _3_0
B1V = _3_1
C0V = _3_I

L0 = 1/(1+lam**2)
L1 = lam/(1+lam**2)
L2 = lam**2/(1+lam**2)

LU = array([
    diag([1,1,1,1,0]),
    diag([0,0,0,0,1]),
    outer(A1U,A1U),
    outer(A1U,B1U),
    outer(A1U,C0U),
    outer(B1U,A1U),
    outer(B1U,B1U),
    outer(B1U,C0U),
    outer(C0U,A1U),
    outer(C0U,B1U),
    outer(C0U,C0U),
])

LS = array([
    1,
    1,
    L2,
    L2,
    L1,
    L2,
    L2,
    L1,
    L1,
    L1,
    L0,
])

LV = array([
    diag([0,0,1]),
    diag([1,1,0]),
    outer(A1V,A1V),
    outer(A1V,B1V),
    outer(A1V,C0V),
    outer(B1V,A1V),
    outer(B1V,B1V),
    outer(B1V,C0V),
    outer(C0V,A1V),
    outer(C0V,B1V),
    outer(C0V,C0V),
])

print "Decomposing Lambda..."
print "U ="
print_matrices(LU)
print "V ="
print_matrices(LV)

LVS = LV*LS.reshape(11,1,1)

recombined_53_H_L = recombine(LU,LVS,15)
assert allclose(recombined_53_H_L,_53_H_L)
#@-node:gcross.20100612231244.1352:qb, d = 5x3
#@-node:gcross.20100612231244.1349:Two-site objects
#@-others

#@-node:gcross.20100612231244.1339:@thin decompose2.py
#@-leo
