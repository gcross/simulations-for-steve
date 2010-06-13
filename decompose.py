#@+leo-ver=4-thin
#@+node:gcross.20100612231244.1369:@thin decompose.py
from __future__ import division
from numpy import *
from numpy.linalg import *

set_printoptions(linewidth=132,precision=15)

import sys

try:
    lam = float(sys.argv[1])
    number_of_sites = int(sys.argv[2])
except:
    print "Usage: %s <lambda> <number of sites> [number of levels]" % sys.argv[0]
    sys.exit()

#@+others
#@+node:gcross.20100612231244.1370:Utility functions
#@+node:gcross.20100612231244.1371:otimes
def otimes(A,B):
    return multiply.outer(A,B).transpose(0,2,1,3).reshape(A.shape[0]*B.shape[0],A.shape[1]*B.shape[1])
#@-node:gcross.20100612231244.1371:otimes
#@+node:gcross.20100612231244.1372:vtimes
def vtimes(A,B):
    return multiply.outer(A,B).ravel()
#@-node:gcross.20100612231244.1372:vtimes
#@+node:gcross.20100612231244.1373:projector
def projector(s):
    return multiply.outer(s,s)

#@-node:gcross.20100612231244.1373:projector
#@+node:gcross.20100612231244.1374:operators
def operators(*ops):
    return reduce(otimes,ops)
#@-node:gcross.20100612231244.1374:operators
#@+node:gcross.20100612231244.1433:svd_decompose
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
#@-node:gcross.20100612231244.1433:svd_decompose
#@+node:gcross.20100612231244.1435:print_matrices
def print_matrices(matrices):
    d = int(sqrt(matrices.shape[1]))
    matrices = matrices.reshape(matrices.shape[0],d,d)
    indent = "    ["
    for m in matrices:
        for i in xrange(d):
            print ("%s(" + "(%.15f):."*d + "()):.") % tuple([indent] + list(m[i]))
            indent = "     "
        print "     ()"
        indent = "    ,"
    print "    ]"
#@-node:gcross.20100612231244.1435:print_matrices
#@-node:gcross.20100612231244.1370:Utility functions
#@+node:gcross.20100612231244.1375:Single site objects
#@+node:gcross.20100612231244.1376:q, d = 5
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
#@-node:gcross.20100612231244.1376:q, d = 5
#@+node:gcross.20100612231244.1377:b, d = 3
_3_0 = array([1,0,0])
_3_1 = array([0,1,0])
_3_I = array([0,0,1])

_3_id = identity(3)
_3_id_proj = diag([1,1,0])
#@-node:gcross.20100612231244.1377:b, d = 3
#@+node:gcross.20100612231244.1378:q, d = 2
_2_0 = array([1,0])
_2_1 = array([0,1])

_2_id = identity(2)
#@-node:gcross.20100612231244.1378:q, d = 2
#@-node:gcross.20100612231244.1375:Single site objects
#@+node:gcross.20100612231244.1379:Two-site objects
#@+node:gcross.20100612231244.1380:bq, d = 3x2
_32_10 = vtimes(_3_0,_2_1)
_32_01 = vtimes(_3_1,_2_0)

s = (1/sqrt(2))*(_32_01-_32_10)
_32_id_proj = otimes(_3_id_proj,_2_id)
_32_Phi = _32_id_proj - projector(s)

print "="*72
print "3 x 2 PHI"
svd_decompose(_32_Phi,3,2)
#@-node:gcross.20100612231244.1380:bq, d = 3x2
#@+node:gcross.20100612231244.1381:bq, d = 3x5
_35_10 = vtimes(_3_0,_5_O1)
_35_01 = vtimes(_3_1,_5_O0)

s = (1/sqrt(2))*(_35_01-_35_10)
_35_id_proj = otimes(_3_id_proj,_5_id_proj)
_35_Phi = _35_id_proj - projector(s)

print "="*72
print "3 x 5 PHI"
svd_decompose(_35_Phi,3,5)
#@-node:gcross.20100612231244.1381:bq, d = 3x5
#@+node:gcross.20100612231244.1382:qb, d = 5x3
_53_I = otimes(projector(_5_I),identity(3)-projector(_3_I)) + otimes(identity(5)-projector(_5_I),projector(_3_I))

s_lam = 1/sqrt(1+lam**2)*((lam/sqrt(2))*(vtimes(_5_C1,_3_0)-vtimes(_5_C0,_3_1))+vtimes(_5_I,_3_I))
_53_Lambda = projector(s_lam) + _53_I

print "="*72
print "5 x 3 LAMBDA"
svd_decompose(_53_Lambda,5,3)
#@-node:gcross.20100612231244.1382:qb, d = 5x3
#@-node:gcross.20100612231244.1379:Two-site objects
#@-others

#@-node:gcross.20100612231244.1369:@thin decompose.py
#@-leo
