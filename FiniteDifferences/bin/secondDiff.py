import random
import numpy as np
import numpy.linalg as nla
import scipy.linalg as sla
import scipy.sparse as sprs
from sympy.abc import w,x,y,z

def secondDiff(type,n=10,sparse=False):
    '''
    secondDiff Create finite difference model matrix.
    D = secondDiff(TYPE,N,SPARSE) creates model matrix TYPE of size N-by-N.

    TYPE is one of the characters 'D', 'M', 'N', or 'C'.
    The command D = secondDiff('D', 100, 1) gives a sparse representation
    D = secondDiff uses the defaults TYPE='D', n=10, and SPARSE=False.
    Change the 3rd argument from 1 to 0 for dense representation!
    If no 3rd argument is given, the default is dense
    If no argument at all, secondDiff will give 10 by 10 matrix D
    '''

    e = np.ones(n)
    e_off = np.ones(n-1)
    D = sprs.diags([e_off,-2*e,e_off],[-1,0,1])

    if str(type) == 'D': 
        D = D
    if str(type) == 'M':
        D = sprs.csr_matrix(D)
        D[0,0] = -1
    if str(type) == 'R': 
        D = sprs.csr_matrix(D)
        D[0,0] = -1
        D[n-1,n-1] = -1
    if str(type) == 'C':
        D = sprs.csr_matrix(D)
        D[0,n-1] = 1
        D[n-1,0] = 1


    if sparse == False:
        return D.todense()
    else:
        return D
    
def pretty_vector(v):
    for i in range(v.size):
        print v[0,i]