# coding: utf-8
import numpy
import scipy
from scipy.linalg import *

h4 = hilbert(4)
eigs_h4 = eig(h4)
print eigs_h4
