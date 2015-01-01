import numpy.linalg, numpy.random, numpy as np, math
from random import random as rand
from numpy.linalg import eig


B = numpy.array([[2,-12],[1,-5]])
y = numpy.array([1,1])
x = y
for i in range(100):
	y = B.dot(x)
	mu = math.sqrt(y.dot(y))
	x = y/mu
print numpy.transpose(eig(B)[1])
print x

[[ 0.9701425   0.24253563]
 [ 0.9486833   0.31622777]]
[ 0.9486833   0.31622777]



A = numpy.random.rand(3,3)
x = numpy.random.rand(3)
for i in range(20):
	y = A.dot(x)
	mu = math.sqrt(y.dot(y))
	x = y/mu
np.transpose(eig(A)[1]),x
