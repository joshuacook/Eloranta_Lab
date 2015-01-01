
# Problem
Eigenproblems abound in the modeling of nature.

Density Functional Theory models yield single particle-like Schrödinger
equations with a nonlinear potential term that accounts for all the many-body
interactions.

Electronic Structure Theory generates Schrodinger like equations.

# Goal

We seek a method for an efficient computational solution to the time-independent
Schrödinger equation.

## the Time-Independent Schrödinger Equation
Given, $\hat{H}$ is the Hamiltonian operator (or Hamiltonian matrix) for a
discrete system), $\psi$ is the wavefunction, and $E$, discrete energy states,
where

\begin{align*}
\hat{H}=T+V(x)=-\frac{\hbar^2}{2m}\Delta +V(x)\tag{1}
\end{align*}

for $\Delta$, the Laplacian.

We have

\begin{align*}
\hat{H}\psi_j(r)=E_j\psi_j(r),\ \ \ j\in\mathbb{N}\tag{2}
\end{align*}

For numerical calculations, we will typically take $H$ to be a matrix of
discrete values describing the position and momentum of particles in the system.

Note that this is an Eigenproblem of the form

$$A\mathbf{x}=\lambda\mathbf{x}$$

### The Imaginary Time Propagation Method
The imaginary time propagation method (ITP) relies on solving the time-dependent
Schrödinger equation:

\begin{align*}
i\hbar\frac{\partial \psi(r,t)}{\partial t}= \hat{H}\psi(r,t)\tag{3}
\end{align*}

in imaginary time where .

### Wick Rotation
We perform a Wick Rotation (setting $t=-i\tau$) to transform Eqn. 3 into a
simple diffusion equation


\begin{align*}
\frac{\partial \psi(r,\tau)}{\partial
\tau}=-\frac{\hat{H}}{\hbar}\psi(r,\tau)\tag{4}
\end{align*}

## Solution to the Diffusion Equation
The formal solution to eqn. 4 is given by

\begin{align*}
\psi(r,\tau)=\exp(-\hat{H}\tau/\hbar)\psi(r,0)\tag{5}
\end{align*}

We expand the initial state $\psi(r,0)$ in terms of the eigenfunctions
$\phi_j(r)$ the correspond to the eigenvalues $E_j$ for

\begin{align*}
\hat{H}\phi_j(r)=E_j\phi_j(r)\tag{6}
\end{align*}

The time evolution starting from the initial state $\psi(r,0)$ can now be
written as

\begin{align*}
\psi(r,\tau)=e^{-\hat{H}\tau/\hbar}\psi(r,0)=e^{-\hat{H}\tau/\hbar}\sum_{j=0}^\i
nfty a_j\psi_j(r)=\sum_{j=0}^\infty a_je^{E_j\tau/\hbar}\phi_j(r)\tag{7}
\end{align*}

## Imaginary Time Propagation as Iterative Solution
As $\tau\to\infty$, $\psi(r,\tau)$ becomes proportional to $\phi_0(r)$. In other
words, iterated $\psi$ functions will converge on the eigenfunction for the base
state of the time-**independent** equation (eq. 6). Here we are solving an
eigenproblem through an interative approximation of a differential equation.

*iterative differential equation methods*

*Olver finite differences*

## Finding the Ground State


	import numpy as np
	import numpy.linalg as la
	import math, time
	import matplotlib.pyplot as plt
	from sys import argv
	import datetime
	%matplotlib inline


	k = 100
   
	eps = 10E-6
	times = np.array([[0.,0.]]) 
	temp_times = times
	H = np.random.rand(k+200,k+200)
	H = H.T.dot(H)
	file = datetime.datetime.now().strftime("%Y%m%d%H%M%S")


	for i in range(k):
			# print i
			i = i+2
			n = i
			err = 1
			conv = 1
			int_H = H[0:n,0:n]
   
			start = time.clock()
   
			phi0 = np.random.rand(n)
			# print la.eig(H)[1].T
			CayleyN = (np.identity(n)-0.5*int_H)
			CayleyP = (np.identity(n)+0.5*int_H)
   
			while(conv > eps/100):
					phi1 = la.solve(CayleyP,CayleyN.dot(phi0))
					mu = math.sqrt(phi1.dot(phi1))
					phi1 = phi1/mu  
					conv = math.sqrt((np.abs(phi1)-np.abs(phi0)).dot(np.abs(phi1)-np.abs(phi0)))
					# err = math.sqrt(2)*math.sqrt(abs(phi1.dot(int_H.dot(int_H)).dot(phi1)- (phi1.dot(int_H.dot(int_H)).dot(phi1))**2))
					# print err
					phi0 = phi1
   
			end = time.clock()
			delta_t = end-start
			temp_times[0][0] = i
			temp_times[0][1] = delta_t
			times = np.concatenate((times,temp_times),axis=0)
			np.savetxt(file,times,fmt='%.4e')
   
		   
	plt.plot(times[:k,0],times[:k,1])
   
	plt.show()



![png](ExcitedStates_files/ExcitedStates_9_0.png)



	np.set_printoptions(precision=6)
	phi0




	array([ 0.098598,  0.100265,  0.095025,  0.095775,  0.09464 ,  0.096599,
			0.096508,  0.100957,  0.103805,  0.10603 ,  0.102552,  0.101034,
			0.10636 ,  0.104147,  0.096384,  0.101009,  0.098411,  0.101118,
			0.096546,  0.100696,  0.09975 ,  0.096243,  0.097646,  0.095947,
			0.106207,  0.101702,  0.100812,  0.095776,  0.093859,  0.097235,
			0.097111,  0.097374,  0.102151,  0.103593,  0.094456,  0.099608,
			0.100065,  0.102219,  0.094074,  0.095701,  0.098194,  0.107129,
			0.097895,  0.097316,  0.099701,  0.098345,  0.100888,  0.100346,
			0.098301,  0.09011 ,  0.11002 ,  0.097654,  0.09853 ,  0.097568,
			0.099519,  0.098629,  0.100716,  0.099087,  0.099384,  0.100641,
			0.100856,  0.105444,  0.101622,  0.101194,  0.095854,  0.100674,
			0.09899 ,  0.101403,  0.098962,  0.103598,  0.099548,  0.097491,
			0.102601,  0.104041,  0.104411,  0.102041,  0.097384,  0.096235,
			0.101188,  0.100598,  0.097123,  0.100782,  0.101379,  0.095636,
			0.102954,  0.098423,  0.095492,  0.094661,  0.097597,  0.093326,
			0.095669,  0.098604,  0.10304 ,  0.091931,  0.096429,  0.101558,
			0.096972,  0.093569,  0.103769,  0.102539,  0.109665])




	la.eig(int_H)[1][:,0]




	array([ 0.098598,  0.100265,  0.095025,  0.095775,  0.09464 ,  0.096599,
			0.096508,  0.100957,  0.103805,  0.106031,  0.102552,  0.101034,
			0.10636 ,  0.104147,  0.096384,  0.101009,  0.098411,  0.101118,
			0.096546,  0.100696,  0.09975 ,  0.096243,  0.097646,  0.095947,
			0.106207,  0.101702,  0.100812,  0.095777,  0.093859,  0.097235,
			0.097111,  0.097374,  0.102151,  0.103593,  0.094456,  0.099608,
			0.100065,  0.102219,  0.094074,  0.095701,  0.098194,  0.107129,
			0.097895,  0.097316,  0.099701,  0.098345,  0.100888,  0.100346,
			0.098301,  0.090111,  0.110021,  0.097654,  0.09853 ,  0.097568,
			0.099519,  0.09863 ,  0.100716,  0.099087,  0.099384,  0.100641,
			0.100856,  0.105444,  0.101621,  0.101194,  0.095854,  0.100674,
			0.09899 ,  0.101403,  0.098962,  0.103598,  0.099548,  0.097491,
			0.102601,  0.104041,  0.104411,  0.102041,  0.097384,  0.096235,
			0.101188,  0.100597,  0.097123,  0.100782,  0.101379,  0.095635,
			0.102954,  0.098423,  0.095493,  0.094661,  0.097597,  0.093326,
			0.095669,  0.098604,  0.10304 ,  0.091931,  0.096429,  0.101558,
			0.096972,  0.093569,  0.103769,  0.102539,  0.109665])




	abs(la.eig(int_H)[1][:,0])-abs(phi0)





	array([ -2.006604e-07,   1.700561e-07,  -1.771677e-07,  -2.059107e-07,
			 3.448798e-08,  -2.433046e-08,   1.131446e-07,   7.423343e-08,
			 5.117748e-08,   3.036601e-07,  -4.612399e-08,   6.389505e-10,
			-1.013733e-07,   1.927032e-07,  -1.619328e-08,   3.747558e-08,
			-1.331242e-07,  -7.244949e-08,   2.192848e-07,   9.029285e-08,
			-1.586413e-07,   1.424134e-07,  -7.180850e-08,   7.261967e-08,
			 4.058222e-09,  -1.176393e-07,   5.143330e-08,   7.148184e-08,
			 2.926113e-07,  -2.592652e-08,  -6.353175e-08,   5.808863e-10,
			-1.635231e-07,  -6.656144e-08,  -3.047453e-08,   1.700310e-07,
			-5.070044e-08,   1.312531e-07,   3.052739e-07,   1.363959e-07,
			 2.177418e-08,  -1.305050e-08,   8.164661e-08,  -4.111662e-09,
			-5.990691e-08,   1.022985e-07,  -9.020645e-08,   2.425612e-07,
			 4.294973e-08,   1.356902e-07,   1.626833e-07,  -1.879414e-07,
			-1.057461e-07,  -1.282285e-07,  -5.011395e-08,   1.993277e-07,
			-1.039867e-07,   4.286127e-08,  -2.186629e-07,   2.953255e-08,
			-1.799635e-08,   5.392744e-08,  -3.411294e-07,   6.603637e-08,
			 5.239992e-08,  -1.596268e-08,  -2.510092e-08,   1.882452e-08,
			-9.418556e-08,   9.788429e-08,   1.082589e-07,  -7.564633e-08,
			-2.450232e-07,   2.311293e-08,  -1.507158e-07,   2.145516e-08,
			-7.755745e-08,   1.188076e-07,  -1.801996e-08,  -2.165753e-07,
			 3.356247e-11,   1.141551e-07,   1.787539e-07,  -1.455793e-07,
			 2.174523e-09,   2.676469e-08,   1.192650e-07,  -3.180931e-08,
			-2.537567e-07,   9.266495e-08,   1.448192e-08,  -1.460201e-07,
			-9.528265e-08,  -4.103974e-08,   1.956326e-07,  -1.420731e-07,
			 4.031461e-08,  -1.365290e-07,   2.868258e-08,  -6.500240e-08,
			-4.943056e-08])




   

### Finding an Excited Eigenstate
Given our basic
