# Imaginary Time Progagation Eigensolver
Numerical Approach {#numerical-approach .unnumbered}
==================

Imaginary Time Propagation {#imaginary-time-propagation .unnumbered}
--------------------------

Matrix eigenvalue problems appear in many different areas of natural
sciences and it is therefore important to improve the computational
efficiency of the current algorithms for solving such problems
numerically.$^1$ The imaginary time propagation method (ITP) relies on
solving the time-dependent Schrödinger equation in imaginary time. We
perform a Wick rotation (setting $t=-i\tau$) to transform the
time-dependent Schrödinger into a diffusion equation:$^2$

$$\label{eq:itpSch}
\frac{\partial \psi(r,\tau)}{\partial
\tau}=-\frac{H}{\hbar}\psi(r,\tau) \implies \psi(r,\tau)=\exp(-H\tau/\hbar)\psi(r,0)$$

where $H$ is the hamiltonian matrix (operator) and $\psi$ is an
approximation to the lowest eigenvector (eigenfunction) at iteration
$\tau$.\

As $\tau\to\infty$, $\psi(r,\tau)$ will converge on the eigenvector
corresponding to the lowest eigenvalue.$^2$ At each time step of the
algorithm, the exponential operator is approximated using the Cayley
unitary form, transforming the eigenvalue problem into a linear problem:

$$\label{eq:CayleyExpansion}
\exp(-H\Delta\tau)\approx\biggP{1+\frac{1}{2}H\Delta\tau}^{-1}\biggP{1-\frac{1}{2}H\Delta\tau} \implies \biggP{1+\frac{1}{2}H\Delta\tau}\psi(r,\tau+\Delta\tau)=\biggP{1-\frac{1}{2}H\Delta\tau}\psi(r,\tau)$$

To start the algorithm, a random vector is typically chosen as the
initial state. The method can also be extended for excited states by
imposing a penalty term in $H$ such that orthogonality is retained. This
method bears some similarities to the well-known power and subspace
iteration methods.$^1$\

An upper limit for the absolute error $\Delta E$, present in the lowest
eigenvalue $E_i(\tau)$ at iteration $\tau$, can be written in terms of
the quantum mechanical standard deviation of $H$ and is used as the
stopping criterion:$^2$

$$\label{eq:error}
\Delta E_i = |E_i -\langle\psi_i(r,\tau)|H|\psi_i(r,\tau)|\leq\sqrt{2}\sqrt{\langle\psi_i(r,\tau)|H^2|\psi_i(r,\tau)\rangle-\langle\psi_i(r,\tau)|H|\psi_i(r,\tau)\rangle^2}$$
