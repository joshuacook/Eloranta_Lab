# Imaginary Time Progagation Eigensolver
Numerical Approach
==================

Imaginary Time Propagation
--------------------------

Matrix eigenvalue problems appear in many different areas of natural
sciences and it is therefore important to improve the computational
efficiency of the current algorithms for solving such problems
numerically.

The implicitly restarted Lanczos method is the go to method for solving eigenproblems. The ITP method has been shown to have better scaling with respect to matrix size as compared to the implicitly restarted Lanczos. 

The method
----------

1. take the Cayley Unitary form of the matrix
2. iteratively solve the equation phi1 = (I+0.5*H)^-1*(I-0.5*H)*phi0
3. Check for convergence using the quantum mechanical standard deviation as a stopping criteria
4. this will yield the ground state eigenvalue
