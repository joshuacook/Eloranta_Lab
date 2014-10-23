# Power Method Solutions

A first project for my work in the Lab of Dr. Jussi Eloranta at California State University at Northridge. 

Herein, I bring myself up to speed with regard to iterative power method solutions to eigenproblems. 


# Simple Power Method

1. Choose a random normalized starting vector, *y*
2. while some convergence criteria is not satisfied
    i. solve the linear equation A*y*=*x*
    ii. normalize *x*
    iii. set *y* to *x*
