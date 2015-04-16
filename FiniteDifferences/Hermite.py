import numpy as np 
indep_vec = np.linspace(0,1,101)
# indep_vec is a discrete representation
# of the closed interval [0,1]

# Define all the variables
mu = 1.008
k = 1
h_ = 6.55E-34
alpha = np.sqrt((k*mu)/h_**2)

#  Build a Linear Function
H_0 = lambda x: 1
H_1 = lambda x: 2*x*np.sqrt(alpha)
H_2 = lambda x: 4*x*np.sqrt(alpha)- 2*x*np.sqrt(alpha)
H_3 = lambda x: 8*(x*np.sqrt(alpha))**3 - 12*x*np.sqrt(alpha)

# Discrete Functions
# H_0: [0,1] to REALS

herm_0 = H_0(indep_vec)

# Discrete Functions
# H_1: [0,1] to REALS
herm_1 = H_1(indep_vec)

# Discrete Functions
# H_2: [0,1] to REALS
herm_2 = H_2(indep_vec)

# Discrete Functions
# H_3: [0,1] to REALS
herm_3 = H_3(indep_vec)

