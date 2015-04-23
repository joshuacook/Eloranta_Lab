from KTBC import KTBC
import numpy as np
import numpy.linalg as nla
import matplotlib.pyplot as plt

def fixedFixedPlot(differential_equation,force_vec,xi,xf,h,mesh_density,should_I_plot=True):

    # independent vector for differential equation
    u_ind_vec = np.linspace(xi,xf,mesh_density)
    
    # difference equation solution
    K = h**2*KTBC('K',h-1)                  # this line calls the KTBC function
    soln = nla.solve(K,force_vec)
    soln = np.insert(soln,0,0)
    soln = np.insert(np.zeros(1),0,soln)
    soln_ind_vec = np.linspace(xi,xf,h+1)

    # plots
    if (should_I_plot == True): 
        plt.plot(u_ind_vec,u(u_ind_vec)) # plot differential equation
        plt.plot(soln_ind_vec,soln) # plot difference equation
    
    return soln_ind_vec, soln, u_ind_vec, u(u_ind_vec)