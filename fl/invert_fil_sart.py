# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 20:11:45 2016

@author: James
"""

import numpy as np
import scipy.sparse as sps

def invert_sart(csc_A,csc_b,max_iterations=50,lam_start=1.0):
    # This code has been checked against Scott Silburn's Matlab code

    shap = csc_A.shape
    lam = lam_start
    colsum = (csc_A.transpose()).dot(sps.csc_matrix(np.ones(shap[0])).transpose())
    lamda = colsum
    #lamda = lamda.multiply(colsum != 0)
    np.reciprocal(lamda.data,out=lamda.data)
    np.multiply(lamda.data,lam,out=lamda.data)
    
    # Initialise output
    sol = sps.csc_matrix(np.zeros((shap[1],1))+np.exp(-1))
    # Create an array to monitor the convergence
    conv = np.zeros(max_iterations)
    
    for i in range(max_iterations):
        # Calculate sol_new = sol+lambda*(x'*(b-Ax))
        tmp = csc_b.transpose()-csc_A.dot(sol)
        tmp2 = csc_A.transpose().dot(tmp)
        #newsol = sol+tmp2*lamda
        newsol = sol+tmp2.multiply(lamda)
        # Eliminate negative values
        newsol = newsol.multiply(newsol > 0.0)
        newsol.eliminate_zeros()
        # Calculate how quickly the code is converging
        conv[i] = (sol.multiply(sol).sum()-newsol.multiply(newsol).sum())/sol.multiply(sol).sum()
        # Set the new solution to be the old solution and repeat
        sol = newsol
        
    return newsol.todense(), conv