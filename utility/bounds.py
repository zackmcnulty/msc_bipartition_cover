'''
    In this file we collect functions for calculating our various bipartition cover bounds as
    functions of the relevant parameters: the species count k, minimmum branch length T_min, and the coverage probability q.
'''

import numpy as np
from utility.coalescent_probabilities import g_ij

#################################################################################
# Original Bound
#################################################################################


def original_bound(k, T_min, q=0.95):
    """
    Calculate the original bound:
    n >= log((1-q)/(k-3)) / log(1 - g_{k-2,1}(T_min))
    """
    if k <= 3:
        return float('inf')
    
    try:
        g_val = g_ij(k-2, 1, T_min)
        if g_val >= 1:  # Avoid log(0) or log(negative)
            return float('inf')
        
        numerator = np.log((1-q)/(k-3))
        denominator = np.log(1 - g_val)
        
        if denominator >= 0:  # log(1-g) should be negative
            return float('inf')
            
        return numerator / denominator
        
    except (ValueError, ZeroDivisionError):
        return float('inf')

#################################################################################
# Caterpillar Bound
#################################################################################

def caterpillar_bound(k, T_min, q=0.95):
    """
    Calculate the improved bound:
    n >= log((1-q)/(k-3)) / (1/(k-3) * sum_{l=2}^{k-2} log(1 - g_{l,1}(T_min)))
    """
    if k <= 3:
        return float('inf')
    
    try:
        # Calculate sum of log terms
        log_sum = 0
        for ell in range(2, k-1):  # l from 2 to k-2
            g_val = g_ij(ell, 1, T_min)
            if g_val >= 1:  # Avoid log(0) or log(negative)
                return float('inf')
            log_term = np.log(1 - g_val)
            if log_term >= 0:  # Should be negative
                return float('inf')
            log_sum += log_term
        
        if log_sum == 0:  # Avoid division by zero
            return float('inf')
            
        numerator = np.log((1-q)/(k-3))
        denominator = log_sum / (k-3)
        
        if denominator >= 0:  # Should be negative
            return float('inf')
            
        return numerator / denominator
        
    except (ValueError, ZeroDivisionError):
        return float('inf')

#################################################################################
# One-Step Bound
#################################################################################

def one_step_bound(k, T_min, q=0.95):
    """
    Calculate the improved bound:
    n >= log((1-q)/(k-3)) / (1/(k-3) * sum_{l=2}^{k-2} log(1 - q_ell))
    """
    if k <= 3:
        return float('inf')
    
    
    # Otherwise, pre-calculate all values of g_ij(T_min) we will need
    g_array = np.full((k, k), np.nan)
    for i in range(1, k):
        for j in range(1, i+1):
            g_array[i, j] = g_ij(i, j, T_min)
    
    try:
        # Calculate sum of log terms
        log_sum = 0
        for ell in range(2, k-1):  # l from 2 to k-2
            
            q_ell = 0
            
            # calculate q_ell
            m_ell = ell // 2
            for r in range(1, m_ell + 1):
                for s in range(1, ell - m_ell + 1):
                    q_ell += g_array[m_ell, r] * g_array[ell - m_ell, s] * g_array[r + s, 1]

            # Add log(1-q_ell) to running sum
            log_term = np.log(1 - q_ell)
            log_sum += log_term
        
        if log_sum == 0:  # Avoid division by zero
            return float('inf')
            
        numerator = np.log((1-q)/(k-3))
        denominator = log_sum / (k-3)
        
        if denominator >= 0:  # Should be negative
            return float('inf')
            
        return numerator / denominator
        
    except (ValueError, ZeroDivisionError):
        return float('inf')

#################################################################################
# Balanced Bound
#################################################################################


def balanced_bound(k, T_min, q=0.95):
    """
    Calculate the improved bound:
    n >= log((1-q)/(k-3)) / (1/(k-3) * sum_{l=2}^{k-2} log(1 - z_ell)) 
    
    as described above. 
    """
    if k <= 3:
        return float('inf')
    
 
    # Pre-calculate g_ij values: entry g[i,j] stores g_ij(T_min)
    g_array = np.full((k - 1, k - 1), 0.0)
    for i in range(1, k - 1):
        for j in range(1, i+1):
            g_array[i, j] = g_ij(i, j, T_min)
            
    # allocate space for distributions Z_i
    # z_array[i,j] = P(Z_i = j)
    z_array = np.zeros((k-1, k-1))
    
    # Set base case:
    z_array[1, 1] = 1
    
    for i in range(2, k-1):
        
        l = i // 2
        r = (i+1) // 2
        
        # calculate sum distribution
        sum_dist = np.convolve(z_array[l, :], z_array[r,:], mode='full')[:k-1]
        
        # calculate new distribution
        z_array[i, :] = sum_dist @ g_array
        

    # sanity checking 
    for i in range(1, z_array.shape[0]):
        assert np.isclose(z_array[i, :].sum(), 1.0)
        
    # calculate our bound using the calculated P(Z_i = 1) 
    try:
        # Calculate sum of log terms
        log_sum = 0
        for ell in range(2, k-1):  # l from 2 to k-2
            
            z_ell = z_array[ell, 1]

            # Add log(1-z_ell) to running sum
            log_term = np.log(1 - z_ell)
            log_sum += log_term
        
        if log_sum == 0:  # Avoid division by zero
            return float('inf')
            
        numerator = np.log((1-q)/(k-3))
        denominator = log_sum / (k-3)
        
        if denominator >= 0:  # Should be negative
            return float('inf')
            
        return numerator / denominator
        
    except (ValueError, ZeroDivisionError):
        return float('inf')