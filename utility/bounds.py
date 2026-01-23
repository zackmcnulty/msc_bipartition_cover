'''
    In this file we collect functions for calculating our various bipartition cover bounds as
    functions of the relevant parameters: the species count k, minimmum branch length T_min, and the coverage probability q.
'''

import numpy as np
from utility.coalescent_probabilities import g_ij

#################################################################################
# Helpers
#################################################################################


def _sum_powers(a_vals, n):
    return float(np.sum(np.power(a_vals, n)))


def _min_n_by_bisection(a_vals, target, n_min=1, upper_bound=None, max_steps=120):
    """
    Return the smallest integer n >= n_min such that sum(a_vals ** n) <= target.
    Assumes a_vals in [0, 1]. Returns inf if the target is unattainable.
    """
    a_vals = np.asarray(a_vals, dtype=float)
    eps = np.finfo(float).tiny

    if a_vals.size == 0:
        return n_min
    if target < 0:
        return float('inf')
    if np.any(a_vals < 0) or np.any(a_vals > 1):
        return float('inf')

    # Avoid exact 1.0 due to underflow; preserves monotonicity in n.
    a_vals = np.minimum(a_vals, 1.0 - eps)

    if _sum_powers(a_vals, n_min) <= target:
        return n_min

    if upper_bound is not None and np.isfinite(upper_bound):
        hi = max(n_min, int(np.ceil(upper_bound)))
    else:
        hi = max(n_min, 1)

    if _sum_powers(a_vals, hi) > target:
        for _ in range(max_steps):
            hi *= 2
            if _sum_powers(a_vals, hi) <= target:
                break
        else:
            return float('inf')

    low = n_min
    while low < hi:
        mid = (low + hi) // 2
        if _sum_powers(a_vals, mid) <= target:
            hi = mid
        else:
            low = mid + 1

    return low


def _original_upper_bound(k, T_min, q, g_val=None):
    if k <= 3 or T_min <= 0:
        return float('inf')

    if g_val is None:
        g_val = g_ij(k - 2, 1, T_min)

    g_val = min(max(g_val, 0.0), 1.0)
    if g_val <= 0 or g_val >= 1:
        return float('inf')

    numerator = np.log((1 - q) / (k - 3))
    denominator = np.log1p(-g_val)
    if denominator >= 0:
        return float('inf')

    upper = numerator / denominator
    if not np.isfinite(upper) or upper < 0:
        return float('inf')

    return float(np.ceil(upper))

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
        g_val = min(max(g_val, 0.0), 1.0)
        if g_val >= 1 or g_val <= 0:  # Avoid log(0) or log(negative)
            return float('inf')
        
        numerator = np.log((1-q)/(k-3))
        denominator = np.log1p(-g_val)
        
        if denominator >= 0:  # log(1-g) should be negative
            return float('inf')

        value = numerator / denominator
        if not np.isfinite(value) or value < 0:
            return float('inf')
            
        return float(np.ceil(value))
        
    except (ValueError, ZeroDivisionError):
        return float('inf')

#################################################################################
# Caterpillar Bound
#################################################################################

def caterpillar_bound(k, T_min, q=0.95):
    """
    Calculate the improved bound using bisection:
    n = min{n : sum_{l=2}^{k-2} (1 - g_{l,1}(T_min))^n <= 1 - q}.
    """
    if k <= 3 or T_min <= 0:
        return float('inf')
    
    try:
        a_vals = []
        g_k2 = None
        for ell in range(2, k-1):  # l from 2 to k-2
            g_val = g_ij(ell, 1, T_min)
            g_val = min(max(g_val, 0.0), 1.0)
            if ell == k - 2:
                g_k2 = g_val
            a_vals.append(1.0 - g_val)

        upper = _original_upper_bound(k, T_min, q, g_val=g_k2)
        return _min_n_by_bisection(a_vals, 1.0 - q, n_min=1, upper_bound=upper)
        
    except (ValueError, ZeroDivisionError):
        return float('inf')

#################################################################################
# One-Step Bound
#################################################################################

def one_step_bound(k, T_min, q=0.95):
    """
    Calculate the improved bound using bisection:
    n = min{n : sum_{l=2}^{k-2} (1 - q_ell)^n <= 1 - q}.
    """
    if k <= 3 or T_min <= 0:
        return float('inf')
    
    
    # Otherwise, pre-calculate all values of g_ij(T_min) we will need
    g_array = np.full((k, k), np.nan)
    for i in range(1, k):
        for j in range(1, i+1):
            g_array[i, j] = g_ij(i, j, T_min)
    
    try:
        a_vals = []
        g_k2 = g_array[k - 2, 1] if k - 2 >= 1 else None
        for ell in range(2, k-1):  # l from 2 to k-2
            
            q_ell = 0
            
            # calculate q_ell
            m_ell = ell // 2
            for r in range(1, m_ell + 1):
                for s in range(1, ell - m_ell + 1):
                    q_ell += g_array[m_ell, r] * g_array[ell - m_ell, s] * g_array[r + s, 1]

            q_ell = min(max(q_ell, 0.0), 1.0)
            a_vals.append(1.0 - q_ell)

        upper = _original_upper_bound(k, T_min, q, g_val=g_k2)
        return _min_n_by_bisection(a_vals, 1.0 - q, n_min=1, upper_bound=upper)
        
    except (ValueError, ZeroDivisionError):
        return float('inf')

#################################################################################
# Balanced Bound
#################################################################################


def balanced_bound(k, T_min, q=0.95):
    """
    Calculate the improved bound using bisection:
    n = min{n : sum_{l=2}^{k-2} (1 - z_ell)^n <= 1 - q}.
    """
    if k <= 3 or T_min <= 0:
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
        a_vals = []
        for ell in range(2, k-1):  # l from 2 to k-2
            
            z_ell = z_array[ell, 1]
            z_ell = min(max(z_ell, 0.0), 1.0)
            a_vals.append(1.0 - z_ell)

        g_k2 = g_array[k - 2, 1] if k - 2 >= 1 else None
        upper = _original_upper_bound(k, T_min, q, g_val=g_k2)
        return _min_n_by_bisection(a_vals, 1.0 - q, n_min=1, upper_bound=upper)
        
    except (ValueError, ZeroDivisionError):
        return float('inf')
