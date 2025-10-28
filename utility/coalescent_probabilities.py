'''
    In this file we collect several helper functions for calculating the coalescent probabilities g_{ij}(T). 
'''

import numpy as np

def falling_factorial(i, k):
    '''
        Calculates falling factorial i * (i-1) * ... *(i-k+1)
    '''
    if k < 0:
        raise ValueError(f'Negative k = {k} passed')
    if k == 0:
        return 1
    
    return falling_factorial(i-1, k-1) * i

def rising_factorial(i, k):
    '''
        Calculates rising factorial i * (i+1) * ... * (i+k-1)
    '''
    if k < 0:
        raise ValueError(f'Negative k = {k} passed')
    if k == 0:
        return 1
    
    return rising_factorial(i+1, k-1) * i


def g_ij(i, j, T):
    """
    Calculate g_ij(T) using the iterative scheme to avoid overflow.
    
    g_ij(T) = sum_{k=j}^i [exp(-binom(k,2)*T) * (2k-1) * (-1)^(k-j) * c_k] / j!
    
    where c_k = j_(k-1) * i_[k] / ((k-j)! * i_(k))
    and c_{k+1} / c_k = (k+j-1) / (k-j+1) * (i-k) / (i+k)
    
    Parameters:
    i (int): Upper index (number of initial lineages)
    j (int): Lower index (target number of lineages)  
    T (float): Time parameter (branch length in coalescent units)
    
    Returns:
    float: The value of g_ij(T)
    """
    if j > i or i <= 0 or j <= 0:
        return 0.0
    
    if i == j:
        return np.exp(-i * (i - 1) * T / 2)
    
    # Calculate initial c_j for k = j
    # c_j = j_(j-1) * i_[j] / ((j-j)! * i_(j))
    #     = j_(j-1) * i_[j] / ( i_(j))
    j_rising_j_minus_1 = rising_factorial(j, j - 1)
    i_falling_j = falling_factorial(i, j)
    i_rising_j = rising_factorial(i, j)
    
    c_k = j_rising_j_minus_1 * i_falling_j / i_rising_j
    
    # Calculate the sum
    result = 0.0
    j_factorial = falling_factorial(j, j)
    
    for k in range(j, i + 1):
        # Calculate the term
        binomial_coeff = k * (k - 1) // 2
        exp_term = np.exp(-binomial_coeff * T)
        sign_term = (-1) ** ((k - j) % 2)
        k_term = (2 * k - 1)
        
        term = exp_term * k_term * sign_term * c_k / (j_factorial)
        result += term
        
        # Update c_k for next iteration using the recurrence relation
        # c_{k+1} / c_k = (k+j-1) / (k-j+1) * (i-k) / (i+k)
        ratio = ((k + j - 1) / (k - j + 1)) * ((i - k) / (i + k))
        c_k *= ratio
    
    return result