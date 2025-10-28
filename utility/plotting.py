'''
    Here we collect functions used to generate the comparison plots for our various bounds:
        * plotting bounds as functions of their relevant parameters
        * Plotting improvement ratios
        * Plotting overestimation ratios using empirically estimated coverage probabilities.
'''

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import bisect

import os
from pathlib import Path

workdir = os.environ.get('WORKDIR')
base_dir = Path(workdir)
figures_dir = base_dir / "figures"

from utility.coalescent_probabilities import g_ij
from utility.msc_sampling import get_empirical_coverage_probs
from utility.bounds import original_bound, caterpillar_bound, one_step_bound, balanced_bound

########################################################################################################
    # Plotting g_ij
########################################################################################################

def plot_multiple_T_values(j, T_values, max_i=20, savepath=None):
    """
    Plot g_ij(T) for multiple T values on the same graph
    """
    fig = plt.figure(figsize=(6, 4))
    
    for T_val in T_values:
        i_range = range(j, max_i + 1)
        i_values = list(i_range)
        g_values = [g_ij(i, j, T_val) for i in i_values]
        plt.plot(i_values, g_values, 'o-', linewidth=2, markersize=4, label=f'T={T_val}')
    
    plt.xlabel('i', fontsize=12)
    plt.ylabel(f'$g_{{i{j}}}(T)$', fontsize=12)
    plt.title(f'$g_{{i{j}}}(T)$ vs i, various T values', fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
    
    if savepath:
        fig.savefig(savepath, dpi=600, bbox_inches='tight')
        plt.close(fig)
        print(f"Plot saved to: {savepath}")


########################################################################################################
    # Comparing Bounds 
########################################################################################################

def plot_bound(bound, T_vals, k_vals, q=0.99, make_plot=False, log_plot=False, figsize=(10,6), savepath=None, bound_name=None):
    '''
        Plots the given bound as both a function of k, drawing one line for each of the provided T values, and
        as a function of T, drawing one line for each provided k. 
        
        k_vals: k_vals to plot
        T_vals: T-vals to plot
        bound: function f(k, T, q) that outputs the given bound for the set parameters
        make_plot: whether to plot results or not
        log_plot: whether to use a log plot or not. 
    '''
    
    fig, axs = plt.subplots(1,2, figsize=figsize)
    
    ax = axs[0]
    
    # Plot bound as a function of k ###################################
    # Calculate bounds for all combinations
    results = []
    for T in T_vals:
        for k in k_vals:
            try:
                bound_val = bound(k, T, q)
                results.append({'k': k, 'T': T, 'bound': bound_val})
            except Exception as e:
                print(f"Error computing bound for k={k}, T={T}: {e}")
                results.append({'k': k, 'T': T, 'bound': np.nan})
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    if make_plot:
        # Set up the plot style
        sns.set_style("whitegrid")
        
        # Create green color palette - darker for larger T values (longer times)
        # Reverse the palette so smaller T values are lighter
        colors = sns.color_palette("Greens", n_colors=len(T_vals))
        colors = colors[::-1]  # Reverse so smaller T (shorter times) are lighter
        
        # Plot each T value as a separate line
        for i, T in enumerate(T_vals):
            T_data = df[df['T'] == T]
            
            if log_plot:
                ax.semilogy(T_data['k'], T_data['bound'], 
                           color=colors[i], linewidth=2, marker='o', markersize=4,
                           label=f'T = {T}')
            else:
                ax.plot(T_data['k'], T_data['bound'], 
                        color=colors[i], linewidth=2, marker='o', markersize=4,
                        label=f'T = {T}')
        
        # Customize the plot
        ax.set_xlabel('Number of Species (k)', fontsize=12)
        ax.set_ylabel('Bound on Number of Gene Trees (n)', fontsize=12)
        
        if bound_name is None:
            ax.set_title(f'Bipartition Cover Bound vs Number of Species (q = {q})', fontsize=14)
        else:
            ax.set_title(f'{bound_name} Bound vs Number of Species (q = {q})', fontsize=14)
        
        # Add legend with custom styling
        legend = ax.legend(title='Branch Length (T)', title_fontsize=11, 
                          fontsize=10, loc='upper left')
        legend.get_title().set_fontweight('bold')
        
        # Set grid and layout
        ax.grid(True, alpha=0.3)

    
    # Plot bound as a function of T ###########################
    
    ax = axs[1]
    
    # Calculate bounds for all combinations
    results = []
    for T in T_vals:
        for k in k_vals:
            try:
                bound_val = bound(k, T, q)
                results.append({'k': k, 'T': T, 'bound': bound_val})
            except Exception as e:
                print(f"Error computing bound for k={k}, T={T}: {e}")
                results.append({'k': k, 'T': T, 'bound': np.nan})
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    if make_plot:
        # Set up the plot style
        sns.set_style("whitegrid")
        
        # Create blue color palette - darker for larger k values (more species)
        # Reverse the palette so smaller k values are lighter
        colors = sns.color_palette("Blues", n_colors=len(k_vals))
        colors = colors[::-1]  # Reverse so smaller k (fewer species) are lighter
        
        # Plot each k value as a separate line
        for i, k in enumerate(k_vals):
            k_data = df[df['k'] == k]
            
            if log_plot:
                ax.semilogy(k_data['T'], k_data['bound'], 
                           color=colors[i], linewidth=2, marker='s', markersize=4,
                           label=f'k = {k}')
            else:
                ax.plot(k_data['T'], k_data['bound'], 
                        color=colors[i], linewidth=2, marker='s', markersize=4,
                        label=f'k = {k}')
        
        # Customize the plot
        ax.set_xlabel('Branch Length (T)', fontsize=12)
        ax.set_ylabel('Bound on Number of Gene Trees (n)', fontsize=12)
        if bound_name is None:
            ax.set_title(f'Bipartition Cover Bound vs Branch Length (q = {q})', fontsize=14)
        else:
            ax.set_title(f'{bound_name} Bound vs Branch Length (q = {q})', fontsize=14)

        
        # Add legend with custom styling
        legend = ax.legend(title='Number of Species (k)', title_fontsize=11, 
                          fontsize=10, loc='upper right')
        legend.get_title().set_fontweight('bold')
        
        # Set grid and layout
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        
        # Show the plot
        plt.show()
        
    # Save plot if path provided
    if savepath:
        fig.savefig(savepath, dpi=600, bbox_inches='tight')
        plt.close(fig)
        print(f"Plot saved to: {savepath}")
    
    return df

########################################################################################################
    # Improvement Ratios 
########################################################################################################


def plot_improvement_ratio(old_bound, new_bound, T_vals, k_vals, q=0.99, make_plot=True, log_plot=False, figsize=(10,6), savepath=None, old_name=None, new_name=None):
    '''
        Plots the improvement ratio old_bound / new_bound as a function of k,
        drawing one line for each of the provided T values.
        
        old_bound: function f(k, T, q) that outputs the old bound
        new_bound: function f(k, T, q) that outputs the new improved bound
        T_vals: T-vals to plot
        k_vals: k_vals to plot
        q: confidence level parameter
        make_plot: whether to plot results or not
        log_plot: whether to use a log plot or not
    '''
    
    fig, axs = plt.subplots(1,2, figsize=figsize)
    ax = axs[0]
    
    # Calculate improvement ratios for all combinations
    results = []
    for T in T_vals:
        for k in k_vals:
            try:
                old_val = old_bound(k, T, q)
                new_val = new_bound(k, T, q)
                
                # Calculate improvement ratio (how many times better the new bound is)
                if new_val > 0:
                    ratio = old_val / new_val
                else:
                    ratio = np.nan
                    
                results.append({
                    'k': k, 
                    'T': T, 
                    'old_bound': old_val,
                    'new_bound': new_val,
                    'improvement_ratio': ratio
                })
            except Exception as e:
                print(f"Error computing bounds for k={k}, T={T}: {e}")
                results.append({
                    'k': k, 
                    'T': T, 
                    'old_bound': np.nan,
                    'new_bound': np.nan,
                    'improvement_ratio': np.nan
                })
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    if make_plot:
        # Set up the plot style
        sns.set_style("whitegrid")
        
        # Create green color palette - darker for larger T values (longer times)
        # Reverse the palette so smaller T values are lighter
        colors = sns.color_palette("Greens", n_colors=len(T_vals))
        colors = colors[::-1]  # Reverse so smaller T (shorter times) are lighter
        
        # Plot each T value as a separate line
        for i, T in enumerate(T_vals):
            T_data = df[df['T'] == T]
            
            if log_plot:
                ax.semilogy(T_data['k'], T_data['improvement_ratio'], 
                           color=colors[i], linewidth=2, marker='o', markersize=4,
                           label=f'T = {T}')
            else:
                ax.plot(T_data['k'], T_data['improvement_ratio'], 
                        color=colors[i], linewidth=2, marker='o', markersize=4,
                        label=f'T = {T}')
        
        # Add horizontal line at y=1 to show "no improvement" baseline
        ax.axhline(y=1, color='red', linestyle='--', alpha=0.7, 
                   label='No improvement (ratio = 1)')
        
        # Customize the plot
        ax.set_xlabel('Number of Species (k)', fontsize=12)
        ax.set_ylabel('Improvement Ratio (Old Bound / New Bound)', fontsize=12)
        
        if old_name is None or new_name is None:
            ax.set_title(f'Bound Improvement Ratio vs Number of Species (q = {q})', fontsize=14)
        else:
            ax.set_title(f'{old_name} / {new_name} vs Number of Species (q = {q})', fontsize=14)
        
        # Add legend with custom styling
        legend = ax.legend(title='Branch Length (T)', title_fontsize=11, 
                          fontsize=10, loc='upper left')
        legend.get_title().set_fontweight('bold')
        
        # Set grid and layout
        ax.grid(True, alpha=0.3)
        
        # Plot as function of T ###################################
        ax = axs[1]

        # Create a different color palette for k values (use blues to distinguish from greens)
        k_colors = sns.color_palette("Blues", n_colors=len(k_vals))
        k_colors = k_colors[::-1]  # Reverse so smaller k values are lighter
        
        # Plot each k value as a separate line
        for i, k in enumerate(k_vals):
            k_data = df[df['k'] == k]
            
            if log_plot:
                ax.semilogy(k_data['T'], k_data['improvement_ratio'], 
                           color=k_colors[i], linewidth=2, marker='s', markersize=4,
                           label=f'k = {k}')
            else:
                ax.plot(k_data['T'], k_data['improvement_ratio'], 
                        color=k_colors[i], linewidth=2, marker='s', markersize=4,
                        label=f'k = {k}')
        
        # Add horizontal line at y=1 to show "no improvement" baseline
        ax.axhline(y=1, color='red', linestyle='--', alpha=0.7, 
                   label='No improvement (ratio = 1)')
        
        # Customize the plot
        ax.set_xlabel('Branch Length (T)', fontsize=12)
        ax.set_ylabel('Improvement Ratio (Old Bound / New Bound)', fontsize=12)
        if old_name is None or new_name is None:
            ax.set_title(f'Bound Improvement Ratio vs Branch Length (q = {q})', fontsize=14)
        else:
            ax.set_title(f'{old_name} / {new_name} vs Branch Length (q = {q})', fontsize=14)
        
        # Add legend with custom styling
        legend = ax.legend(title='Number of Species (k)', title_fontsize=11, 
                          fontsize=10, loc='upper right')
        legend.get_title().set_fontweight('bold')
        
        # Set grid and layout
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
        
                        # Save plot if path provided
        if savepath:
            fig.savefig(savepath, dpi=600, bbox_inches='tight')
            plt.close(fig)
            print(f"Plot saved to: {savepath}")



##########################################################################################
    # Overestimation Plots
##########################################################################################
def make_overestimation_plot_vs_q(empirical_cdfs, bound_function, log_plot=True, 
                           figsize=(15, 5), savepath=None, bound_name=None, tree_name=None):
    """
    Create overestimation plots comparing empirical coverage to theoretical bounds.
    
    Parameters:
    -----------
    empirical_cdfs: dictionary mapping (k, T_min) -> ndarray
        Contains the empirical cdfs estimated via function get_empirical_coverage_probabilities
    bound_function : function  
        Function that takes (k, T_min, q) and returns theoretical bound
    log_plot : bool
        Whether to use log scale for y-axis
    figsize : tuple
        Figure size for subplots
    save_path : str, optional
        Filepath to save the figure
    """
    # Extract T_vals and k_vals from empirical_cdfs
    T_vals = sorted(list(set([t for k,t in empirical_cdfs.keys()])))
    k_vals = sorted(list(set([k for k,t in empirical_cdfs.keys()])))
    
    # Set up subplots - one for each T_min value
    fig, axes = plt.subplots(1, len(T_vals), figsize=figsize, sharey=True)
    if len(T_vals) == 1:
        axes = [axes]
    
    # Color palette for different k values
    colors = sns.color_palette("Set1", n_colors=len(k_vals))
    
    
    for t_idx, T_min in enumerate(T_vals):
        ax = axes[t_idx]

        for k_idx, k in enumerate(k_vals):

            key = (k, T_min)
            empirical_coverage_probs = empirical_cdfs[key]

            # Calculate theoretical bounds for these probability levels
            # Use the scaling relationship to avoid recomputing bounds
            q0 = 0.5
            initial_bound = bound_function(k, T_min, q0)

            def scaled_bound(q):
                return initial_bound * np.log((1-q) / (k-3)) / np.log((1-q0) / (k-3))

            # Filter out q values too close to 1 and 0
            valid_mask = (empirical_coverage_probs < 0.99) * (empirical_coverage_probs > 0.1)
            q_vals = empirical_coverage_probs[valid_mask]
            n_vals = np.arange(1, len(empirical_coverage_probs) + 1)
            n_vals_valid = n_vals[valid_mask]

            theoretical_bounds = scaled_bound(q_vals)

            # Compute overestimation ratios
            overestimate_ratios = theoretical_bounds / n_vals_valid

            # Plot
            if log_plot:
                ax.semilogy(q_vals, overestimate_ratios, 
                           color=colors[k_idx], marker='o', markersize=3,
                           label=f'k = {k}', linewidth=2)
            else:
                ax.plot(q_vals, overestimate_ratios, 
                       color=colors[k_idx], marker='o', markersize=3,
                       label=f'k = {k}', linewidth=2)

        # Add horizontal line at ratio = 1 (perfect prediction)
        ax.axhline(y=1, color='black', linestyle='--', alpha=0.5, linewidth=1)

        # Customize subplot
        ax.set_xlabel('Probability (q)', fontsize=12)
        ax.set_title(f'$T_{{min}} = {T_min}$', fontsize=14)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)

        # Add legend to rightmost plot
        if t_idx == len(T_vals) - 1:
            ax.legend(title='Species', bbox_to_anchor=(1.05, 1), loc='upper left')
                

        # Set y-label on leftmost plot
    axes[0].set_ylabel('$n_b / n_s$', fontsize=12)
    
    # Overall title
    if bound_name is None or tree_name is None:
        fig.suptitle('Bound Overestimation vs Empirical Coverage Probabilities', fontsize=16, y=0.98)
    else:
        fig.suptitle(f'{bound_name} Bound Overestimation for {tree_name} Trees', fontsize=16, y=0.98)
        
    
    plt.tight_layout()
    
    if savepath:
        fig.savefig(savepath, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {savepath}")
    
    plt.show()




def make_overestimation_plot_vs_T_k(empirical_cdfs, bound_function, q=0.95,
                                    log_plot=True, figsize=(15, 5), savepath=None, bound_name=None, tree_name=None):
    """
    Create overestimation plots comparing empirical coverage to theoretical bounds. Plots over-estimation
    for the given quantile q.
    
    Parameters:
    -----------
    empirical_cdfs: dictionary mapping (k, T_min) -> ndarray
        Contains the empirical cdfs estimated via function get_empirical_coverage_probabilities
    bound_function : function  
        Function that takes (k, T_min, q) and returns theoretical bound
    log_plot : bool
        Whether to use log scale for y-axis
    figsize : tuple
        Figure size for subplots
    save_path : str, optional
        Path to save the figure
    """
    # Extract T_vals and k_vals from empirical_cdfs
    T_vals = sorted(list(set([t for k,t in empirical_cdfs.keys()])))
    k_vals = sorted(list(set([k for k,t in empirical_cdfs.keys()])))
    
    # Set up subplots - one for each T_min value
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
    if len(T_vals) == 1:
        axes = [axes]
    
    # Color palette for different k values
    colors = sns.color_palette("Set1", n_colors=len(k_vals))
    
    total_iterations = len(T_vals) * len(k_vals)
    results = []
    
    
    for t_idx, T_min in enumerate(T_vals):
        for k_idx, k in enumerate(k_vals):

            key = (k, T_min)
            empirical_coverage_probs = empirical_cdfs[key] 
            
            # Compute empirical quantile
            n_vals = np.arange(1, len(empirical_coverage_probs) + 1)
            empirical_quantile = n_vals[bisect.bisect_right(empirical_coverage_probs, q)]

            # Calculate theoretical bounds for these probability levels
            theoretical_bound = bound_function(k, T_min, q=q)

            # Compute overestimation ratios
            overestimation_ratio = theoretical_bound / empirical_quantile

            results.append({
                'k': k, 
                'T': T_min, 
                'theoretical': theoretical_bound,
                'empirical': empirical_quantile,
                'overestimation_ratio': overestimation_ratio
            })

    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Set up the plot style
    sns.set_style("whitegrid")

    # Plot as function of k ##################################################
    ax = axes[0]
    
    # Create green color palette - darker for larger T values (longer times)
    # Reverse the palette so smaller T values are lighter
    colors = sns.color_palette("Greens", n_colors=len(T_vals))
    colors = colors[::-1]  # Reverse so smaller T (shorter times) are lighter

    # Plot each T value as a separate line
    for i, T in enumerate(T_vals):
        T_data = df[df['T'] == T]

        if log_plot:
            ax.semilogy(T_data['k'], T_data['overestimation_ratio'], 
                       color=colors[i], linewidth=2, marker='o', markersize=4,
                       label=f'T = {T}')
        else:
            ax.plot(T_data['k'], T_data['overestimation_ratio'], 
                    color=colors[i], linewidth=2, marker='o', markersize=4,
                    label=f'T = {T}')

    # Add horizontal line at y=1 to show "no improvement" baseline
    ax.axhline(y=1, color='red', linestyle='--', alpha=0.7, 
               label='No Overestimation (ratio = 1)')

    # Customize the plot
    ax.set_xlabel('Number of Species (k)', fontsize=12)
    ax.set_ylabel('Overestimation Ratio (Old Bound / New Bound)', fontsize=12)
    #ax.set_title(f'Bound Overestimation Ratio vs Number of Species (q = {q})', fontsize=14)


    # Add legend with custom styling
    legend = ax.legend(title='Branch Length (T)', title_fontsize=11, 
                      fontsize=10, loc='upper left')
    legend.get_title().set_fontweight('bold')

    # Set grid and layout
    ax.grid(True, alpha=0.3)

    # Plot as function of T ###################################
    ax = axes[1]

    # Create a different color palette for k values (use blues to distinguish from greens)
    k_colors = sns.color_palette("Blues", n_colors=len(k_vals))
    k_colors = k_colors[::-1]  # Reverse so smaller k values are lighter

    # Plot each k value as a separate line
    for i, k in enumerate(k_vals):
        k_data = df[df['k'] == k]

        if log_plot:
            ax.semilogy(k_data['T'], k_data['overestimation_ratio'], 
                       color=k_colors[i], linewidth=2, marker='s', markersize=4,
                       label=f'k = {k}')
        else:
            ax.plot(k_data['T'], k_data['overestimation_ratio'], 
                    color=k_colors[i], linewidth=2, marker='s', markersize=4,
                    label=f'k = {k}')

    # Add horizontal line at y=1 to show "no improvement" baseline
    ax.axhline(y=1, color='red', linestyle='--', alpha=0.7, 
               label='No Overestimation (ratio = 1)')

    # Customize the plot
    ax.set_xlabel('Branch Length (T)', fontsize=12)
    ax.set_ylabel('Overestimation Ratio (Old Bound / New Bound)', fontsize=12)
    #ax.set_title(f'Bound Overestimation Ratio vs Branch Length (q = {q})', fontsize=14)

    # Add legend with custom styling
    legend = ax.legend(title='Number of Species (k)', title_fontsize=11, 
                      fontsize=10, loc='upper right')
    legend.get_title().set_fontweight('bold')

    # Set grid and layout
    ax.grid(True, alpha=0.3)
    
    # Overall title
    if bound_name is None or tree_name is None:
        fig.suptitle('Bound Overestimation vs Empirical Coverage Probabilities', fontsize=16, y=0.95)
    else:
        fig.suptitle(f'{bound_name} Bound Overestimation for {tree_name} Trees', fontsize=16, y=0.95)

    plt.tight_layout()
    plt.show()

    # Save plot if path provided
    if savepath:
        fig.savefig(savepath, dpi=600, bbox_inches='tight')
        plt.close(fig)
        print(f"Plot saved to: {savepath}")



def make_all_overestimation_plots(species_tree_generator, tree_name, k_vals, T_vals, q, all_bounds, names,
                                  num_samples=1000, max_genes=1000, tolerance=0.01):
    '''
        Helper function for joining all our plotting functions above.
    '''
    # Get empirical coverage probabilities
    empirical_cdfs = get_empirical_coverage_probs(species_tree_generator, k_vals, T_vals, num_samples=num_samples, max_genes=max_genes, tolerance=tolerance)

    # Plot overestimation as a function of T_min, k, q
    for name, bound in zip(names, all_bounds):
        print('#' * 80)
        print(f'Results for {name} bound')
        print('#' * 80)
        savepath = figures_dir / 'overestimation_ratios' / f'{tree_name.lower()}_tree_{name.lower()}_bound_vs_T_k.png'
        make_overestimation_plot_vs_T_k(empirical_cdfs, bound, 
                                   log_plot=True, 
                                   q = q,
                                   figsize=(15, 5), savepath=savepath,
                                   bound_name=name, tree_name=tree_name)
        savepath = figures_dir / 'overestimation_ratios' / f'{tree_name.lower()}_tree_{name.lower()}_bound_vs_q.png'
        make_overestimation_plot_vs_q(empirical_cdfs, bound, 
                                       log_plot=True, 
                                       figsize=(15, 5), savepath=savepath,
                                       bound_name=name, tree_name=tree_name)
    