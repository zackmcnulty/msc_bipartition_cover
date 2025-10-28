'''
    Here we collect code for generating many samples from the MSC model and estimating the empirical CDF of the 
    number of gene trees required to achieve a bipartition cover bound. 
'''
import dendropy
from dendropy.simulate import treesim
from dendropy import Tree, TreeList
from tqdm import tqdm
import numpy as np
from scipy.stats import ecdf

def get_T_min(species_tree : Tree) -> float:
    '''
        Finds the minimum length of an internal branch within a species tree.
    '''
    internal_lengths = [edge.length for edge in species_tree.postorder_edge_iter() 
                   if edge.length is not None and not edge.is_terminal()]
    return min(internal_lengths)

def add_flips(bipartition_set):
    '''
        Adds the flipped versions of every bipartition when they are stored as bitmasks strings
        (e.g. AB | CD -> CD | AB is equivalent to 1100 flipping to 0011)
    '''
    mapping = {'0': '1', '1': '0'}
    flips = set()
    for x in bipartition_set:
        flipped = ''.join(mapping[char] for char in x)
        flips.add(flipped)
        
    bipartition_set.update(flips)

def sample_gene_trees_until_cover(species_tree: dendropy.Tree, num_samples=1, max_genes=1000) -> int:
    '''
        Samples gene trees from the given species_tree under the MSC model until
        we get a bipartition cover of species_tree. Repeats this process num_samples times.
        Returns a list arr where arr[i] is the number of gene trees
        required on the ith trial to get a cover, or returns np.nan if the first max_genes trees do
        not form a cover.
    '''
    # create mapping from gene names to species names; e.g. dendropy needs to know where
    # in the species tree each gene lineage begins. Since we have no gene lineage names to
    # begin with, we use this helper function to create names for the gene lineages and 
    # map them to the appropriate species.
    # genes_per_species = 1 # only one gene lineage per leaf/species
    # gene_to_species_map = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
    #                             containing_taxon_namespace=species_tree.taxon_namespace,
    #                             num_contained=genes_per_species)
    gene_to_species_map = {}
    for taxon in species_tree.taxon_namespace:
        gene_to_species_map[taxon] = taxon  # Map each taxon to itself
        
    gene_to_species_map = dendropy.TaxonNamespaceMapping(mapping_dict=gene_to_species_map)

    # Extract all species tree bipartitions
    species_bipartitions = species_tree.encode_bipartitions()
    
    cover_counts = np.zeros(num_samples)
    
    for trial in range(num_samples):
        
        missing_bipartitions = {str(s) for s in species_bipartitions}
        num_genes = 0

        
        while len(missing_bipartitions) > 0:
            num_genes += 1

            # cover not found in allotted time
            if num_genes > max_genes:
                cover_counts[trial] = np.nan
                break

            # Sample new gene tree (simpler approach)
            gene_tree = treesim.contained_coalescent_tree(
                containing_tree=species_tree,
                gene_to_containing_taxon_map=gene_to_species_map
            )

            # Get bipartitions from this gene tree; Since bipartitions can be
            # represented in two different ways, we add both to our collection 
            gene_bipartitions = {str(s) for s in gene_tree.encode_bipartitions()}
            add_flips(gene_bipartitions)

            # Remove any bipartitions we found from missing set
            missing_bipartitions -= gene_bipartitions

        else:
            cover_counts[trial] = num_genes
    
    return cover_counts

########################################################################################
    # Simulating empirical bipartition coverage probabilities
########################################################################################

class InsufficientSamplingError(RuntimeError):
    """Raised when too many simulations hit the max_genes limit, 
    suggesting insufficient sampling for reliable CDF estimation."""
    pass

def get_empirical_coverage_probs(species_tree_generator, k_vals, T_vals, num_samples=100, max_genes=1000, tolerance=0.01):
    '''
        Estimates the distribution of the number of gene trees required to obtain a bipartition cover through Monte Carlo simulation. 
        
        Parameters
        ----------
            species_tree_generator : function
                    Function that takes (k, T_min) and returns a species tree topology
            k_vals : list[int]
                List of species counts to test
            T_vals : list[float]
                List of minimum branch lengths to test
            num_samples : int
                Number of simulation trials per (k, T_min) combination
            max_genes : int
                Maximum gene trees to sample before giving up
            tolerance: float
                The maximum probability allowed for P(count >= max_genes). If this probability is too large
                it suggests many simulations are being pre-maturely stopped at max_genes, and our empirical
                CDF estimate is probably off. 
        
        Returns dictionary with keys (k, T_min) and values ndarrays x of length max_genes representing the empirical CDF of the bipartition cover count.
        Namely, x[n] = P(cover_count <= n).
    '''
    
    total_iterations = len(T_vals) * len(k_vals)
    all_empirical_cdfs = dict()
    
    with tqdm(total=total_iterations, desc="Processing combinations") as pbar:
    
        for t_idx, T_min in enumerate(T_vals):
            for k_idx, k in enumerate(k_vals):
                pbar.set_description(f"T_min={T_min}, k={k}")
                pbar.update(1)

                # generate species tree
                species_tree = species_tree_generator(k, T_min)

                # Collect empirical data
                empirical_counts = sample_gene_trees_until_cover(species_tree, max_genes=max_genes, num_samples=num_samples)

                # Compute empirical CDF
                f = ecdf(empirical_counts)
                #n_vals = np.arange(min(empirical_counts), max(empirical_counts) + 1)
                n_vals = np.arange(1, max_genes + 1)
                empirical_coverage_probs = f.cdf.evaluate(n_vals)
                
                # Raise error if it appears max_genes is too small (too much mass at max_genes)
                if 1 - empirical_coverage_probs[-2] >= tolerance:
                    raise InsufficientSamplingError(
                        f"Exceeded tolerance: {1 - empirical_coverage_probs[-2]:.3f} ≥ {tolerance}. "
                        f"Too many runs reached max_genes={max_genes}. Consider increasing max_genes."
                    )
                
                key = (k, T_min)
                all_empirical_cdfs[key] = empirical_coverage_probs
                
    return all_empirical_cdfs