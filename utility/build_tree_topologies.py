'''
    In this file we define functions used to construct dendropy Trees for various different species tree topologies.
'''
import dendropy
from dendropy.simulate import treesim
from dendropy import Tree, TreeList
from utility.msc_sampling import get_T_min

def create_caterpillar_tree(k, branch_length=1.0):
    """Create a caterpillar tree with k taxa"""
    taxa = [f"T{i}" for i in range(k)]
    taxon_namespace = dendropy.TaxonNamespace(taxa)
    
    # Build newick string for caterpillar: ((((T0,T1),T2),T3),T4)
    newick = taxa[0]
    for i in range(1, k):
        newick = f"({newick},{taxa[i]}):{branch_length}"
    newick += ";"
    
    tree = dendropy.Tree.get(data=newick, schema="newick", 
                           taxon_namespace=taxon_namespace)
    return tree

def build_balanced_newick(taxon_list, branch_length):
    
    '''
        Helper method for create_balanced_tree
    '''
    if len(taxon_list) == 1:
        return taxon_list[0]
    elif len(taxon_list) == 2:
        return f"({taxon_list[0]}:{branch_length},{taxon_list[1]}:{branch_length})"
    else:
        mid = len(taxon_list) // 2
        left = build_balanced_newick(taxon_list[:mid], branch_length)
        right = build_balanced_newick(taxon_list[mid:], branch_length)
        return f"({left},{right}):{branch_length}"

    
def create_balanced_tree(k, branch_length=1.0):
    """Create a balanced tree with k taxa (works best when k is power of 2)"""
    if k == 1:
        return dendropy.Tree.get(data="T0;", schema="newick")
    
    taxa = [f"T{i}" for i in range(k)]

    newick = build_balanced_newick(taxa, branch_length) + ";"
    tree = dendropy.Tree.get(data=newick, schema="newick")
    return tree

def create_yule_tree(num_species, T_min, birth_rate=1):
    '''
        Generates a tree under the Yule model, but renormalizes the edge lengths so that the 
        minimum length is T_min. 
    '''
    yule_tree = treesim.birth_death_tree(
    birth_rate=birth_rate, 
    death_rate=0.0, 
    num_extant_tips=num_species
    )
    
    T = get_T_min(yule_tree)
    for edge in yule_tree.postorder_edge_iter():
        edge.length *= T_min/T
        
    return yule_tree