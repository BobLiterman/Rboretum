from Bio import Phylo
import numpy as np
from io import StringIO
from ete3 import Tree as EteTree
import tempfile

# From Magazzù Giuseppe (https://github.com/biopython/biopython/issues/1974)
def to_ete3(tree):
    with tempfile.NamedTemporaryFile(mode="w") as tmp:
        Phylo.write(tree, tmp, 'newick')
        tmp.flush()
        return EteTree(tmp.name,format=1)
    
# RelTime and rooting code from https://github.com/adamhockenberry/dca-weighting/tree/master/Code/supporting_functions.py
def final_root(clade1, clade2, max_distance, tree):
    """
    Given the target clades, this actually re-roots the tree between them at the midpoint.
    Input(s):
    clade1 - Bio.Phylo clade object that belongs to tree
    clade2 - Bio.Phylo clade object that belongs to tree
    max_distance - the numerical distance between the two clades
    tree - Bio.Phylo tree object to be re-rooted
    Output(s):
    tree - the tree object, re-rooted
    """
    
    #Depth to go from the ingroup tip toward the outgroup tip 
    root_remainder = 0.5 * (max_distance - (tree.root.branch_length or 0)) 
    #This better be the case
    assert root_remainder >= 0 
    #Crawl between the nodes to find the middle branch
    for node in tree.get_path(clade2): 
        root_remainder -= node.branch_length 
        if root_remainder < 0: 
            outgroup_node = node 
            outgroup_branch_length = -root_remainder
            break 
    else: 
        raise ValueError("Somehow, failed to find the midpoint!")
    #Specifying the outgroup_branch_length directly with this flag lead to some
    #error-prone behavior so I'm doing it in two steps. Must be a bug and/or mis-understanding
    #in Bio.Phylo
    tree.root_with_outgroup(outgroup_node, outgroup_branch_length=0.0)
    assert outgroup_node == tree.root.clades[1]
    tree.root.clades[0].branch_length = tree.root.clades[0].branch_length + root_remainder
    tree.root.clades[1].branch_length = tree.root.clades[1].branch_length - root_remainder
    return tree

def RelTime_AJH(my_tree, modifier=1.):
    """
    This was my first interpretation of the algorithm and it both works and gives really similar
    but not quite exact results for the 3-taxon tree. Both algorithms work identically on the tree
    given here:
    
    https://www.pnas.org/content/109/47/19333
    """
    my_tree.add_features(rate=modifier)
    if len(my_tree.get_children()) == 0:
        return
    else:
        l_child = my_tree.children[0]
        r_child = my_tree.children[1]
        l_side_dists = [i.dist*len(i.get_leaves()) for i in l_child.traverse()]
        r_side_dists = [i.dist*len(i.get_leaves()) for i in r_child.traverse()]
        
        #########################################################
        ###Distinction occurs here
        total_modifier = np.sum(l_side_dists+r_side_dists)/len(my_tree.get_leaves())
        l_modifier = np.sum(l_side_dists)/len(l_child.get_leaves())/total_modifier * modifier
        r_modifier = np.sum(r_side_dists)/len(r_child.get_leaves())/total_modifier * modifier
        #########################################################
        
        my_tree = RelTime_AJH(l_child, l_modifier)
        my_tree = RelTime_AJH(r_child, r_modifier)
    return

def getRelTimeTree(tree_string):
    
    # Read in tree from R as string
    pre_rooted_tree = Phylo.read(StringIO(tree_string),format='newick',rooted=True)
    
    # Root tree at midpoint between root clades from https://github.com/adamhockenberry/dca-weighting/tree/master/Code/supporting_functions.py
    rerooted_tree = final_root(clade1="Root_B",clade2="Root_A",max_distance=pre_rooted_tree.distance('Root_A','Root_B'),tree=pre_rooted_tree)
    
    # Run RelTime from https://github.com/adamhockenberry/dca-weighting/tree/master/Code/supporting_functions.py
    rerooted_tree.root.branch_length = 0.0
    RelTime_AJH(rerooted_tree)
    ete_tree = to_ete3(rerooted_tree)
    
    return(ete_tree.write())
