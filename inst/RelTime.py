from Bio import Phylo
import numpy as np
from io import StringIO
import ete3
import tempfile

# From Magazzù Giuseppe (https://github.com/biopython/biopython/issues/1974)
def to_ete3(tree):
    with tempfile.NamedTemporaryFile(mode="w") as tmp:
        Phylo.write(tree, tmp, 'newick')
        tmp.flush()
        return ete3.Tree(tmp.name,format=1)
    
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

def rel_time_AJH(tree):
    """
    This closely (exactly?) follows the original implementation but note that the explanation
    of the original algorithm fails to mention what happens with zero length branches which of course
    give zero division errors
    """
    depth_dict = tree.depths(unit_branch_lengths=True)
    for key in tree.get_terminals():
        del depth_dict[key]

    inv_depth_dict = {}
    for key,val in depth_dict.items():
        try:
            inv_depth_dict[val].append(key)
        except KeyError:
            inv_depth_dict[val] = [key]


    for depth in range(max(list(inv_depth_dict.keys())), -1, -1):
        for clade in inv_depth_dict[depth]:
            temp = clade.depths()
            lens = [temp[term]-clade.branch_length for term in clade.get_terminals()]
            for ds_clade in clade.clades:
                ds_lens = [temp[term]-clade.branch_length for term in ds_clade.get_terminals()]
                if np.mean(lens) > 0:
                    ds_clade.rate = np.mean(ds_lens) / np.mean(lens)
                else:
                    ds_clade.rate = 0.
                for all_ds in ds_clade.get_terminals() + ds_clade.get_nonterminals():
                    if all_ds == ds_clade:
                        pass
                    else:
                        all_ds.rate = all_ds.rate*ds_clade.rate
    return tree

# RelTime-like implementation from Dr. Adam J. Hockenberry
def RelTime_AJH(my_tree, modifier=1.):
    my_tree.add_features(rate=modifier)
    if len(my_tree.get_children()) == 0:
        return
    else:
        l_child = my_tree.children[0]
        r_child = my_tree.children[1]
        l_side_dists = [i.dist*len(i.get_leaves()) for i in l_child.traverse()]
        r_side_dists = [i.dist*len(i.get_leaves()) for i in r_child.traverse()]
        

        total_modifier = np.sum(l_side_dists+r_side_dists)/len(my_tree.get_leaves())
        l_modifier = np.sum(l_side_dists)/len(l_child.get_leaves())/total_modifier * modifier
        r_modifier = np.sum(r_side_dists)/len(r_child.get_leaves())/total_modifier * modifier

        my_tree = RelTime_AJH(l_child, l_modifier)
        my_tree = RelTime_AJH(r_child, r_modifier)

    return

def getRelTimeTree_Orig(tree_string):
    
    # Read in tree from R as string
    pre_rooted_tree = Phylo.read(StringIO(tree_string),format='newick',rooted=True)
    
    # Root tree at midpoint between root clades from https://github.com/adamhockenberry/dca-weighting/tree/master/Code/supporting_functions.py
    rerooted_tree = final_root(clade1="Root_B",clade2="Root_A",max_distance=pre_rooted_tree.distance('Root_A','Root_B'),tree=pre_rooted_tree)
    
    # Run RelTime from https://github.com/adamhockenberry/dca-weighting/tree/master/Code/supporting_functions.py
    rerooted_tree.root.branch_length = 0.0
    rerooted_tree = rel_time_AJH(rerooted_tree)
    for node in rerooted_tree.get_terminals() + rerooted_tree.get_nonterminals():
        if node == rerooted_tree.root:
            continue
        node.branch_length = node.branch_length/node.rate
    rerooted_tree.root.branch_length = None
    ete_tree = to_ete3(rerooted_tree)
    
    return(ete_tree.write())
    
def getRelTimeTree_New(tree_string):
    
    # Read in tree from R as string       
    r_tree = ete3.Tree(tree_string,format=1)
    RelTime_AJH(r_tree)
    
    # Convert tree to chronogram
    for i in r_tree.get_descendants():
      i.dist = i.dist/i.rate
    
    return(r_tree.write())
