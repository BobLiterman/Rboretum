## Rboretum: Comparisons and Visualizations of Phylogenetic Tree Topology and Alignment Support in R  
**Dr. Bob Literman**  
**Dr. Rachel Schwartz**  
University of Rhode Island  

v.1.0 (11/16/2019)  

### DEPENDENCIES

R version >= 3.6.1

- ggtree
- ape (>= 5.3)
- phytools (>= 0.6.99)
- phangorn (>= 2.5.5)
- testit (>= 0.9)
- tidyverse (>= 1.2.1)
- purrr
- tools (>= 3.6.1)
- reticulate (>= 1.13)
- scales (>= 1.0.0)
- viridis (>= 0.5.1)
- cowplot (>= 1.0.0)
- Python3
- Biopython

### INSTALLATION  

```
devtools::install_github("BobLiterman/Rboretum")
```

### RBORETUM FUNCTIONS  

#### read.rooted()  
- **Function Description:**  
Reads an unrooted or rooted tree into R as a phylo object, rooted at a specified outgroup  

- **Arguments:**  
      **1) tree_path:** Path to tree file (must be readable by ape::read.tree or ape::read.nexus)  
      **2) root_taxa:** Character vector containing desired outgroup taxa  
      
- **Returns:**  
Rooted phylo object

#### trim.tree()  
- **Function Description:**  
Trims a phylo object or a set of multiPhylo objects down to a desired set of taxa

- **Arguments:**  
      **1) tree:** phylo or multiPhylo object  
      **2) taxa:** Character vector of desired tip labels to keep (or discard if remove=TRUE)  
      **3) remove:** *OPTIONAL:* If TRUE, tip labels specified by 'taxa' are removed from all trees rather than retained [Default: FALSE, prune 'tree' down to 'taxa']  
      
- **Returns:**  
phylo or multiPhylo object with desired taxa retained (or removed)

#### check.tip()  
- **Function Description:**  
Checks whether tree(s) contain specific taxa  

- **Arguments:**  
      **1) tree:** phylo or multiPhylo object  
      **2) taxa:** Character vector of desired tip labels to query
      
- **Returns:**  
TRUE if all taxa in tree/all trees; otherwise, FALSE

#### check.root()  
- **Function Description:**  
Checks whether specified taxa are part of a root split

- **Arguments:**  
      **1) tree:** phylo object  
      **2) taxa:** Character vector of desired tip labels to query
      
- **Returns:**  
TRUE if taxa are part of a root split; otherwise, FALSE

#### convert.tips()  
- **Function Description:**  
Converts tip labels on a phylo object based on a user-supplied name conversion dataframe

- **Arguments:**  
      **1) tree:** phylo object  
      **2) name_df:** Dataframe containing name equivalency columns  
      **3) from:** Column ID with current tip labels  
      **4) to:** Column ID with desired tip labels  
      
- **Returns:**  
phylo object with converted tip labels

#### get.splits()  
- **Function Description:**  
Breaks down a rooted phylo object into all splits

- **Arguments:**  
      **1) tree:** phylo object  
      
- **Returns:**  
Dataframe with split information, node numbers, and bootstrap values (if present)  

#### get.clades()  
- **Function Description:**  
Extracts monophyletic clades from a rooted phylo object, including both root clades  
- **Arguments:**  
      **1) tree:** phylo object  
      
- **Returns:**  
Character vector with all clades, with taxa separated by semicolons  

#### check.shared()  
- **Function Description:**  
Given a multiPhylo object, returns TRUE if all trees share at least three taxa   
- **Arguments:**  
      **1) trees:** multiPhylo object  
      
- **Returns:**  
TRUE if three tip labels are shared in all trees; otherwise, FALSE  

#### get.shared()  
- **Function Description:**  
Given a multiPhylo object, returns taxa that occur in all trees  
- **Arguments:**  
      **1) trees:** multiPhylo object  
      
- **Returns:**  
Character vector of taxa shared in all trees  

#### same.taxa()  
- **Function Description:**  
Given a multiPhylo object, returns TRUE if all trees share all taxa
- **Arguments:**  
      **1) trees:** multiPhylo object  
      
- **Returns:**  
TRUE if all taxa occur in all trees; otherwise, FALSE  

#### same.topology()  
- **Function Description:**  
Given a multiPhylo object, returns TRUE if all trees have identical topology, after pruning to common species set if necessary  
- **Arguments:**  
      **1) trees:** multiPhylo object  
      
- **Returns:**  
TRUE if trees have the same topology, after pruning to common species set if necessary ; otherwise, FALSE  

#### check.comparable()  
- **Function Description:**  
Checks if two trees (1) share at least three species, and (2) have a unique topology  
- **Arguments:**  
      **1) tree_1:** phylo object  
      **2) tree_2:** phylo object
      
- **Returns:**  
TRUE if trees share >= 3 species and have a unique topology; otherwise, FALSE  

#### get.comparable()  
- **Function Description:**  
Given a multiPhylo object, returns a dataframe with information about which trees are comparable (share >= 3 species and have a unique topology)  
- **Arguments:**  
      **1) trees:** multiPhylo object  
      **2) tree_names:** *OPTIONAL:* Character vector containing names for trees in multiPhylo. If missing, tree names will be 'Tree_1', 'Tree_2', etc.  
      **3) return_only_comparable:** *OPTIONAL:* If TRUE, return table will only contain information about comparable trees. [Default: FALSE, return all comparisons]  
      
- **Returns:**  
Dataframe containing pairwise tree comparisons, and if comparable, a semicolon-separated list of shared species  

#### compare.clades()
- **Function Description:**  
Given a multiPhylo object, returns a dataframe with information about all splits,  and what trees they do (and do not) occur in    
- **Arguments:**  
      **1) trees:** multiPhylo object  
      **2) tree_names:** *OPTIONAL:* Character vector containing names for trees in multiPhylo. If missing, tree names will be 'Tree_1', 'Tree_2', etc.  
      **3) return_only_shared:** *OPTIONAL:* If TRUE, return table will only contain information clades present in all trees. [Default: FALSE, return information about all clades]  
      
- **Returns:**  
Dataframe containing all clades from multiPhylo, the clade size, and which trees the clade is found in. If 'return_only_shared'==TRUE, returns a two-column dataframe with common clades their size  

#### alignment.signal()  
- **Function Description:**  
Processes a multiple sequence alignment (FASTA, Phylip, or Nexus) and extracts phylogenetic signal information from each site that has data for at least three taxa.  

**NOTE 1:** Signal can be extracted for only a subset of taxa included in the total alignment by specifying fewer taxa via 'species_info'  

**NOTE 2:** To use this function, the user must first source the python script through reticulate (linked to Python3 with Biopython installed):  
```
source_python(paste(system.file(package="Rboretum"), "Split_Processor.py", sep="/"))
```
- **Arguments:**  
      **1) alignment_path:** Path to alignment file  
      **2) species_info:** Either (1) a tree containing species of interest; or (2) a character vector with the desired taxa IDs. *NAMES MUST MATCH IDs IN ALIGNMENT FILE*  
      **3) informative_gaps:** *OPTIONAL:* If TRUE, gaps (-) in the alignment  will be treated as informative characters. [Default: FALSE, gaps are treated as missing data]  
      **4) alignment_name:** *OPTIONAL:* A name for the alignment. If missing, the file name (without extension) will be used.    
      
- **Returns:**  
Dataframe containing phylogenetic signal contained in each site of the alignment, including the pattern of variation (e.g. invariant, singletons, biallelic, triallelic, etc.), and which clades are supported  

#### tree.support()
- **Function Description:**  
Maps output from alignment.signal() onto a phylo object (i.e. How many sites from this alignment support this tree, and which nodes?). Support from multiple alignments can be added to a single dataframe by including the 'existing_splits' argument  

- **Arguments:**  
      **1) signal:** Output from alignment.signal() run with only+all taxa from 'tree'  
      **2) tree:** Rooted phylo object  
      **3) max_missing:** *OPTIONAL:* Only process signal from sites containing 'max_missing' species with missing data [Default: 0]  
      **4) alignment_name:** *OPTIONAL:* A name for the alignment. If missing, the name from alignment.signal() will be used  
      **5) include_gap:** *OPTIONAL:* If FALSE, do not process signal from sites containing gaps (-) [Default: TRUE, process gap signal]  
      **6) include_biallelic:** *OPTIONAL:* If FALSE, do not process signal from sites with biallelic variation [Default: TRUE, process biallelic sites]  
      **7) include_triallelic:** *OPTIONAL:* If FALSE, do not process signal from sites with triallelic variation [Default: TRUE, process triallelic sites]  
      **8) include_quadallelic:** *OPTIONAL:* If FALSE, do not process signal from sites with quadallelic variation [Default: TRUE, process quadallelic sites]  
      **9) include_pentallelic:** *OPTIONAL:* If FALSE, do not process signal from sites with pentallelic variation [Default: TRUE, process pentallelic sites]  
      **10) only_gap:** *OPTIONAL:* If TRUE, only process data from gap positions [Default: FALSE, process non-gap sites]  
      **11) existing_splits:** *OPTIONAL:* Output from tree.support() run with a different alignment (or other options) against the same tree. Will add a new column to the existing dataframe assuming that the alignment ID is unique  
      
- **Returns:**  
Dataframe containing all splits from 'tree', as well as support values from the specified alignment   

#### clade.support()
- **Function Description:**  
Queries output from alignment.signal() for support for a specific clade (i.e. How many sites this alignment support this clade?).  

- **Arguments:**  
      **1) signal:** Output from alignment.signal() run with only+all taxa from 'tree'  
      **2) clade:** Character vector containing taxa from clade to query  
      **3) max_missing:** *OPTIONAL:* Only process signal from sites containing 'max_missing' species with missing data [Default: 0]  
      **4) include_gap:** *OPTIONAL:* If FALSE, do not process signal from sites containing gaps (-) [Default: TRUE, process gap signal]  
      **5) include_biallelic:** *OPTIONAL:* If FALSE, do not process signal from sites with biallelic variation [Default: TRUE, process biallelic sites]  
      **6) include_triallelic:** *OPTIONAL:* If FALSE, do not process signal from sites with triallelic variation [Default: TRUE, process triallelic sites]  
      **7) include_quadallelic:** *OPTIONAL:* If FALSE, do not process signal from sites with quadallelic variation [Default: TRUE, process quadallelic sites]  
      **8) include_pentallelic:** *OPTIONAL:* If FALSE, do not process signal from sites with pentallelic variation [Default: TRUE, process pentallelic sites]  
      **9) only_gap:** *OPTIONAL:* If TRUE, only process data from gap positions [Default: FALSE, process non-gap sites]  
      **10) as_root:** *OPTIONAL:* If TRUE, limit signal to biallelic sites with no taxa missing [Default: FALSE]  
      
- **Returns:**  
Integer count of sites from alignment.signal() that support the queried clade  

#### basic.treePlot()
- **Function Description:**  
A ggtree wrapper that makes plotting trees 'your way' easier!  

- **Arguments:**  
    **1) tree:** ggtree object  
    **2) branch_length:** *OPTIONAL:* If TRUE, plot tree with branch lengths [Default: FALSE, plot as cladogram]  
    **3) branch_weight:** *OPTIONAL:* Adjust thickness of tree branches  
    **4) node_label:** *OPTIONAL:*  Choice of node labels for tree:  
    - 'bs': Print bootstrap values [Default]  
    - 'node': Print ggtree node id  
    - 'none': No node labels  
    
    **5) node_size:** *OPTIONAL:* Adjust size of node labels  
    **6) node_nudge:** *OPTIONAL:* Adjust ggtree node label 'nudge_x'  
    **7) taxa_size:** *OPTIONAL:* Adjust tip label size  
    **8) taxa_italic:** *OPTIONAL:* If TRUE, print tip labels in italic [Default: FALSE]  
    **9) taxa_align:** *OPTIONAL:* Set tip label alignment [Default: No alignment]  
    - 'left'  
    - 'right'  
    
    **10) taxa_offset:** *OPTIONAL:* Set ggtree tip label offset  
    **11) xmax:** *OPTIONAL:* Add padding to right side of plot (e.g. if tip labels run off the plot)  
    **12) reverse_x:** *OPTIONAL:*  If TRUE, plot tree with tip labels on the left [Default: FALSE, plot tree with tip labels on the right]  
    **13) rename_tips:** *OPTIONAL:* Two column dataframe where column 1 contains current tip labels and column 2 contains desired tip labels  
    
#### support.treePlot()  
- **Function Description:**  
Plots tree like basic.treePlot, but also plots support information on the nodes. Support can come from tree.support(), compare.clades() or both.  

**Note:** Must include either tree_support or clade_support argument; can include both    
 
- **Additional Arguments from those in basic.treePlot:**  
    **1) tree_support:** *OPTIONAL:* Output from tree.support() run for the same tree and one or more alignments. Will add geom_nodepoint labels to nodes with total support across alignments  
    **2) clade_support:** *OPTIONAL:* Output from compare.clades() run with the same tree as part of a multiPhylo. Will color geom_nodepoint geoms based on how many trees in the multiPhylo contain that clade    
    **3) support_scales:** *OPTIONAL:* How to choose size of geom_nodepoint  
    
    - Single numeric argument: geom_nodepoint geoms will all be scaled to this size [Default: 10]  
    - c(x,y): Node geoms for non-zero support nodes will be sized based on total support values that have been rescaled to be between x and y (both numeric values; ie. c(1,10))  
    - 'log': Node geoms for non-zero support nodes will be sized based on total support values that have been log transformed  
      
    **4) node_alpha:** *OPTIONAL:* Adjust geom_nodepoint alpha value [Default: 0.9]  
    **5) node_color:** *OPTIONAL:* Adjust geom_nodepoint color if clade_support is missing [Default: 'red']  
    **6) legend_shape_size:**  Adjust size of legend shapes if clade_support is provided  
    **7) legend_font_size:** Adjust size of legend font if clade_support is provided  
    **8) legend_title_size:** Adjust size of legend title if clade_support is provided  
    **9) node_label:**  
    - 'none': No node labels [Default]  
    - 'bs': Print bootstrap values  
    - 'node': Print ggtree node id  
    - 'support': Print total site support for node (not scaled)  
    

#### pies.treePlot()
- **Function Description:**  
Plots tree like basic.treePlot, but also plots support information on the nodes, broken down by alignment as a pie chart. 

**Note:** Most arguments are identical to support.treePlot, but:  
 
 *1.* 'tree_support' must contain information about 2+ alignments to make a pie  
 *2.* Clade support cannot be mapped (no clade_support)  
 *3.* X-axis cannot be reversed (no reverse_x)  

- **Additional Arguments from those in support.treePlot:**  
  **1) legend_position:** *OPTIONAL:* Due to complications in making a legend for inset pie charts, legend coordinates must be set in the form of c(xmin,xmax,ymin,ymax) or it will likely not be visible in the plot  
  **2) pie_xnudge:** *OPTIONAL:* Set ggtree pie label hjust  
  **3) pie_ynudge:** *OPTIONAL:* Set ggtree pie label vjust  
#### tandem.treePlot()
- **Function Description:**  
Simple wrapper that allows plotting two ggtree/ggplot objects side by side in a single plot  

- **Arguments**  
    **1) plot_1:** ggtree/ggplot object to be plotted on the left   
    **2) plot_2:** ggtree/ggplot object to be plotted on the right  
    
- **Returns:**  
A side-by-side plot with plot_1 and plot_2  

#### tableCount()
- **Function Description:**  
Simple helper function that queries output from table() but returns 0 rather than NA if item is not found

- **Arguments:**  
      **1) search_table:** Output from table()  
      **2) name:** Item to search for in table  
      
- **Returns:**  
Integer count of occurrences of  item from table  
#### semiVector()
- **Function Description:**  
Simple helper function that converts a semicolon-separated string into a character vector  

- **Arguments:**  
      **1) string_to_split:** Semicolon-separated character string  
      
- **Returns:**  
Character vector  
