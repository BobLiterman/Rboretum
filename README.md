## Rboretum: Comparisons and Visualizations of Phylogenetic Tree Topology and Alignment Support in R  
**Dr. Bob Literman**  
**Dr. Rachel Schwartz**  
University of Rhode Island  

v.1.0 (01/01/2020)  
** **
### Program Dependencies
- R version >= 3.6.1
- Python version >= 3.6+ w/Biopython  
** **
### R Package Dependencies  
- ape (>= 5.3)  
- phytools (>= 0.6.99)  
- phangorn (>= 2.5.5)  
- testit (>= 0.9)  
- tidyverse (>= 1.2.1)  
- tools (>= 3.6.1)  
- reticulate (>= 1.13)  
- scales (>= 1.0.0)  
- viridis (>= 0.5.1)  
- cowplot (>= 1.0.0)  
- ggtree (>= 1.16.6)  
- grDevices (>= 3.6.1)  
- viridisLite (>= 0.3.0)  
- rlist (>= 0.4.6.1)  
- BiocManager (>= 1.30.10)  
** **
### Rboretum Installation
```
library(devtools)
devtools::install_github("BobLiterman/Rboretum")
```
** **
### Loading Library
```
library(Rboretum)
source_python(paste(system.file(package="Rboretum"), "Split_Processor.py", sep="/"))
```
** **
Are you visiting this repo after seeing the **totally rad** talk at SSB 2020 by @ErenAda ? Welcome! Just so you know...Rboretum is very much a work in progress, **and raising early issues would be very helpful**. What do you wish you could easily do in R? Let me know and I'll try to make it happen for you.

For basic usage, please see a **brief tutorial script in 'inst' folder**...all to be expanded soon!

Feel free to reach out to me on Twitter (@BobLiterman) or over e-mail (literman@uri.edu)
** **
### Rboretum Function Summary    

#### Tree Functions
1) **readRooted**: Read in tree(s) as rooted phylo or multiPhylo
2) **treeTrimmer**: Prune tree(s) down to desired set of taxa  
3) **checkTips**: Check whether taxon labels are present in tree(s), and if so, optionally, whether they are monophyletic or root groups
4) **getTreeSplits**: Break down rooted phylo or multiPhylo objects into respective sets of splits  
5) **getTreeClades**: Extract monophyletic groups from rooted phylo or multiPhylo objects, and/or tally clade distribution among multiple trees
6) **getSharedTaxa**: Given a multiPhylo, return list of shared taxon labels
7) **getUniqueTopologies**: Reduce multiPhylo down to unique topologies    
8) **summarizeMultiPhylo**: Print summary information about a multiPhylo object  
9) **convertLabels**: Convert taxon labels for trees and Rboretum files  
10) **extractNodeAges**:  Extract node age information from rooted phylo or multiPhylo objects that are scaled to time  
11) **treeNamer**:  Adds dummy names to multiPhylo trees (e.g. Tree_1, Tree_2, etc.)  

#### Alignment Support Functions  
1) **getAlignmentSignal**: For one or more multiple-sequence alignments (NEXUS, FASTA, phylip, phylip-relaxed), get the pattern of sequence variation (e.g. invariant, biallelic, triallelic, etc.) at each site as it relates to a specified set of taxa, as well as missing taxa and gap detection, and parsimony-based split analysis  
2) **getTreeSupport**: Process output from getAlignmentSignal to query support for specific data subsets (e.g. only biallelic sites, or no missing data allowed), for either a specific tree, set of trees, or clade of interest  

#### Tree Plotting Functions  
1) **treePlotter**: ggtree wrapper to plot tree(s) with lots of options  
2) **tandemPlotter**: Plot trees/ggplots side by side  
** **
### Rboretum Function Arguments
#### NOTE:  
**Required arguments are in bold**  
*Optional arguments are in italics*  
**Default options are bold**  
** **
1) **readRooted**: Returns phylo or multiPhylo  
    - **to_root**: Character vector of tree locations. Options include:   
      - Path to single tree  
        - '/path/to/tree.nwk'  
      - Paths to multiple trees  
        - c('/path/to/tree_1.nwk','/path/to/tree_2.nwk')
      - Path to directory containing tree(s)  
        - '/path/to/tree/directory'  
    - **root_taxa**: Character vector of root taxon labels. All trees will be rooted with these taxa  
    - *tree_names*: If reading in multiple trees, character vector of tree names [**Default: No tree names**]  
    - *dummy_names*: TRUE/**FALSE**; Add dummy tree names to multiPhylo [e.g. Tree_1, Tree_2, etc.]
    - *prefix*: If a directory is given, a character vector of tree file prefixes (e.g. Read in only files that start with 'RAxML') [**Default: All files in directory**]  
    - *suffix*: If a directory is given, a character vector of tree file suffixes (e.g. Read in only files that end with 'nexus') [**Default: All files in directory**]    
** **    
2) **treeTrimmer**: Returns phylo or multiPhylo   
    - **tree**: Tree(s) to trim. Options include:  
      - phylo object  
      - multiPhylo object
    - **taxa**: Character vector of taxon labels to retain or remove from tree(s)  
    - *remove*: TRUE/**FALSE**; Remove 'taxa' from tree(s) rather than trimming tree(s) down to only 'taxa'  
** **    
3) **checkTips**: Returns TRUE or FALSE 
    - **tree**: Tree(s) to check. Options include:  
      - phylo object  
      - multiPhylo object
    - **taxa**: Character vector of taxon labels. *Return TRUE if  'taxa' are present in all trees in 'tree', else; FALSE*
    - *check_mono*: TRUE/**FALSE**; Also check if 'taxa' are monophyletic in all supplied trees  
    - *check_root*: TRUE/**FALSE**; Also check if 'taxa' are part of a root split in all supplied trees  
** **
4) **getTreeSplits**: Returns dataframe  
    - **tree**: Tree(s) to split. Options include:  
      - phylo object  
      - multiPhylo object  
** **
5) **getTreeClades**: *Default*: Returns character vector of monophyletic clades from tree(s)  
    - **tree**: Tree(s) from which to extract clades. Options include:  
      - phylo object  
      - multiPhylo object    
    - *include_root*: TRUE/**FALSE**; Also return both sides of root split  
    - *print_counts*: TRUE/**FALSE**; Print a summary table of unique splits and how many/which trees contain them, but return clades  
    - *return_counts*: TRUE/**FALSE**; Instead of clades, return summary table of unique splits and how many/which trees contain them  
** **
6) **getSharedTaxa**: Returns a character vector of taxon labels present in all trees
    - **trees**: multiPhylo object  
** **
7) **getUniqueTopologies**: *Default*: Returns phylo or multiPhylo  
    - **trees**: multiPhylo object  
    - *print_table*: TRUE/**FALSE**; Print summary table of unique topologies  
    - *return_table*: TRUE/**FALSE**; Return summary table of unique topologies  
** **
8) **summarizeMultiPhylo**: No return  
    - **trees**: multiPhylo object to summarize
** **
9) **convertLabels**: Returns object with converted taxon labels  
    - **to_convert**:  Object to convert. Options include:  
      - phylo object  
      - multiPhylo object  
      - Character vector of taxon IDs  
      - List of character vectors containing taxon IDs  
      - Character vector of semicolon-separated clades  
    - **name_df**: Dataframe containing name equivalency columns   
    - *from*: Column name with current taxon tip IDs [**Default: Column 1**]    
    - *to*: Column name with desired taxon tip IDs [**Default: Column 2**]     
** **
10) **extractNodeAges**: Returns a dataframe  
      - **tree**: Tree(s) from which to extract dates. Options include:  
        - phylo object  
        - multiPhylo object  
      - *return_summary*: TRUE/**FALSE**; If a multiPhylo is provided, return summary values for each clade age rather than separate values for each tree  
** **
11) **treeNamer**: Returns named multiPhylo  
    - **trees**: multiPhylo of trees to name  
** **
12) **getAlignmentSignal**: Returns a dataframe  
    - **alignment_path**: Character vector of alignment locations. Options include:   
      - Path to single alignment  
      - Paths to multiple alignments  
      - Path to directory containing alignment(s)  
    - **species_info**: Analyze alignment(s) for the following taxa. Options include:  
      - phylo object  
      - multiPhylo object  
      - Character vector of taxon labels  
    - *use_gaps*: TRUE/**FALSE**; Treat gaps (-) in alignment as potentially informative data  
    - *alignment_name*: Character vector containing a unique name for each alignment file  [**Default**: Use file name]  
    - *prefix*: If a directory is given, a character vector of alignment file prefixes [**Default: All files in directory**]  
    - *suffix*: If a directory is given, a character vector of alignment file suffixes [**Default: All files in directory**] 
    - *existing_signal*: Output from getAlignmentSignal run against the same taxa, but with different alignments. Will append results to dataframe if names are unique.  
** **
13) **getTreeSupport**: *Default*: Returns dataframe  
    - **signal**: Output from getAlignmentSignal run against same taxa as present in 'tree',  or at least with all taxa from 'clade'    
    - **tree**: Tree(s) on which to map support. Options include:  
      - phylo object  
      - multiPhylo object 
      
      OR 
    - **clade**: Character vector containing taxon labels. Get support for this clade (**Note**: supercedes 'tree' if both are provided)  
    - *separate_signal*: **TRUE**/FALSE; If 'signal' contains data from multiple alignments, return results separated by alignment. If FALSE, return total support.  
    - *return_integer*: TRUE/**FALSE**; Return integer vector of total support for each clade across alignments, rather than a dataframe  
    - *include_root*: TRUE/**FALSE**; Get support for root clade  
    - *dataset_name*: Character vector of names for each dataset [**Default: Alignment name + 'm<MAX_MISSING>'**]  
    - *max_missing*: Number of taxon allowed to have missing data in a column and still extract signal [**Default: 0**]  
    - *include_gap*: TRUE/**FALSE**; Use signal from columns where any taxon has a gap (-)  
    - *only_gap*: TRUE/**FALSE**; Only use signal from columns with gaps (-) present    
    - *include_singleton*: **TRUE**/FALSE; Use signal from columns where some taxon are singletons  
    - *include_biallelic*: **TRUE**/FALSE; Use signal from biallelic sites  
    - *include_triallelic*: **TRUE**/FALSE; Use signal from triallelic sites  
    - *include_quadallelic*: **TRUE**/FALSE; Use signal from quadallelic sites  
    - *include_pentallelic*: **TRUE**/FALSE; Use signal from pentallelic sites (column contains A,C,T,G,-)  
    - *return_table*: TRUE/**FALSE**; Return table of support counts for all clades rather than support for particular trees or clades  
    - *existing_support*: Output from getTreeSupport run against the same taxa, but with different alignments/options. Will append results to dataframe if names are unique.
** **
14) **treePlotter**: Returns plot object  
    - **tree**: Tree to plot. Options include:  
      - phylo object  
      - multiPhylo object 
    - *clade_support*: Output of getTreeClades(return_counts=TRUE); Will color node labels based on clade prevalence in multiPhylo   
    - *tree_support*: Output of getTreeSupport(); Will add alignment support information to tree(s) for node labeling or geom_size scaling   
    - *geom_size*: If clade_support or tree_support are provided, how should the geoms be sized? Options include:  
      - Single numeric value [**Default: 4**]  
      - c(min,max): With tree_support, scale total tree_support to fall with min:max   
      - 'log': With tree_support, log-transform total tree_support values  
    - *scale_range*: c(min,max); With tree_support, min and max of values to scale. Raw values below 'min' or above 'max' will be scaled to in-scale min and max values [**Default: Scale all values**]
    - *use_pies*: TRUE/**FALSE**; If tree_support contains data from multiple alignments, plot pie chart on node rather than geom_nodepoint  
    - *pie_xnudge*: xnudge for pie geom [**Default: 0**]  
    - *pie_ynudge*: ynudge for pie geom [**Default: 0**]  
    - *pie_legend_position*: Numerical vector of length four specifying legend xmin,xmax,ymin,ymax respectively. [**Default = c(1,1,1,1)**; **NOTE**: This is a necessary oddity of generating a legend for nodepie data]
    - *branch_length*: TRUE/**FALSE**; Plot tree with branch lengths  
    - *branch_weight*: ggtree branch thickness [**Default: 1**] 
    - *node_label*: ggtree node label choice. Options include:  
      - **"bs": Bootstrap values**  
      - "none": No labels  
      - "node": ggtree node ID  
      - "support": With tree_support, the total raw site support count for each node  
      - "clade": With clade_support, the percent of trees from multiPhylo that contain that node  
    - *node_label_font_size*: Font size for node labels [**Default: 5**]
    - *node_label_fontface*: Fontface for node label. Options include:  
      - **"plain"**  
      - "bold"  
      - "italic"  
      - "bold.italic"  
    - *node_label_color*: Color for node label [**Default: 'black'**]
    - *node_label_box*: TRUE/**FALSE**; Print node label as a label in a box with a white background  
    - *node_label_nudge*: Node label nudge [**Default: 0**]  
    - *taxa_font_size*: Font size for taxon labels [**Default: 5**] 
    - *taxa_fontface*: Fontface for taxon label. Options include:  
      - **"plain"**  
      - "bold"  
      - "italic"  
      - "bold.italic"   
    - *taxa_offset*: Offset value for taxon labels [**Default: 0**]
    - *xmax*: ggplot xlim upper limit (e.g if long tip labels run off plot) [**Default: Don't extend plot**]
    - *reverse_x*: TRUE/**FALSE**; Plot tree in reverse orientation (tips on left)   
    - *to_color*: Tips or clade to color. Options include:  
      - Character vector of taxa, all to be colored the same color  
      - List of taxon vectors, each of which will have their own color. List can be named for use with a legend (set highlight_legend == TRUE)  
    - *colors*:  Colors for clade highlighting. Must be hex or valid R colors.  
    - *highlight_legend*: TRUE/**FALSE**; Include legend for group coloration  
    - *color_branches*: TRUE/**FALSE**; Color branches rather than tip labels [**NOTE**: Not available if clade_support is given]  
    - *plot_title*: Character vector of plot titles [**Default: No name for phylo; tree name for multiPhylo**]
    - *legend_shape_size*: ggplot2 size for legend icons [**Default: 5**]  
    - *legend_font_size*: ggplot2 size for legend font [**Default: 10**]  
    - *legend_title_size*: ggplot size for legend title [**Default: 10**]  
    - *geom_alpha*: ggplot2 alpha value for geom_nodepoint (or pies) [**Default: 0.9**]   
    - *geom_color*: ggplot2 color value for geom_nodepoint if clade_support not provided [**Default: 'red'**]

** **
15) **tandemPlotter** Returns plot object  
    - Separate ggtree/ggplot objects, or a list of object to plot side by side.  

