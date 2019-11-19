## Rboretum: Comparisons and Visualizations of Phylogenetic Tree Topology and Alignment Support in R  
**Dr. Bob Literman**  
**Dr. Rachel Schwartz**  
University of Rhode Island  

v.1.0 (11/16/2019)  

### DEPENDENCIES
- ggtree
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
- Python3
- Biopython

### INSTALLATION  
1) Clone repo  

```
git clone https://github.com/BobLiterman/Rboretum.git
```
2) In R:  

```
library(devtools)
devtools::install(/path/to/git/clone)
```

### RBORETUM FUNCTIONS  

#### read.rooted()  
- **Function Description:**  
Reads unrooted trees (or rooted trees) into R as a phylo object rooted at a specified outgroup  

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

#### get.splits()
#### get.clades()

#### check.shared()
#### get.shared()
#### same.taxa()
#### same.topology()
#### check.comparable()
#### get.comparable()
#### compare.clades()

#### alignment.signal()
#### tree.support()
#### clade.support()

#### basic.treePlot()
#### support.treePlot()
#### pies.treePlot()
#### tandem.treePlot()

#### tableCount()
#### semiVector()


