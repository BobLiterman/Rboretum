## Rboretum: Comparisons and Visualizations of Phylogenetic Tree Topology and Alignment Support in R  
**Dr. Bob Literman**  
**Dr. Rachel Schwartz**  
University of Rhode Island  

v.1.0 (01/01/2020)  
** **
### Program Dependencies
- R version >= 3.6.1
- Python version >= 3.6+ w/[Biopython](https://biopython.org/wiki/Download)  
** **
### R Package Dependencies  
- ape (>= 5.4.1)  
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
- ggimage
- grDevices (>= 3.6.1)  
- viridisLite (>= 0.3.0)  
- rlist (>= 0.4.6.1)  
- BiocManager (>= 1.30.10)  
** **
### Rboretum Installation
```
library(devtools)
devtools::install_github("BobLiterman/Rboretum",build_vignettes=TRUE)
```

Note that if you have any problems relating to the installation of `reticulate` see 
the following [issue](https://github.com/BobLiterman/Rboretum/issues/10).

** **
### Loading Library
```
library(Rboretum)
sourceRboretum()
```
** **
Rboretum is very much a work in progress, **and raising early issues would be very helpful**. What do you wish you could easily do in R? Let me know and I'll try to make it happen for you.

For basic usage, please:  
- See the associated vignette.
  ```
  vignette(package="Rboretum", "vignette")
  ```
- Visit the ReadTheDocs (under construction)

Feel free to reach out to me on Twitter (@BobLiterman) or over e-mail (Robert.Literman@fda.hhs.gov)
** **