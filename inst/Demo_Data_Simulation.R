library(Rboretum)

mammal_tree <- Rboretum::read.rooted(tree_path='RAxML_bipartitions.HomSap_SISRS_Biallelic_NoMissing_CDS',root_taxa = c("BalMys","BisBis","BosTau","BubBub","CapAeg","CapHir","ElaDav","GirTip","HipAmp","OdoVir","OkaJoh","OviAri"))

##### DATA SIMULATION #####
rates <- discrete.gamma(1,7)

align1 <- simSeq(mammal_tree, l = 5000, type="USER", levels=c('A','C','T','G','N'),bf = c(.2375,.2375,.2375,.2375,.05), rate=rates[1])
align2 <- simSeq(mammal_tree, l = 5000, type="USER", levels=c('a','c','t','g','N'),bf = c(.2125,.2125,.2125,.2125,.15), rate=rates[3])
align3 <- simSeq(mammal_tree, l = 5000, type="USER", levels=c('a','c','t','g','-','N'),bf = c(.2125,.2125,.2125,.2125,.10,.05), rate=rates[5])
align4 <- simSeq(mammal_tree, l = 5000, type="USER", levels=c('a','c','t','g','N'),bf=c(.2375,.2375,.2375,.2375,.05), rate=rates[7])

write.phyDat(align1,format = "phylip",file = 'Alignment_1.phylip')
write.phyDat(align2,format = "fasta",file = 'Alignment_2.fasta')
write.phyDat(align3,format = "nexus",file = 'Alignment_3.nexus')
write.phyDat(align4,format = "phylip",file = 'Alignment_4.phylip')
