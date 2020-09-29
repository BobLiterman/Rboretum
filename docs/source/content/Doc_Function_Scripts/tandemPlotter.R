library(Rboretum)

plot_1 <- ggtree(ape::rtree(10),branch.length = 'none') + geom_tiplab(offset=.2)
plot_2 <- ggtree(ape::rtree(10),branch.length = 'none') + geom_tiplab(offset=-.4) + scale_x_reverse()
plot_3 <- ggtree(ape::rtree(25),branch.length = 'none') + geom_tiplab(offset=.2)

tandemPlotter(plot_1,plot_2)
tandemPlotter(plot_1,plot_2,vertical = TRUE)
tandemPlotter(tandemPlotter(plot_1,plot_2),plot_3,vertical = TRUE)

