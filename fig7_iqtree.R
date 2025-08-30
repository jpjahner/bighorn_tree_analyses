####################################################
## fig7_iqtree.R
####################################################

## JPJ 17 viii 25

## USAGE: Rscript fig7_iqtree.R
## PURPOSE: to make the phylogenetic tree based on the iqtree analysis


## load data and packages

library(ape)
library(ggplot2)
library(ggtree)


## figure

tree <- read.tree("bh2tree.contree")
tree_noroot <- drop.tip(tree, "ST_ST_1500")

tree_gg <- ggtree(tree_noroot, layout="rectangular", branch.length="none", right=TRUE) +
				geom_tiplab(size=3, hjust=-0.1) +
				geom_nodelab(size=2.5, hjust=-0.2, color="black") +
				theme(legend.position = "none")

lins <- c("Stone (outgroup)", rep("California", 12), rep("Desert", 10), rep("Sierra", 2), rep("Rocky", 16))
lins_frame <- data.frame(lins)
rownames(lins_frame) <- tree$tip.label


phylo_lins <- gheatmap(tree_gg, lins_frame, offset=7, width=0.15, colnames=F) +
					scale_fill_manual(values=c("#007FDA", "#FF441E", "#B300BE", "#FFCB53", "black"), name="Lineage")

pdf(file="fig7_iqtree.pdf", height=10, width=6)
#quartz(height=10, width=6)
phylo_lins
dev.off()

