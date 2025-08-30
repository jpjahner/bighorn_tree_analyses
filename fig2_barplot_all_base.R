####################################################
## fig2_barplot_all_base.R
####################################################

## JPJ 8 viii 25

## USAGE: Rscript fig2_barplot_all_base.R
## PURPOSE: to make the base barplot figure of all individuals (Fig. 2)


## load data and packages

library(MetBrewer)
bpcols <- met.brewer("Archambault", 7)

ids <- read.delim("bighorn_tree_10inds_190.txt", header=TRUE)
	dim(ids)
	head(ids)

barplot_ids <- read.delim("bighorn_tree_10inds_190_barplot.txt", header=TRUE)
	dim(barplot_ids)
	head(barplot_ids)

q2 <- read.csv("q2_bighorn_tree.txt", header=TRUE)
	dim(q2)
	head(q2)

q3 <- read.csv("q3_bighorn_tree.txt", header=TRUE)
	dim(q3)
	head(q3)

q4 <- read.csv("q4_bighorn_tree.txt", header=TRUE)
	dim(q4)
	head(q4)

q5 <- read.csv("q5_bighorn_tree.txt", header=TRUE)
	dim(q5)
	head(q5)

q6 <- read.csv("q6_bighorn_tree.txt", header=TRUE)
	dim(q6)
	head(q6)

q7 <- read.csv("q7_bighorn_tree.txt", header=TRUE)
	dim(q7)
	head(q7)

##  q8 <- read.csv("q8_bighorn_tree.txt", header=TRUE)
##  	dim(q8)
##  	head(q8)

##  q9 <- read.csv("q9_bighorn_tree.txt", header=TRUE)
##   	dim(q9)
##  	head(q9)
	
##  q10 <- read.csv("q10_bighorn_tree.txt", header=TRUE)
##  	dim(q10)
##  	head(q10)




## make figure

clusters <- cbind(
					q2[1:190, 2], q2[191:380, 2],
					q3[1:190, 2], q3[191:380, 2], q3[381:570, 2],
					q4[1:190, 2], q4[191:380, 2], q4[381:570, 2], q4[571:760, 2],
					q5[1:190, 2], q5[191:380, 2], q5[381:570, 2], q5[571:760, 2], q5[761:950, 2],
					q6[1:190, 2], q6[191:380, 2], q6[381:570, 2], q6[571:760, 2], q6[761:950, 2], q6[951:1140, 2],
					q7[1:190, 2], q7[191:380, 2], q7[381:570, 2], q7[571:760, 2], q7[761:950, 2], q7[951:1140, 2], q7[1141:1330, 2]
					)

clusters_barplot <- matrix(NA, dim(clusters)[1], dim(clusters)[2])	
for (i in 1:dim(clusters_barplot)[1])
	{
	clusters_barplot[i,] <- clusters[ids[,1]==barplot_ids[i,1],]
	}

t_clusters <- t(clusters_barplot)

colors= c(
			bpcols[3], bpcols[1],
			bpcols[3], bpcols[5], bpcols[1],
			bpcols[5], bpcols[2], bpcols[1], bpcols[3],
			bpcols[3], bpcols[2], bpcols[5], bpcols[7], bpcols[1],
			bpcols[3], bpcols[5], bpcols[2], bpcols[4], bpcols[1], bpcols[7],
			bpcols[2], bpcols[7], bpcols[4], bpcols[1], bpcols[5], bpcols[3], bpcols[6]
			)

#quartz(width=10, height=6)
pdf(file="fig2_barplot_all_base.pdf", width=10, height=6)
par(mar=c(4,0,0.2,0))
barplot(t_clusters, col=colors, beside=F, names.arg=barplot_ids[,1], cex.names=0.2, las=2, yaxt="n",  space=0, lwd=0.75)
dev.off()




