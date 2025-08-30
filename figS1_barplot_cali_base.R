####################################
## figS1_barplot_cali_base.R
####################################

## JPJ 11 viii 25

## USAGE: Rscript figS1_barplot_cali_base.R
## PURPOSE: to create the base cali subset entropy barplot (Fig. S1)


## load data and packages

library(MetBrewer)
cali_cols <- met.brewer("Hokusai3", 5)

cali_ids <- read.delim("cali_sub_60.txt", header=F)
	dim(cali_ids)
	head(cali_ids)


## figure

cq2 <- read.csv("q2_cali_sub.txt", header=TRUE)
	dim(cq2)
	head(cq2)

cq3 <- read.csv("q3_cali_sub.txt", header=TRUE)
	dim(cq3)
	head(cq3)

cq4 <- read.csv("q4_cali_sub.txt", header=TRUE)
	dim(cq4)
	head(cq4)

cq5 <- read.csv("q5_cali_sub.txt", header=TRUE)
	dim(cq5)
	head(cq5)



clusters <- cbind(
					cq2[1:60, 2], cq2[61:120, 2],
					cq3[1:60, 2], cq3[61:120, 2], cq3[121:180, 2],
					cq4[1:60, 2], cq4[61:120, 2], cq4[121:180, 2], cq4[181:240, 2],
					cq5[1:60, 2], cq5[61:120, 2], cq5[121:180, 2], cq5[181:240, 2], cq5[241:300, 2]
					)
					
t_clusters <- t(clusters)

colors= c(
			cali_cols[5], cali_cols[2],
			cali_cols[2], cali_cols[5], cali_cols[3],
			cali_cols[5], cali_cols[2], cali_cols[4], cali_cols[3],
			cali_cols[3], cali_cols[1], cali_cols[5], cali_cols[4], cali_cols[2]
			)

pdf(file="figS1_barplot_cali_base.pdf", width=12, height=5)			
#quartz(width=12, height=5)
par(mar=c(4,0,0.2,0))
barplot(t_clusters, col=colors, beside=F, names.arg=cali_ids[,1], cex.names=.5, las=2, yaxt="n",  space=0,lwd=1.5)
dev.off()





