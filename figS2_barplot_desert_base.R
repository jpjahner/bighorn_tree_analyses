####################################
## figS1_barplot_desert_base.R
####################################

## JPJ 11 viii 25

## USAGE: Rscript figS1_barplot_desert_base.R
## PURPOSE: to create the base desert subset entropy barplot (Fig. S2)


## load data and packages

library(MetBrewer)
des_cols <- met.brewer("Greek", 5)

desert_ids <- read.delim("desert_sub_55.txt", header=F)
	dim(desert_ids)
	head(desert_ids)


## figure

dq2 <- read.csv("q2_desert_sub.txt", header=TRUE)
	dim(dq2)
	head(dq2)

dq3 <- read.csv("q3_desert_sub.txt", header=TRUE)
	dim(dq3)
	head(dq3)

dq4 <- read.csv("q4_desert_sub.txt", header=TRUE)
	dim(dq4)
	head(dq4)

dq5 <- read.csv("q5_desert_sub.txt", header=TRUE)
	dim(dq5)
	head(dq5)



clusters <- cbind(
					dq2[1:55, 2], dq2[56:110, 2],
					dq3[1:55, 2], dq3[56:110, 2], dq3[111:165, 2],
					dq4[1:55, 2], dq4[56:110, 2], dq4[111:165, 2], dq4[166:220, 2],
					dq5[1:55, 2], dq5[56:110, 2], dq5[111:165, 2], dq5[166:220, 2], dq5[221:275, 2]
					)
					
t_clusters <- t(clusters)

colors= c(
			des_cols[3], des_cols[2],
			des_cols[2], des_cols[3], des_cols[1],
			des_cols[1], des_cols[3], des_cols[4], des_cols[2],
			des_cols[4], des_cols[5], des_cols[3], des_cols[2], des_cols[1]
			)

pdf(file="figS2_barplot_desert_base.pdf", width=12, height=5)			
#quartz(width=12, height=5)
par(mar=c(4,0,0.2,0))
barplot(t_clusters, col=colors, beside=F, names.arg=desert_ids[,1], cex.names=.5, las=2, yaxt="n",  space=0,lwd=1.5)
dev.off()




