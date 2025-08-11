####################################
## figS1_barplot_rocky_base.R
####################################

## JPJ 11 viii 25

## USAGE: Rscript figS1_barplot_rocky_base.R
## PURPOSE: to create the base rocky subset entropy barplot (Fig. S3)


## load data and packages

library(MetBrewer)
rocky_cols <- met.brewer("VanGogh1", 5)

rocky_ids <- read.delim("rocky_sub_75.txt", header=F)
	dim(rocky_ids)
	head(rocky_ids)


## figure

rq2 <- read.csv("q2_rocky_sub.txt", header=TRUE)
	dim(rq2)
	head(rq2)

rq3 <- read.csv("q3_rocky_sub.txt", header=TRUE)
	dim(rq3)
	head(rq3)

rq4 <- read.csv("q4_rocky_sub.txt", header=TRUE)
	dim(rq4)
	head(rq4)

rq5 <- read.csv("q5_rocky_sub.txt", header=TRUE)
	dim(rq5)
	head(rq5)



clusters <- cbind(
					rq2[1:75, 2], rq2[76:150, 2],
					rq3[1:75, 2], rq3[76:150, 2], rq3[151:225, 2],
					rq4[1:75, 2], rq4[76:150, 2], rq4[151:225, 2], rq4[226:300, 2],
					rq5[1:75, 2], rq5[76:150, 2], rq5[151:225, 2], rq5[226:300, 2], rq5[301:375, 2]
					)
					
t_clusters <- t(clusters)

colors= c(
			rocky_cols[1], rocky_cols[3],
			rocky_cols[3], rocky_cols[2], rocky_cols[1],
			rocky_cols[1], rocky_cols[4], rocky_cols[2], rocky_cols[3],
			rocky_cols[1], rocky_cols[5], rocky_cols[2], rocky_cols[3], rocky_cols[4]
			)

pdf(file="figS3_barplot_rocky_base.pdf", width=12, height=5)			
#quartz(width=12, height=5)
par(mar=c(4,0,0.2,0))
barplot(t_clusters, col=colors, beside=F, names.arg=rocky_ids[,1], cex.names=.4, las=2, yaxt="n",  space=0,lwd=1.5)
dev.off()




