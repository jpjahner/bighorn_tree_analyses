###########################################
## fig4_pca.R
###########################################

## JPJ 8 viii 25

## USAGE: Rscript fig4_pca.R
## PURPOSE: to make the pca figure (Fig. 4)


## load data and packages

library(MetBrewer)
bpcols <- met.brewer("Archambault", 7)

gprobs <- read.csv("gprob3_bighorn_tree.txt", header=TRUE)
	dim(gprobs)
	gprobs[1:10,1:10]

gprobs_noname <- gprobs[,-1]
	dim(gprobs_noname)
	gprobs_noname[1:10,1:10]

new_pops <- read.delim("bighorn_tree_10inds_pops_190_renamed.txt", header=FALSE)
	dim(new_pops)
	head(new_pops)
	
	
## analyses (full and subsets)

full_pca <- prcomp(gprobs_noname, center=TRUE, scale=FALSE)
	summary(full_pca)

	## Importance of components:
	##                          PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8
	## Standard deviation     26.28 20.8360 9.89855 9.83749 9.58758 9.02884 8.08547 7.62842
	## Proportion of Variance  0.19  0.1194 0.02696 0.02662 0.02529 0.02243 0.01799 0.01601
	## Cumulative Proportion   0.19  0.3094 0.33638 0.36300 0.38829 0.41072 0.42870 0.44471


rocky_pops <- as.matrix(new_pops[109:183,])
rocky_pca <- prcomp(gprobs_noname[109:183,], center=TRUE, scale=FALSE)
	summary(rocky_pca)

	## Importance of components:
	##                            PC1     PC2      PC3     PC4    PC5     PC6     PC7     PC8
	## Standard deviation     17.7715 14.5491 10.90282 9.78406 9.0973 8.38222 7.69246 6.68910
	## Proportion of Variance  0.1156  0.0775  0.04352 0.03505 0.0303 0.02572 0.02167 0.01638
	## Cumulative Proportion   0.1156  0.1931  0.23665 0.27170 0.3020 0.32773 0.34939 0.36578


desSie_pops <- as.matrix(new_pops[c(61:108,184:190),])
desSie_pca <- prcomp(gprobs_noname[c(61:108,184:190),], center=TRUE, scale=FALSE)
	summary(desSie_pca)

	## Importance of components:
	##                            PC1     PC2      PC3      PC4      PC5     PC6     PC7     PC8
	## Standard deviation     18.4728 18.2260 17.01314 15.05244 13.86591 8.17979 7.74796 7.50635
	## Proportion of Variance  0.1053  0.1025  0.08933  0.06992  0.05933 0.02065 0.01853 0.01739
	## Cumulative Proportion   0.1053  0.2078  0.29715  0.36708  0.42641 0.44706 0.46559 0.48297


cali_pops <- as.matrix(new_pops[1:60,])
cali_pca <- prcomp(gprobs_noname[1:60,], center=TRUE, scale=FALSE)
	summary(cali_pca)

	## Importance of components:
	##                             PC1      PC2     PC3     PC4     PC5     PC6     PC7     PC8
	## Standard deviation     13.80237 11.33508 8.80694 7.82713 7.33711 7.14695 6.78331 6.34065
	## Proportion of Variance  0.09753  0.06578 0.03971 0.03136 0.02756 0.02615 0.02356 0.02058
	## Cumulative Proportion   0.09753  0.16331 0.20302 0.23438 0.26194 0.28809 0.31165 0.33223



## figure

#quartz(width=8, height=8)
pdf(file="fig4_pca.pdf", width=8, height=8)
par(mar=c(5,5,2,1), mfrow=c(2,2))

## full pca
plot(full_pca$x[,1], full_pca$x[,2], type="n", xlab="PC 1 (19.0%)", ylab="PC 2 (11.9%)", cex.lab=1.5, cex.axis=1.25, las=1)
points(full_pca$x[new_pops[,1]=="JD", 1], full_pca$x[new_pops[,1]=="JD", 2], pch=21, cex=2, bg="#007FDA")
points(full_pca$x[new_pops[,1]=="KL", 1], full_pca$x[new_pops[,1]=="KL", 2], pch=22, cex=2, bg="#007FDA")
points(full_pca$x[new_pops[,1]=="SO", 1], full_pca$x[new_pops[,1]=="SO", 2], pch=23, cex=2, bg="#007FDA")
points(full_pca$x[new_pops[,1]=="GI", 1], full_pca$x[new_pops[,1]=="GI", 2], pch=24, cex=2, bg="#007FDA")
points(full_pca$x[new_pops[,1]=="SI", 1], full_pca$x[new_pops[,1]=="SI", 2], pch=25, cex=2, bg="#007FDA")
points(full_pca$x[new_pops[,1]=="HM", 1], full_pca$x[new_pops[,1]=="HM", 2], pch=21, cex=2, bg="grey60", col="#007FDA", lwd=1.5)
points(full_pca$x[new_pops[,1]=="WM", 1], full_pca$x[new_pops[,1]=="WM", 2], pch=23, cex=2, bg="#FF441E")
points(full_pca$x[new_pops[,1]=="LM", 1], full_pca$x[new_pops[,1]=="LM", 2], pch=21, cex=2, bg="#FF441E")
points(full_pca$x[new_pops[,1]=="MM", 1], full_pca$x[new_pops[,1]=="MM", 2], pch=22, cex=2, bg="#FF441E")
points(full_pca$x[new_pops[,1]=="KO", 1], full_pca$x[new_pops[,1]=="KO", 2], pch=24, cex=2, bg="#FF441E")
points(full_pca$x[new_pops[,1]=="CC", 1], full_pca$x[new_pops[,1]=="CC", 2], pch=25, cex=2, bg="#FF441E")
points(full_pca$x[new_pops[,1]=="ML", 1], full_pca$x[new_pops[,1]=="ML", 2], pch=21, cex=2, bg="#FFCB53")
points(full_pca$x[new_pops[,1]=="CA", 1], full_pca$x[new_pops[,1]=="CA", 2], pch=21, cex=2, bg="#B300BE")
points(full_pca$x[new_pops[,1]=="RM", 1], full_pca$x[new_pops[,1]=="RM", 2], pch=22, cex=2, bg="#B300BE")
points(full_pca$x[new_pops[,1]=="SR", 1], full_pca$x[new_pops[,1]=="SR", 2], pch=23, cex=2, bg="#B300BE")
points(full_pca$x[new_pops[,1]=="SS", 1], full_pca$x[new_pops[,1]=="SS", 2], pch=24, cex=2, bg="#B300BE")
points(full_pca$x[new_pops[,1]=="AI", 1], full_pca$x[new_pops[,1]=="AI", 2], pch=25, cex=2, bg="#B300BE")
points(full_pca$x[new_pops[,1]=="LS", 1], full_pca$x[new_pops[,1]=="LS", 2], pch=21, cex=2, bg="grey60", col="#B300BE", lwd=1.5)
points(full_pca$x[new_pops[,1]=="WK", 1], full_pca$x[new_pops[,1]=="WK", 2], pch=22, cex=2, bg="grey60", col="#B300BE", lwd=1.5)
points(full_pca$x[new_pops[,1]=="GT", 1], full_pca$x[new_pops[,1]=="GT", 2], pch=23, cex=2, bg="grey60", col="#B300BE", lwd=1.5)
mtext("A", adj=-0.25, cex=2)
mtext("All herds", cex=1.5, line=0.25)
box(lwd=1.25)


## rocky pca
plot(rocky_pca$x[,1], rocky_pca$x[,2], type="n", xlab="PC 1 (11.6%)", ylab="PC 2 (7.8%)", cex.lab=1.5, cex.axis=1.25, las=1)
points(rocky_pca$x[rocky_pops[,1]=="CA", 1], rocky_pca$x[rocky_pops[,1]=="CA", 2], pch=21, cex=2, bg="#B300BE")
points(rocky_pca$x[rocky_pops[,1]=="RM", 1], rocky_pca$x[rocky_pops[,1]=="RM", 2], pch=22, cex=2, bg="#B300BE")
points(rocky_pca$x[rocky_pops[,1]=="SR", 1], rocky_pca$x[rocky_pops[,1]=="SR", 2], pch=23, cex=2, bg="#B300BE")
points(rocky_pca$x[rocky_pops[,1]=="SS", 1], rocky_pca$x[rocky_pops[,1]=="SS", 2], pch=24, cex=2, bg="#B300BE")
points(rocky_pca$x[rocky_pops[,1]=="AI", 1], rocky_pca$x[rocky_pops[,1]=="AI", 2], pch=25, cex=2, bg="#B300BE")
points(rocky_pca$x[rocky_pops[,1]=="LS", 1], rocky_pca$x[rocky_pops[,1]=="LS", 2], pch=21, cex=2, bg="grey60", col="#B300BE", lwd=1.5)
points(rocky_pca$x[rocky_pops[,1]=="WK", 1], rocky_pca$x[rocky_pops[,1]=="WK", 2], pch=22, cex=2, bg="grey60", col="#B300BE", lwd=1.5)
points(rocky_pca$x[rocky_pops[,1]=="GT", 1], rocky_pca$x[rocky_pops[,1]=="GT", 2], pch=23, cex=2, bg="grey60", col="#B300BE", lwd=1.5)
mtext("B", adj=-0.25, cex=2)
mtext("Rocky herds", cex=1.5, line=0.25)
legend("topright", legend=c("CA", "RM", "SR", "SS", "AI", "LS", "WK", "GT"), pch=c(21,22,23,24,25,21,22,23), pt.bg=c(rep("#B300BE",5), rep("grey60",3)), col=c(rep("black",5), rep("#B300BE",3)), pt.cex=1.75, cex=1.1)
box(lwd=1.25)


## desert + sierra pca
plot(desSie_pca$x[,1], desSie_pca$x[,2], type="n", xlab="PC 1 (10.5%)", ylab="PC 2 (10.3%)", cex.lab=1.5, cex.axis=1.25, las=1)
points(desSie_pca$x[desSie_pops[,1]=="WM", 1], desSie_pca$x[desSie_pops[,1]=="WM", 2], pch=21, cex=2, bg="#FF441E")
points(desSie_pca$x[desSie_pops[,1]=="LM", 1], desSie_pca$x[desSie_pops[,1]=="LM", 2], pch=22, cex=2, bg="#FF441E")
points(desSie_pca$x[desSie_pops[,1]=="MM", 1], desSie_pca$x[desSie_pops[,1]=="MM", 2], pch=23, cex=2, bg="#FF441E")
points(desSie_pca$x[desSie_pops[,1]=="KO", 1], desSie_pca$x[desSie_pops[,1]=="KO", 2], pch=24, cex=2, bg="#FF441E")
points(desSie_pca$x[desSie_pops[,1]=="CC", 1], desSie_pca$x[desSie_pops[,1]=="CC", 2], pch=25, cex=2, bg="#FF441E")
points(desSie_pca$x[desSie_pops[,1]=="ML", 1], desSie_pca$x[desSie_pops[,1]=="ML", 2], pch=21, cex=2, bg="#FFCB53")
mtext("C", adj=-0.25, cex=2)
mtext("Desert + Sierra herds", cex=1.5, line=0.25)
legend("bottomleft", legend=c("WM", "LM", "MM", "KO", "CC", "ML"), pch=c(21,22,23,24,25,21), pt.bg=c(rep("#FF441E",5), "#FFCB53"), pt.cex=1.75, cex=1.1)
box(lwd=1.25)


## cali pca
plot(cali_pca$x[,1], cali_pca$x[,2], type="n", xlab="PC 1 (9.8%)", ylab="PC 2 (6.6%)", cex.lab=1.5, cex.axis=1.25, las=1)
points(cali_pca$x[cali_pops[,1]=="JD", 1], cali_pca$x[cali_pops[,1]=="JD", 2], pch=21, cex=2, bg="#007FDA")
points(cali_pca$x[cali_pops[,1]=="KL", 1], cali_pca$x[cali_pops[,1]=="KL", 2], pch=22, cex=2, bg="#007FDA")
points(cali_pca$x[cali_pops[,1]=="SO", 1], cali_pca$x[cali_pops[,1]=="SO", 2], pch=23, cex=2, bg="#007FDA")
points(cali_pca$x[cali_pops[,1]=="GI", 1], cali_pca$x[cali_pops[,1]=="GI", 2], pch=24, cex=2, bg="#007FDA")
points(cali_pca$x[cali_pops[,1]=="SI", 1], cali_pca$x[cali_pops[,1]=="SI", 2], pch=25, cex=2, bg="#007FDA")
points(cali_pca$x[cali_pops[,1]=="HM", 1], cali_pca$x[cali_pops[,1]=="HM", 2], pch=21, cex=2, bg="grey60", col="#007FDA", lwd=1.5)
mtext("D", adj=-0.25, cex=2)
mtext("California herds", cex=1.5, line=0.25)
legend("bottomright", legend=c("JD", "KL", "SO", "GI", "SI", "HM"), pch=c(21,22,23,24,25,21), pt.bg=c(rep("#007FDA",5), "grey60"), col=c(rep("black",5), "#007FDA"), pt.cex=1.75, cex=1.1)
box(lwd=1.25)

dev.off()












