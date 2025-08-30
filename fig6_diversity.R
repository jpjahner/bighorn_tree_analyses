###############################################
## fig6_diversity.R
###############################################

## JPJ 8 viii 25

## USAGE: Rscript fig6_diversity.pdf
## PURPOSE: to create the genetic diversity figure (Fig. 6)


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


## analyses

pop_list <- unique(new_pops[,1])
afreqs <- matrix(NA, length(pop_list), dim(gprobs_noname)[2])
for (i in 1:length(pop_list)){
	print(i)
	for (j in 1:dim(gprobs_noname)[2]){
		afreqs[i,j] <- mean(gprobs_noname[new_pops[,1]==pop_list[i],j]) / 2
	}
}

het <- 2 * afreqs * (1 - afreqs)
mean_het <- matrix(NA, length(pop_list), 1)
for (i in 1:length(pop_list)){
	mean_het[i,1] <- mean(het[i,])
}


## figure

pdf(file="fig6_diversity.pdf", height=5, width=6)
layout(matrix(c(1,1,1,1,1,2), nrow=1, ncol=6))
par(mar=c(5,5,1,1))
## het
plot(0, type="n", xlim=c(0,3), ylim=c(0.12,0.22), xlab="Lineage", ylab="Expected heterozygosity", cex.lab=2, cex.axis=1.5, xaxt="n")
axis(1, at=c(0.5, 1.5, 2.5), labels=c("California", "Desert + Sierra", "Rocky"), cex.axis=1.5)
points(0.5, mean_het[pop_list=="JD", 1], pch=21, cex=3, bg="#007FDA")
points(0.5, mean_het[pop_list=="KL", 1], pch=22, cex=3, bg="#007FDA")
points(0.5, mean_het[pop_list=="SO", 1], pch=23, cex=3, bg="#007FDA")
points(0.5, mean_het[pop_list=="GI", 1], pch=24, cex=3, bg="#007FDA")
points(0.5, mean_het[pop_list=="SI", 1], pch=25, cex=3, bg="#007FDA")
points(0.5, mean_het[pop_list=="HM", 1], pch=21, cex=3, bg="grey60", col="#007FDA", lwd=1.5)
points(1.5, mean_het[pop_list=="WM", 1], pch=23, cex=3, bg="#FF441E")
points(1.5, mean_het[pop_list=="LM", 1], pch=21, cex=3, bg="#FF441E")
points(1.5, mean_het[pop_list=="MM", 1], pch=22, cex=3, bg="#FF441E")
points(1.5, mean_het[pop_list=="KO", 1], pch=24, cex=3, bg="#FF441E")
points(1.5, mean_het[pop_list=="CC", 1], pch=25, cex=3, bg="#FF441E")
points(1.5, mean_het[pop_list=="ML", 1], pch=21, cex=3, bg="#FFCB53")
points(2.5, mean_het[pop_list=="CA", 1], pch=21, cex=3, bg="#B300BE")
points(2.5, mean_het[pop_list=="RM", 1], pch=22, cex=3, bg="#B300BE")
points(2.5, mean_het[pop_list=="SR", 1], pch=23, cex=3, bg="#B300BE")
points(2.5, mean_het[pop_list=="SS", 1], pch=24, cex=3, bg="#B300BE")
points(2.5, mean_het[pop_list=="AI", 1], pch=25, cex=3, bg="#B300BE")
points(2.5, mean_het[pop_list=="LS", 1], pch=21, cex=3, bg="grey60", col="#B300BE", lwd=1.5)
points(2.5, mean_het[pop_list=="WK", 1], pch=22, cex=3, bg="grey60", col="#B300BE", lwd=1.5)
points(2.5, mean_het[pop_list=="GT", 1], pch=23, cex=3, bg="grey60", col="#B300BE", lwd=1.5)

## legends
par(mar=c(0,0,0,0))
plot(0, type="n", xlab="", ylab="", axes=FALSE)
legend("top", box.lty=0, legend=c("CA", "RM", "SR", "SS", "AI", "LS", "WK", "GT", "ML", "WM", "LM", "MM", "KO", "CC", "JD", "KL", "SO", "GI", "SI", "HM"), pch=c(21,22,23,24,25,21,22,23,21,22,23,24,25,21,21,22,23,24,25,21), pt.bg=c(rep("#B300BE",5), rep("grey60",3), "#FFCB53", rep("#FF441E",5), rep("#007FDA",5), "grey60"), col=c(rep("black",5), rep("#B300BE",3), rep("black",11), "#007FDA"), pt.cex=2.5, cex=1.6)
dev.off()




