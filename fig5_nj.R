##################################################
## fig5_nj.R
####################################################

## JPJ 8 viii 25

## USAGE: Rscript fig5_nj.R
## PURPOSE: to make the neighbor joining tree (Fig. 5) and export pairwise Fsts (Table S4)


## load data and packages

library(ape)

gprobs <- read.csv("gprob3_bighorn_tree.txt", header=TRUE)
	dim(gprobs)
	gprobs[1:10,1:10]

gprobs_noname <- gprobs[,-1]
	dim(gprobs_noname)
	gprobs_noname[1:10,1:10]
	
new_pops <- read.delim("bighorn_tree_10inds_pops_190_renamed.txt", header=FALSE)
	dim(new_pops)
	head(new_pops)


## calculate pairwise fst

uniq_pops <- c("CA", "RM", "SR", "SS", "AI", "LS", "WK", "GT", "ML", "WM", "LM", "MM", "KO", "CC", "JD", "KL", "SO", "GI", "SI", "HM")
	## matches order for Table S4

afreqs <- matrix(0, length(uniq_pops), dim(gprobs_noname)[2])
for (i in 1:length(uniq_pops))
	{
	sub_dat <- subset(gprobs_noname, new_pops[,1]==as.character(uniq_pops[i]))
	for (j in 1:dim(gprobs_noname)[2])
		{
		af <- mean(sub_dat[,j]) / 2
		afreqs[i,j] <- af
		}
	print(i)
	}
afreqs[,1:10]


hudsonFst2<-function(p1=NA, p2=NA, n1=NA, n2=NA)
	{
    numerator<- p1 * (1 - p1) + p2 * (1 - p2)
    denominator<- p1 * (1 - p2) + p2 * (1 - p1)
    fst<- 1 - numerator/denominator
    out <- cbind(numerator, denominator, fst)
    return(out)
	}

fst_mean_out <- matrix(0, length(uniq_pops), length(uniq_pops))

for (i in 1:length(uniq_pops))
	{
	for (j in 1:length(uniq_pops))
		{
		if (i==j) { fst_mean_out[i,j] <- 0 }
		else if (i>j)
			{
			locus_vector <- vector()
			for (k in 1:dim(gprobs_noname)[2])
				{
				if 		(afreqs[i,k]==0 && afreqs[j,k]==0) { locus_vector <- append(locus_vector, 0) }
				else if (afreqs[i,k]==1 && afreqs[j,k]==1) { locus_vector <- append(locus_vector, 0) }
				else
					{
					loc_fst <- hudsonFst2(p1=afreqs[i,k], p2=afreqs[j,k])
					locus_vector <- append(locus_vector, loc_fst[3])
					}
				}
			mean_fst <- mean(locus_vector)
			fst_mean_out[i,j] <- mean_fst
			fst_mean_out[j,i] <- mean_fst
			}
		print(i); print(j)
		}
	}

rownames(fst_mean_out) <- uniq_pops
colnames(fst_mean_out) <- uniq_pops

write.table(fst_mean_out, file="fst_mean_out.txt", quote=F)


## nj tree

tr <- nj(as.dist(fst_mean_out))
#quartz(height=6,width=6)
pdf(file="fig5_nj.pdf", height=6, width=6)
par(mar=c(0,0,0,0))
plot(tr, type="unrooted", use.edge.length=TRUE, edge.width=3, label.offset=2)
tiplabels(pch=c(21,22,23,24,25,21,22,23,21,21,22,23,24,25,21,22,23,24,25,21), bg=c(rep("#B300BE",5), rep("grey60",3), "#FFCB53", rep("#FF441E",5), rep("#007FDA",5), "grey60"), col=c(rep("black", 5), rep("#B300BE",3), rep("black", 11), "#007FDA"), cex=2.5, lwd=2.5)
legend("bottomright", legend=c("CA", "RM", "SR", "SS", "AI", "LS", "WK", "GT"), pch=c(21,22,23,24,25,21,22,23), pt.bg=c(rep("#B300BE",5), rep("grey60",3)), col=c(rep("black",5), rep("#B300BE",3)), pt.cex=1.75, cex=1.1)
legend("topleft", legend=c("WM", "LM", "MM", "KO", "CC", "ML"), pch=c(21,22,23,24,25,21), pt.bg=c(rep("#FF441E",5), "#FFCB53"), pt.cex=1.75, cex=1.1)
legend("topright", legend=c("JD", "KL", "SO", "GI", "SI", "HM"), pch=c(21,22,23,24,25,21), pt.bg=c(rep("#007FDA",5), "grey60"), col=c(rep("black",5), "#007FDA"), pt.cex=1.75, cex=1.1)
dev.off()









