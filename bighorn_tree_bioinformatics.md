# bighorn tree bioinformatics

NOTES: demultiplexing happened previously based on methods of Jahner et al. (2019) Evolutionary Applications. Starting with alignment of fastqs


GRABBING AN INTERACTIVE JOB:
```{bash}
salloc --account=evolgen --time=3:00:00
```


## populations

* California bighorn
    * Gilpin, British Columbia (CB_GI)
    * Hart Mountain, Oregon (CB_HMP)
    * Junction / Deer Park, British Columbia (CB_JD)
    * Kamloops Lake, British Columbia (CB_KL)
    * Sinlahekin, Washington (CB_SI)
    * South Okanagan, British Columbia (CB_SO)
* Desert bighorn (Mexican 1, Nelson 3, Peninsular 1)
    * Carrizo Canyon, California (PB_CC)
    * Kofa, Arizona (MB_KO)
    * Lone Mountain, Nevada (DB_212)
    * Muddy Mountains, Nevada (DB_268)
    * White Mountains, California (DB_WM)
* Rocky Mountain bighorn
    * Antelope Island, Utah (RB_AI)
    * Cadomin, Alberta (RB_HI)
    * Georgetown, Colorado (RB_GT)
    * Lower Salmon, Idaho (RB_LS)
    * Ram Mountain, Alberta (RB_RM)
    * Sheep River, Alberta (RB_SR)
    * South Salmo, British Columbia (RB_SS)
    * Whiskey Mountain, Wyoming (RB_WY)
* Sierra Nevada bighorn
    * Langley Mountain, California (SN_LM)


## alignment (using new bighorn genome)

### grabbing reference from NCBI
```{bash}
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/042/477/335/GCF_042477335.2_ARS-UI_OviCan_v2/*.fna.gz .

gunzip GCF_042477335.2_ARS-UI_OviCan_v2_genomic.fna.gz

grep -c "^>" GCF_042477335.2_ARS-UI_OviCan_v2_genomic.fna 
    ## 42
```



### remove 60 line endings
```{bash}
perl /home/jjahner/perl_scripts/remove60_fasta.pl GCF_042477335.2_ARS-UI_OviCan_v2_genomic.fna
```

### change reference names to something shorter

```{bash}
perl /home/jjahner/perl_scripts/rename_scaff.pl no60_GCF_042477335.2_ARS-UI_OviCan_v2_genomic.fna

mv renamed_no60_GCF_042477335.2_ARS-UI_OviCan_v2_genomic.fna.txt bh2_genome.fna
```


### make index for bwa

```{bash}
sbatch slurm_bh2tree_bwa_index.sh
```


### align 

```{bash}
sbatch slurm_bh2tree_runbwa.sh
```


### convert sam to bam

```{bash}
sbatch slurm_bh2tree_sam2bam.sh
```



### renaming south salmo files
Originally labeled as California bighorn, but actually Rocky Mountain
```{bash}
rename CB_SS RB_SS aln_CB_SS*
```



### make bam list

```{bash}
ls *.sorted.bam > bam_list.txt
```



## calling variants

```{bash}
sbatch slurm_variants.sh
```



## filtering

Number of loci in vcf
```{bash}
grep "^scaffold" -c variants_rawfiltered_27jun25.vcf
    ## 745105
```


making id file for reheadering:

```{bash}
ls *.sorted.bam > bams.txt

sed -s "s/aln_//" bams.txt | sed -s "s/.sorted.bam//" > bighorn_tree_ids_col.txt
```

reheader
```{bash}
module load arcc/1.0 gcc/14.2.0 bcftools/1.20 vcftools/0.1.17

bcftools reheader -s bighorn_tree_ids_col.txt variants_rawfiltered_27jun25.vcf -o rehead_variants_rawfiltered_27jun25.vcf
```

### first filter investigation
```{bash}
perl run_first_filter.pl rehead_variants_rawfiltered_27jun25.vcf

grep "Sites" first_filter_out/*
    ## first_filter_out/stdout_maf1_miss6:After filtering, kept 75664 out of a possible 745105 Sites
    ## irst_filter_out/stdout_maf1_miss7:After filtering, kept 46926 out of a possible 745105 Sites
    ## first_filter_out/stdout_maf1_miss8:After filtering, kept 22232 out of a possible 745105 Sites
    ## first_filter_out/stdout_maf1_miss9:After filtering, kept 7981 out of a possible 745105 Sites
    ## first_filter_out/stdout_maf2_miss6:After filtering, kept 62675 out of a possible 745105 Sites
    ## first_filter_out/stdout_maf2_miss7:After filtering, kept 38644 out of a possible 745105 Sites
    ## first_filter_out/stdout_maf2_miss8:After filtering, kept 18639 out of a possible 745105 Sites
    ## first_filter_out/stdout_maf2_miss9:After filtering, kept 6649 out of a possible 745105 Sites
    ## first_filter_out/stdout_maf3_miss6:After filtering, kept 53991 out of a possible 745105 Sites
    ## first_filter_out/stdout_maf3_miss7:After filtering, kept 33360 out of a possible 745105 Sites
    ## first_filter_out/stdout_maf3_miss8:After filtering, kept 16244 out of a possible 745105 Sites
    ## first_filter_out/stdout_maf3_miss9:After filtering, kept 5823 out of a possible 745105 Sites
    ## first_filter_out/stdout_maf4_miss6:After filtering, kept 47528 out of a possible 745105 Sites
    ## first_filter_out/stdout_maf4_miss7:After filtering, kept 29521 out of a possible 745105 Sites
    ## first_filter_out/stdout_maf4_miss8:After filtering, kept 14410 out of a possible 745105 Sites
    ## first_filter_out/stdout_maf4_miss9:After filtering, kept 5200 out of a possible 745105 Sites
    ## first_filter_out/stdout_maf5_miss6:After filtering, kept 42748 out of a possible 745105 Sites
    ## first_filter_out/stdout_maf5_miss7:After filtering, kept 26563 out of a possible 745105 Sites
    ## first_filter_out/stdout_maf5_miss8:After filtering, kept 12978 out of a possible 745105 Sites
    ## first_filter_out/stdout_maf5_miss9:After filtering, kept 4697 out of a possible 745105 Sites
```
Move forward with maf3 miss7


### create ids file
```{bash}
vcftools --vcf variants_maf3_miss7.recode.vcf --missing-indv

cut -f 1 out.imiss > bighorn_tree_10inds_190.txt

sed "s/INDV/ind/" bighorn_tree_10inds_190.txt > bighorn_tree_10inds_190_good_head.txt

grep -v "ind" bighorn_tree_10inds_190_good_head.txt | cut -f 1,2 -d "_" > bighorn_tree_10inds_pops_190.txt

```


### rename pops file

```{bash}
sed 's/CB_GI/GI/g' bighorn_tree_10inds_pops_190.txt |\
sed 's/CB_HMP/HM/g' |\
sed 's/CB_JD/JD/g' |\
sed 's/CB_KL/KL/g' |\
sed 's/CB_SI/SI/g' |\
sed 's/CB_SO/SO/g' |\
sed 's/DB_212/LM/g' |\
sed 's/DB_268/MM/g' |\
sed 's/DB_WM/WM/g' |\
sed 's/MB_KO/KO/g' |\
sed 's/PB_CC/CC/g' |\
sed 's/RB_AI/AI/g' |\
sed 's/RB_GT/GT/g' |\
sed 's/RB_HI/CA/g' |\
sed 's/RB_LS/LS/g' |\
sed 's/RB_RM/RM/g' |\
sed 's/RB_SR/SR/g' |\
sed 's/RB_SS/SS/g' |\
sed 's/RB_WY/WK/g' |\
sed 's/SN_LM/ML/g' > bighorn_tree_10inds_pops_190_renamed.txt
```



### generate mpgl

```{bash}
perl /home/jjahner/perl_scripts/vcf2mpglV1.3TLP.pl variants_maf3_miss7.recode.vcf
    ## Number of loci: 33084; number of individuals 190
```

### calculate coverage

```{bash}
sbatch slurm_bh2tree_calc_cov.sh
```


### filter out over-assembled loci

```{R}
cut -d " " -f 1 variants_maf3_miss7.recode.mpgl > loc_names_33084.txt

module load arcc/1.0 gcc/14.2.0 r/4.4.0
R

dat <- read.csv("variants_maf3_miss7.recode.mpgl_coverage.csv", header=F)
	dim(dat)
	dat[1:10,1:10]

dat_noname <- dat[,-1]
	dim(dat_noname)
	dat_noname[1:10,1:10]

loc_names <- read.delim("loc_names_33084.txt", header=F)
	head(loc_names)

avg_33084 <- vector()
in_out_33084_10 <- vector()
in_out_33084_12 <- vector()
in_out_33084_14 <- vector()
in_out_33084_16 <- vector()

for (i in 1:33084)
	{
	avg <- mean(dat_noname[,i])
	avg_33084 <- append(avg_33084, avg)
	
	if (avg <= 10)	{ in_out_33084_10 <- append(in_out_33084_10, 1) }
	else			{ in_out_33084_10 <- append(in_out_33084_10, 0) }

	if (avg <= 12)	{ in_out_33084_12 <- append(in_out_33084_12, 1) }
	else			{ in_out_33084_12 <- append(in_out_33084_12, 0) }
	
	if (avg <= 14)	{ in_out_33084_14 <- append(in_out_33084_14, 1) }
	else			{ in_out_33084_14 <- append(in_out_33084_14, 0) }
	
	if (avg <= 16)	{ in_out_33084_16 <- append(in_out_33084_16, 1) }
	else			{ in_out_33084_16 <- append(in_out_33084_16, 0) }
	
	}



sum(in_out_33084_10) / 33084
	## 31172 (94.2%)

sum(in_out_33084_12) / 33084
	## 31920 (96.5%)
	
sum(in_out_33084_14) / 33084
	## 32320 (97.7%)
	
sum(in_out_33084_16) / 33084
	## 32523 (98.3%)

 choosing to kill all locs with mean cov/ind >= 14



sub_14 <- dat_noname[,in_out_33084_14==1]
	dim(sub_14)
sub_14_avg <- subset(avg_33084, in_out_33084_14==1)	

kill_locs <- subset(loc_names, in_out_33084_14==0)
	dim(kill_locs)
	head(kill_locs)

write.table(kill_locs, file="high_cov_loc_list_to_be_removed.txt", quote=F, row.names=F, col.names=F)

quit()
```

```{bash}
sed "s/:/\t/" high_cov_loc_list_to_be_removed.txt > high_cov_loc_list_to_be_removed_tabdelim.txt

module load arcc/1.0 gcc/14.2.0 bcftools/1.20 vcftools/0.1.17

vcftools --vcf variants_maf3_miss7.recode.vcf --exclude-positions high_cov_loc_list_to_be_removed_tabdelim.txt --recode --recode-INFO-all --out variants_maf3_miss7_noHighCov
	## After filtering, kept 190 out of 190 Individuals
	## After filtering, kept 32596 out of a possible 33360 Sites
```



### create mpgl file

```{bash}
perl /home/jjahner/perl_scripts/vcf2mpglV1.3TLP.pl variants_maf3_miss7_noHighCov.recode.vcf
```



### calculate coverage

```{bash}
sbatch slurm_bh2tree_calc_cov.sh
```

IN R
```{bash}
module load arcc/1.0 gcc/14.2.0 r/4.4.0
R
j <- read.csv("variants_maf3_miss7_noHighCov.recode.mpgl_coverage.csv", header=F)
k <- j[,-1]
mean_vect <- vector()
for (i in 1:190) { mean_vect <- append(mean_vect, mean(as.numeric(k[i,]))) }
mean(mean_vect)
	## 4.146418
```




### generate pntest file

```{bash}
perl /home/jjahner/perl_scripts/gl2genestV1.3.pl variants_maf3_miss7_noHighCov.recode.mpgl mean
```


### make a transposed genotype matrix

```{bash}
module load arcc/1.0 gcc/14.2.0 r/4.4.0
R

read.table("pntest_mean_variants_maf3_miss7_noHighCov.recode.txt", header=F)->gl
read.table("bighorn_tree_10inds_190.txt", header=T)->ids
read.table("bighorn_tree_10inds_pops_190.txt", header=F)->pops

t(gl)->tgl
cbind(ids, pops, tgl)->tidsgl
write.table(tidsgl, file="bighorn_tree_gl_matrix_190_32320.txt", sep=" ", row.names=F, col.names=F , quote=F)
```





## entropy full

### LDA for starting values

```{R}
module load arcc/1.0 gcc/14.2.0 r/4.4.0
R

g <- read.table("pntest_mean_variants_maf3_miss7_noHighCov.recode.txt", header=F)

# dim(g)
# 32320  190

names <- read.table("bighorn_tree_10inds_190.txt", header=T)
pops <- read.table("bighorn_tree_10inds_pops_190.txt", header=F)

nind <- dim(g)[2]
nloci <- dim(g)[1]

## pca on the genotype matrix
#pcg<-prcomp(x=t(g),center=TRUE,scale=FALSE)

## OR calculate N x N genotype covariance matrix

gmn<-apply(g,1,mean, na.rm=T)
gmnmat<-matrix(gmn,nrow=nloci,ncol=nind)
gprime<-g-gmnmat ## remove mean
gcovarmat<-matrix(NA,nrow=nind,ncol=nind)
for(i in 1:nind){
    for(j in i:nind){
        if (i==j){
            gcovarmat[i,j]<-cov(gprime[,i],gprime[,j], use="pairwise.complete.obs")
        }
        else{
            gcovarmat[i,j]<-cov(gprime[,i],gprime[,j], use="pairwise.complete.obs")
            gcovarmat[j,i]<-gcovarmat[i,j]
        }
    }
}


## pca on the genotype covariance matrix
pcgcov<-prcomp(x=gcovarmat,center=TRUE,scale=FALSE)
pcgcov->pcg


## LDA

library(MASS)
k1<-kmeans(pcg$x[,1:5],1,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k2<-kmeans(pcg$x[,1:5],2,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k3<-kmeans(pcg$x[,1:5],3,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k4<-kmeans(pcg$x[,1:5],4,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k5<-kmeans(pcg$x[,1:5],5,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k6<-kmeans(pcg$x[,1:5],6,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k7<-kmeans(pcg$x[,1:5],7,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k8<-kmeans(pcg$x[,1:5],8,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k9<-kmeans(pcg$x[,1:5],9,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k10<-kmeans(pcg$x[,1:5],10,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k11<-kmeans(pcg$x[,1:5],11,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k12<-kmeans(pcg$x[,1:5],12,iter.max=10,nstart=10,algorithm="Hartigan-Wong")

ldak1<-lda(x=pcg$x[,1:5],grouping=k1$cluster,CV=TRUE)
ldak2<-lda(x=pcg$x[,1:5],grouping=k2$cluster,CV=TRUE)
ldak3<-lda(x=pcg$x[,1:5],grouping=k3$cluster,CV=TRUE)
ldak4<-lda(x=pcg$x[,1:5],grouping=k4$cluster,CV=TRUE)
ldak5<-lda(x=pcg$x[,1:5],grouping=k5$cluster,CV=TRUE)
ldak6<-lda(x=pcg$x[,1:5],grouping=k6$cluster,CV=TRUE)
ldak7<-lda(x=pcg$x[,1:5],grouping=k7$cluster,CV=TRUE)
ldak8<-lda(x=pcg$x[,1:5],grouping=k8$cluster,CV=TRUE)
ldak9<-lda(x=pcg$x[,1:5],grouping=k9$cluster,CV=TRUE)
ldak10<-lda(x=pcg$x[,1:5],grouping=k10$cluster,CV=TRUE)
ldak11<-lda(x=pcg$x[,1:5],grouping=k11$cluster,CV=TRUE)
ldak12<-lda(x=pcg$x[,1:5],grouping=k12$cluster,CV=TRUE)


write.table(round(ldak1$posterior,5),file="ldak1.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak2$posterior,5),file="ldak2.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak3$posterior,5),file="ldak3.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak4$posterior,5),file="ldak4.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak5$posterior,5),file="ldak5.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak6$posterior,5),file="ldak6.txt",quote=F,row.names=F,col.names=F)                                      
write.table(round(ldak7$posterior,5),file="ldak7.txt",quote=F,row.names=F,col.names=F)                                      
write.table(round(ldak8$posterior,5),file="ldak8.txt",quote=F,row.names=F,col.names=F)                                      
write.table(round(ldak9$posterior,5),file="ldak9.txt",quote=F,row.names=F,col.names=F)                                      
write.table(round(ldak10$posterior,5),file="ldak10.txt",quote=F,row.names=F,col.names=F)                                      
write.table(round(ldak11$posterior,5),file="ldak11.txt",quote=F,row.names=F,col.names=F)                                      
write.table(round(ldak12$posterior,5),file="ldak12.txt",quote=F,row.names=F,col.names=F)                                      
                                 
```



### Making .mpgl files for entropy

```{bash}
grep "_" bighorn_tree_10inds_190.txt > bighorn_tree_10inds_190_nohead.txt

perl /home/jjahner/perl_scripts/create_entropy_top_2rows.pl bighorn_tree_10inds_190_nohead.txt 

cat entropy_2rows.txt ../variants_maf3_miss7_noHighCov.recode.mpgl > bighorn_tree_entropy.mpgl
```

NOTE: need to manually add to the top 190 32320 1 


### launch

```{bash}
perl run_entropy.pl bighorn_tree_entropy.mpgl

    ## slurm-16859220 - slurm-16859269
    ## my $n_reps = 5;      ## number of replicate chains
    ## my $max_k = 10;       ## maximum k you want to consider
    ## my $ent_ploidy = 2;  ## proposed ploidy
    ## my $ent_l = 100000;   ## length of chain
    ## my $ent_b = 50000;    ## number of iterations to discard for burn-in
    ## my $ent_t = 10;      ## thinning interval
    ## my $ent_s = 20;      ## Dirichlet initialization scalar
    ## my $ent_e = 0.01;    ## per-locus error rate
    ## my $ent_m = 1;       ## input format (0 = read counts; 1 = genotype likelihoods)
    ## my $ent_w = 1;       ## output includes allele frequencies? [0/1]
    ## my $ent_D = 0;       ## dic or waic [0/1]
    ## my @lazy_caterer = ('1','2','7','11','16','22','29','37','46','56'); ## seeding for replicates (-r)
```










