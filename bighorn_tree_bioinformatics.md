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
```

### calculate coverage

```{bash}
sbatch slurm_wyo_bighorn_calc_cov.sh
```










