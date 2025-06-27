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





