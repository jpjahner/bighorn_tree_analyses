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
