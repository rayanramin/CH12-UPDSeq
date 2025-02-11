### to perform the following steps you need these tools and R packages
# bwa 
# samtools
# picard
# macs2
# bam-readcount
# R packages: tidyverse, remotes, bioseq, RSamtools, GenomicRanges

### we have not tested this on other systems but you can try installing these tools using the conda/mamba package manager

#### creating the conda environment
conda create -n ch12_updseq
conda activate ch12_updseq
mamba install bam-readcount
mamba install r-base r-tidyverse r-stringr r-remotes r-bioseq 
mamba install bioconda::bioconductor-rsamtools bioconda::bioconductor-genomicranges bioconda::samtools 
R -e ".libPaths(.libPaths()[grep('.conda',.libPaths())]); remotes::install_github('fkeck/bioseq',upgrade='never')"

### alignment step
### Replace [the values in brackets] with the appropriate values for your system/environment
ref=[path/to/reference/genome(mm9)]
jardir=[path/to/picard.jar]
name=[Sample_Name]
R1=[path/to/R1.fastq/file]
R2=[path/to/R2.fastq/file]
tmpdir=[path/to/tmpdir]

outfiles_1=($name.mm9.aligned.1.bam)
outfiles_1b=($name.mm9.sorted.1b.bam)
outfiles_2a=($name.mm9.marked.2a.bam)
mfiles=($name.mm9.dup_metrics.txt)
outfiles_2b=($name.mm9.sorted.dedup.bam)
###
bwa mem -R '@RG\tID:1\tSM:'$name'\tPL:illumina\tLB:lib1\tPU:unit1' $ref $R1 $R2 | samtools view -bS - > $outfiles_1
java -Xmx32g -Djava.io.tmpdir=$tmpdir -jar $jardir'picard.jar' SortSam I=$outfiles_1 O=$outfiles_1b SORT_ORDER=coordinate TMP_DIR=$tmpdir
java -Djava.io.tmpdir=$tmpdir -Xmx32g -jar $jardir'picard.jar' MarkDuplicates I=$outfiles_1b O=$outfiles_2a CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$mfiles TMP_DIR=$tmpdir
samtools view -bF 0X400  $outfiles_1b > $outfiles_2b
###
# QC filter step (see reference 82)
###
# GCmm filter step
input=$outfiles_2a
output_3=($name.mm9.marked.3.GCmm.bam)
filter_reads_GCmm.sh -i $input -o $output_3
###
# MACS2 peak calling step
macs2 callpeak \
-t [treatment/stimulated Sample bam file] \
-c [control/unstimulated Sample bam file] \
-n [output file name] \
--bw 700 \
--keep-dup all
###


### ### ### 
### UI and Double Mismatch analysis:
## The analysis can be performed using the following R scripts, the required inputs are:
## 1- the input bam files in a list, txt file with 1 bam file per line
## 2- bed file containing the coordinates of the regions of interest
## 3- reference genome fasta file
## optional inputs:
## 4- snp database file (i.e. mm9 snp128.txt)
## 5- prefix for the output files (default: output)
## 6- output directory (needs to be created before running the script) (default: current directory)
## 7-contexts to analyze (default: NC, WRCY & nWRCY)


#########
# ## To download mm9 snps
# # wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/database/snp128.txt.gz
# # gunzip snp128.txt.gz
snpdb='path/to/mm9/snp128.txt'
#########

## bed file containing the coordinates of the peaks used in the manuscript
# cat macs2_peaks_switch_regions.bed

# chr12   114663630   114664335   Set1_Sm
# chr12   114499853   114500188   Set1_Sa1
# chr12   114501040   114501639   Set1_Sa2
# chr12   114663520   114664931   Set2_Sm
# chr12   114499843   114500269   Set2_Sa
# chr12   114663650   114664705   Set3_Sm

## bed file containing the coordinates of the peaks
# cat Peaks_IGH.bed

# chr12   114663520       114664931       Sm
# chr12   114499843       114500269       Sa1
# chr12   114501040       114501639       Sa2
# chr12   114656721       114660941       Cm
# chr12   114492992       114498447       Ca


####
## put the path to the bam files in a text file
ref=/path/to/ref/genome

Rscript ./RunUI_Analysis.R \
 --regions=macs2_peaks_switch_regions.bed \
 --bams=bamlist.txt \
 --ref=$ref \
 --snpdb=$snpdb \
 --prefix=outputs1

Rscript ./MisMatchExtractor.R \
 --regions=Peaks_IGH.bed \
 --bams=bamlist.txt \
 --ref=$ref \
 --snpdb=$snpdb \
 --prefix=outputs2
