# Scripts to generate CN estimates for Multicopy Genes and VNTRs and run PheWAS analysis
This repository contains scripts for generating Copy number estimates for Multicopy Genes and VNTRs and run PheWAS with various phenotypes, including data generation, data normalization, QC steps and PheWAS analysis. All coordinates used in this pipeline were in UCSC Human Genome hg38.

## Defining CNV regions
VNTR regions were defined from the Simple Repeats track in the UCSC Genome Browser. We selected those VNTRs with motif >=10 bp and span >= 100 bp. Overlapping VNTRs were merged into single regions.
```
cat simpleRepeat.txt| \
  cut -f 2-4,6,7,17 | \
  awk '$1 ~ /^chr[0-9X]+$/ && $3-$2 >= 100 && $4 >= 10' |
  bedtools sort |
  bedtools merge -c 4,6,5 -o distinct,distinct,distinct |
  awk 'BEGIN{OFS = "\t"}{print $1,$2,$3,$1":"$2"-"$3,$4,$5,$6}' > VNTR_100bp_10motif.bed
```

Multicopy genes were defined from the results of CNVnator analysis. We utilized the top most variable genes in 625 [Human Genome Diversity Panel](https://www.internationalgenome.org/data-portal/data-collection/hgdp) cohort. The Exon coordinates of these RefSeq genes were downloaded from UCSC browser and 100bp padding was added each side of the exon.

```
cat Refseq_exon_hg38.txt |
  awk 'BEGIN{OFS ="\t"}{print $1";"$4,$2-100,$3+100}' |
  bedtools sort |
  bedtools merge -i stdin |
  sed -e 's/;/\t/' |
  awk 'OFS = "\t"{print $1,$3,$4,$2}' > Refseq_exon_hg38.bed
```
 
## Defining technical control regions
For control regions for the analysis of Multicopy genes, we utilized a set of “invariant genes”, defined as the 200 least variable (lowest Standard Deviation) genes in [HGDP cohort](https://www.internationalgenome.org/data-portal/data-collection/hgdp) and having pLI scores >0.9. 

Control regions for VNTRs were the 1000 bp flanks on each side. Any part of the flank overlapping other VNTRs were trimmed.
```
cut -f 1-4 VNTR_100bp_10motif.bed |
  bedtools flank -i stdin -g hg38.chrom.sizes -b 1000 |
  bedtools subtract -a stdin -b VNTR_100bp_10motif.bed |
  bedtools slop -i stdin -g hg38.chrom.sizes -b 1 |
  bedtools intersect -wa -wb -a stdin -b VNTR_100bp_10motif.bed | cut -f 1-8|
  awk '$4==$8'| cut -f 1-8 |
  bedtools slop -i stdin -g hg38.chrom.sizes -b -1 |
  awk '$3==$6 || $2==$7'| cut -f 1-4 |
  bedtools sort > VNTR_100bp_10motif_1kb_flanks.bed
```
 
## Read Depth generation and normalization
The read depth for each Multicopy gene, Invariant gene, VNTR and their flanks were generated using [mosdepth](https://github.com/brentp/mosdepth).
```
mosdepth -b "$PREFIX"_region.bed -f $GENOME -n $PREFIX $CRAM
```
The raw read depth were normalized using the script GCbinsAndNorm.r.
 
## Quality Control
- PCA: was generated using prcomp function in R and outliers were removed based on first 10 PCs by manual inspection
- Density plots: were generated using density function in R and samples were outlier from the distribution were removed
 
## PheWAS analysis
- [REGENIE](https://rgcgithub.github.io/regenie/) on binary and quantitative traits on each TOPMed subcohort and Ancestry (runRegenie_binary.sh, runRegenie_quantitative.sh).
- Merging of REGENIE output using [METAL](https://genome.sph.umich.edu/wiki/METAL_Documentation) (METAL_example_script.sh).
