Class 13: Genome Informatics
================

The purpose of this lab session is to introduce a set of tools used in high- throughput sequencing and the process of investigating interesting gene variance in Genomics. High-throughput sequencing is now routinely applied to gain insight into a wide range of important topics in biology and medicine.

There are a number of gene variants associated with childhood asthma. A study from Verlaan et al. (2009) shows that 4 candidate SNPs demonstrate significant evidence for association.

What are those 4 candidate SNPs? rs12936231, rs8067378, rs9303277, and rs7216389

What three genes do these variants overlap or effect? ZPBP2, GSDMB, and ORMDL3

What is the location of rs8067378 and what are the different alleles for rs8067378? Chromosome 17: 39895095 (forward strand) A/G variants (43% G)

What proportion of the Mexican Ancestry in Los Angeles sample population (MXL) are homozygous for the asthma associated SNP (G|G)?

``` r
mxl <- read.csv("373531-SampleGenotypes-Homo_sapiens_Variation_Sample_rs8067378.csv")
head(mxl)
```

    ##   Sample..Male.Female.Unknown. Genotype..forward.strand. Population.s.
    ## 1                  NA19648 (F)                       A|A ALL, AMR, MXL
    ## 2                  NA19649 (M)                       G|G ALL, AMR, MXL
    ## 3                  NA19651 (F)                       A|A ALL, AMR, MXL
    ## 4                  NA19652 (M)                       G|G ALL, AMR, MXL
    ## 5                  NA19654 (F)                       G|G ALL, AMR, MXL
    ## 6                  NA19655 (M)                       A|G ALL, AMR, MXL
    ##   Father Mother
    ## 1      -      -
    ## 2      -      -
    ## 3      -      -
    ## 4      -      -
    ## 5      -      -
    ## 6      -      -

``` r
table(mxl$Genotype..forward.strand.)
```

    ## 
    ## A|A A|G G|A G|G 
    ##  22  21  12   9

``` r
# proportion
table(mxl$Genotype..forward.strand.)/nrow(mxl)
```

    ## 
    ##      A|A      A|G      G|A      G|G 
    ## 0.343750 0.328125 0.187500 0.140625

14 % of the individuals are homozygous for the G allele at this SNP locus.

Quality scores in fastq files
-----------------------------

``` r
library(seqinr)
library(gtools)
phred <- asc( s2c("DDDDCDEDCDDDDBBDDDCC@") ) - 33
phred
```

    ##  D  D  D  D  C  D  E  D  C  D  D  D  D  B  B  D  D  D  C  C  @ 
    ## 35 35 35 35 34 35 36 35 34 35 35 35 35 33 33 35 35 35 34 34 31

Galaxy server IP: 149.165.171.141

Population Scale Analysis
-------------------------

One sample is obviously not enough to know what is happening in a population. You are interested in assessing genetic differences on a population scale. So, you processed about ~230 samples and did the normalization on a genome level. Now, you want to find whether there is any association of the 4 asthma-associated SNPs (rs8067378...) on *ORMDL3* expression.

``` r
population <- read.table("https://bioboot.github.io/bggn213_S19/class-material/rs8067378_ENSG00000172057.6.txt", header = TRUE)
head(population)
```

    ##    sample geno      exp
    ## 1 HG00367  A/G 28.96038
    ## 2 NA20768  A/G 20.24449
    ## 3 HG00361  A/A 31.32628
    ## 4 HG00135  A/A 34.11169
    ## 5 NA18870  G/G 18.25141
    ## 6 NA11993  A/A 32.89721

Determine the sample size for each genotype and their corresponding median expression levels for each of these genotypes.

``` r
table(population$geno)
```

    ## 
    ## A/A A/G G/G 
    ## 108 233 121

``` r
summary(population)
```

    ##      sample     geno          exp        
    ##  HG00096:  1   A/A:108   Min.   : 6.675  
    ##  HG00097:  1   A/G:233   1st Qu.:20.004  
    ##  HG00099:  1   G/G:121   Median :25.116  
    ##  HG00100:  1             Mean   :25.640  
    ##  HG00101:  1             3rd Qu.:30.779  
    ##  HG00102:  1             Max.   :51.518  
    ##  (Other):456

``` r
# look at the summary statistics for G/G individuals only
inds <- population$geno == "G/G"
summary(population[inds,])
```

    ##      sample     geno          exp        
    ##  HG00099:  1   A/A:  0   Min.   : 6.675  
    ##  HG00109:  1   A/G:  0   1st Qu.:16.903  
    ##  HG00112:  1   G/G:121   Median :20.074  
    ##  HG00116:  1             Mean   :20.594  
    ##  HG00118:  1             3rd Qu.:24.457  
    ##  HG00120:  1             Max.   :33.956  
    ##  (Other):115

``` r
# look at the summary statistics for A/A individuals only
inds2 <- population$geno == "A/A"
summary(population[inds2,])
```

    ##      sample     geno          exp       
    ##  HG00096:  1   A/A:108   Min.   :11.40  
    ##  HG00100:  1   A/G:  0   1st Qu.:27.02  
    ##  HG00101:  1   G/G:  0   Median :31.25  
    ##  HG00102:  1             Mean   :31.82  
    ##  HG00104:  1             3rd Qu.:35.92  
    ##  HG00105:  1             Max.   :51.52  
    ##  (Other):102

The median expression in A/A individuals is somewhat higher.

``` r
# look at the summary statistics for A/G individuals only
inds3 <- population$geno == "A/G"
summary(population[inds3,])
```

    ##      sample     geno          exp        
    ##  HG00097:  1   A/A:  0   Min.   : 7.075  
    ##  HG00103:  1   A/G:233   1st Qu.:20.626  
    ##  HG00106:  1   G/G:  0   Median :25.065  
    ##  HG00110:  1             Mean   :25.397  
    ##  HG00114:  1             3rd Qu.:30.552  
    ##  HG00115:  1             Max.   :48.034  
    ##  (Other):227

A/G individuals seem to have an intermediate expression phenotype.

Make a boxplot showing each of the three genotypes.

``` r
# y ~ grp
boxplot(exp ~ geno, data = population, notch = TRUE)
```

![](Class13_files/figure-markdown_github/unnamed-chunk-11-1.png)
