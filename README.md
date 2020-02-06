## genotypeErrors
The rules of transmission genetics are applied at the genomic scale in this novel method to estimate genotyping error rates at biallelic loci.

Description of contents:

#### error_rates.R
  R code that uses counts for each genotype-phase combination in a three-generation pedigree with five individuals to calculate genotyping error rates
  
#### genoCombinationUtil.py
  Converts .vcf file with genotypes and .bed file with phases to a .txt file with genotype-phase combinations for input into R.
  
#### genotypeSims.py
  Simulates genotype-phase combinations from a three-generation pedigree assuming a neutral site frequency spectrum.
