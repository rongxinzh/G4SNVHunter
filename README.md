# G4SNVHunter <img src="./vignettes/images/logo.png" align = "right" width = "150" />

G4SNVHunter is an R package leveraging the G4Hunter algorithm to systematically 
identify single nucleotide variants (SNVs), including single nucleotide 
polymorphisms (SNPs), with the potential to disrupt the formation of
G-quadruplex (G4) structures.

**Please note that this is a development version. For the stable release, 
please install the version from the [main](
https://github.com/rongxinzh/G4SNVHunter/tree/main) branch.**

# Installation

You can install directly from GitHub using:
```r
# install.packages("devtools")
devtools::install_github("rongxinzh/G4SNVHunter")
```

Please note that your R version needs to be &#8805; 4.3.

If your operating system is Windows, you will need to install 
[Rtools](https://ohdsi.github.io/Hades/rSetup.html#Installing_RTools) on your 
system first.

# Quick Start

First, you need to library the G4SNVHunter package,

```r
library(G4SNVHunter)
```

In this example, we need to use the human genome sequence of chromosome 21 
(hg19 version).

```r
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
}
library(BSgenome.Hsapiens.UCSC.hg19)

# Load sequence for chromosome 21 (hg19)
hg19 <- BSgenome.Hsapiens.UCSC.hg19

chr21_seq <- DNAStringSet(hg19$chr21)

# Chromosome names are needed for analysis
names(chr21_seq) <- "chr21"

```

We have prepared some sample variant data, which can be loaded as,

```r
data(snv_gr)
```

Before assessing the impact of SNV on G4, we need to predict the G4-prone 
regions,

```r
chr21_G4 <- G4HunterDetect(chr21_seq)
```

Then, the impact of SNVs on the formation of G4 structures can be evaluated as,

```r
snv_eff <- SNVImpactG4(chr21_G4, snv_gr, alt_col = "alt")
```

Since not all SNVs can significantly impact the formation of G4s, we need to 
filter the output. While the filtering parameters can be set flexibly, a set of 
recommended parameters is:

```r
filtered_snv_eff <- filterSNVImpact(snv_eff, 
                                    raw_score_threshold = 1.5,
                                    mut_score_threshold = 1.2)

print(filtered_snv_eff)
```

Generally, G4 structures with an absolute G4Hunter score above 1.5 are 
considered to have a high confidence in forming stable structures. In contrast, 
those with scores below 1.2, or more conservatively, below 1.0, are typically 
less likely to form stable structures.

# Vignette

For full use of our package, please refer to our [vignette](
https://rongxinzh.github.io/G4SNVHunter/G4SNVHunter.html)
(highly recommended). 

# Bug Reports

If you encounter any issues, have questions, or would like to make suggestions, 
feel free to report them at our [bug tracker](
https://github.com/rongxinzh/G4SNVHunter/issues).

# Contact

For additional inquiries, contact us at: 
Email: rongxinzhang@outlook.com
