# G4SNVHunter <img src="./vignettes/images/logo.png" align = "right" width = "150" />

G4SNVHunter is an R package leveraging the G4Hunter algorithm to systematically 
identify single nucleotide variants (SNVs), 
along with other small-scale variants such as indels and MNVs,
that have the potential to disrupt G-quadruplex (G4) formation propensity.


# Installation
### Option 1: Install from GitHub

You can install the package directly from GitHub,
```r
# install.packages("devtools")
devtools::install_github("rongxinzh/G4SNVHunter")
```

To run the sample code in our [vignette](
https://rongxinzh.github.io/G4SNVHunter/G4SNVHunter.html
), set the `dependencies` parameter to `TRUE`,
```r
# install.packages("devtools")
devtools::install_github("rongxinzh/G4SNVHunter", dependencies = TRUE)
```

**NOTE**

* Your R version must be &#8805; 4.3.
* If you are using Windows, please install 
[Rtools](https://ohdsi.github.io/Hades/rSetup.html#Installing_RTools) 
before proceeding.

### Option 2: Install from Bioconductor

The package is also available on Bioconductor as a development (devel) version. 
To install it, follow these steps,
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version = 'devel')

BiocManager::install("G4SNVHunter")
```

If you want to run the sample code in our [vignette](
https://rongxinzh.github.io/G4SNVHunter/G4SNVHunter.html
), 
set `dependencies` to `TRUE`,

```r
BiocManager::install("G4SNVHunter", dependencies = TRUE)
```

**NOTE**

* The Bioconductor installation (`version = 'devel'`) requires the latest 
version of R.
* Make sure to initialize Bioconductor devel using 
`BiocManager::install(version = 'devel')` before installing the package.

# Quick Start

First, you need to library the G4SNVHunter package,

```r
library(G4SNVHunter)
```

In this example, we need to use the human genome sequence of chromosome 16 
(hg19 version).

```r
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
}
library(BSgenome.Hsapiens.UCSC.hg19)

# Load sequence for chromosome 16 (hg19)
hg19 <- BSgenome.Hsapiens.UCSC.hg19

chr16_seq <- DNAStringSet(hg19$chr16)

# Chromosome names are needed for analysis
names(chr16_seq) <- "chr16"

```

We have prepared some sample variant data, which can be loaded as,

```r
vcf_file <- system.file("extdata", 
                        "example_variants_chr16.vcf", 
                        package = "G4SNVHunter")
variants <- loadVariant(vcf_file, file_type = "vcf")
seqlevels(variants) <- paste0("chr", seqlevels(variants))
```

Before assessing the impact of SNV on G4, we need to predict the G4-prone 
regions,

```r
chr16_G4 <- G4HunterDetect(chr16_seq)
```

Then, the impact of variants on the formation of G4 structures 
can be evaluated as,

```r
result <- G4VarImpact(G4 = chr16_G4, 
                      variants = variants, 
                      ref_col = "REF",
                      alt_col = "ALT")
```

Since not all variants can significantly impact the formation of G4s, we need 
to filter the output. While the filtering parameters can be set flexibly, 
a set of recommended parameters is:

```r
filtered_var_eff <- filterVarImpact(result, 
                                    raw_score_threshold = 1.5,
                                    mut_score_threshold = 1.2)

print(filtered_var_eff)
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
