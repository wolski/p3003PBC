# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R package (`p3003PBC`) containing educational materials for the FGCZ (Functional Genomics Center Zurich) Protein Bioinformatics Course. The package includes lecture slides and R exercises covering experimental design, hypothesis testing, linear models, multiplicity correction, sample size estimation, and GSEA.

## Package Structure

- **R/**: Core R functions
  - `utilities.R`: Contains `sample_stats()` function for statistical simulations
  - `BHExtracted.R`: Implementation of Benjamini-Hochberg p-value adjustment
  - `hello.R`: Template function

- **vignettes/**: Course lecture materials as R Markdown (xaringan) presentations
  - `ExpDesign.Rmd`: Experimental design lecture
  - `Hypothesis_Testing_p3003.Rmd`: Hypothesis testing lecture
  - `IntroductionToLinearModels.Rmd`: Linear models lecture
  - `Multiplicity_p3003.Rmd`: Multiple testing correction lecture
  - `SampleSizeEstimation.Rmd`: Sample size estimation lecture
  - `GSEA_with_msig.Rmd`: Gene set enrichment analysis lecture
  - `OtherSoftware.Rmd`: Overview of other proteomics software
  - `BriefIntroductionTooR.Rmd`: R programming introduction
  - `renderAllVignettest.R`: Script to render all vignette presentations

- **inst/**: Additional resources
  - `images/`: Lecture slide images
  - `fgcz_formatting/`: FGCZ branding files (CSS, headers, footers)
  - `BaseR_JG/`: Base R workshop materials and example data

## Common Development Commands

### Building and Installing the Package
```r
# Install package with source code
devtools::install()

# Or from command line
R CMD INSTALL --no-multiarch --with-keep.source .
```

### Rendering Lecture Materials
```r
# Render individual vignette
rmarkdown::render("vignettes/ExpDesign.Rmd")

# Or use the batch rendering script
source("vignettes/renderAllVignettest.R")
```

### Package Development
```r
# Load package for development
devtools::load_all()

# Check package
devtools::check()

# Run tests (if available)
devtools::test()
```

## Technology Stack

- **Presentation Framework**: xaringan (remark.js-based slides)
- **Output Format**: HTML presentations with custom CSS (metropolis theme)
- **Statistical Packages**: coin, lme4, lmerTest, multcomp, pwr, asympTest
- **Data Manipulation**: tidyverse, dplyr, tidyr, purrr
- **Visualization**: ggplot2, ggfortify, GGally, ggbeeswarm
- **Proteomics**: MSqRob (from statOmics), limma (Bioconductor)

## Vignette Rendering Details

All vignettes use xaringan::moon_reader with:
- 16:9 aspect ratio
- metropolis theme
- GitHub-style code highlighting
- Incremental slide numbering
- Custom CSS from `trug-ggplot2.css`

## Important Notes

- This is an educational package, not a production tool
- Vignettes are designed as presentation slides, not analysis reports
- Images are referenced from `../inst/images/` in vignettes
- Package dependencies include both CRAN and Bioconductor packages (see DESCRIPTION:42-44)
