# List of all vignettes
vignettes <- c(
  "ExpDesignShort.Rmd",
  "FASTA_ProteinInference.Rmd",
  "GSEA_Slides.Rmd",
  "Multiplicity_Simple.Rmd",
  "ExpDesign.Rmd"
#  "ExpDesignProteomics.Rmd",
#  "Hypothesis_Testing_p3003.Rmd",
#  "IntroductionToLinearModels.Rmd",
#  "Multiplicity_p3003.Rmd",
#  "SampleSizeEstimation.Rmd",
#  "GSEA_with_msig.Rmd",
#  "OtherSoftware.Rmd",
#  "BriefIntroductionTooR.Rmd",
#  "ImpactOfFilteringOnProteinFDR.Rmd",
#  "BioCED2019.Rmd"
)

# Render all vignettes to HTML
for (v in vignettes) {
  rmarkdown::render(v)
}

# Render index page last (after all vignettes are generated)
rmarkdown::render("index_vignettes.Rmd")

# Print all HTML files to PDF using pagedown
# Use full file path to avoid relative path issues with images
html_files <- sub("\\.Rmd$", ".html", vignettes)

# Create pdf subfolder if it doesn't exist
if (!dir.exists("pdf")) {
  dir.create("pdf")
}

for (h in html_files) {
  if (file.exists(h)) {
    full_path <- normalizePath(h)
    pdf_name <- file.path("pdf", sub("\\.html$", ".pdf", h))
    pagedown::chrome_print(full_path, output = pdf_name, timeout = 20)
  }
}
