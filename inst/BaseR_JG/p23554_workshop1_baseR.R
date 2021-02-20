#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2021
#
#


# R
# this is a comment and will not be excecuted (even if put into the console)
# because of the hash-tag at the beginning!
print("hello world") # this is a print statement

# setting the working directory
# useful in RStudio is also -> Session -> Set working directory -> To source file location
setwd("~/") # this is the path to your home directory
getwd() # get workding directory
if(FALSE){
#  configuration (downloading packages from CRAN repository)
install.packages("BiocManager") # installer for bioconductor packages
install.packages("devtools") # installer for github packages
install.packages("gplots") # with this you can draw nice heatmaps with dendograms
install.packages("matrixStats") # this is used for colMedians

# packages from bioconductor
BiocManager::install("limma") # this is used for plotDesities here
}
# attaching packages
library("gplots")
library("limma")
library("matrixStats")

# listing
list.files()



#
#
#   read in data from a spread sheet & prepare for further analysis
#
#

# read in a whole proteinGroups file (output from MaxQuant tool)
dat <- read.table("p2691_combocourse_proteinGroups.txt", sep = "\t", header = TRUE)

# list all columns of this object
colnames(dat)

# get protein accessions (uniprot form, indistinguishables separated by semi-colon)
myproteins <- dat$Majority.protein.IDs

# get quantitative values
lookUpString <- "Intensity\\."

id_x <- grep(x = colnames(dat), pattern = lookUpString)

# check: do you get the right ones?
colnames(dat)[id_x]

#do it stepwise
better1 <- gsub(x = colnames(dat)[id_x], pattern = "Intensity\\.", replacement = "")
betterColnames <- gsub(x = better1, pattern = "20180313_", replacement = "Sample_")
# do it all in one
betterColnames <- gsub(x = gsub(x = colnames(dat)[id_x], pattern = "Intensity\\.", replacement = ""), pattern = "20180313_", replacement = "Sample_")

# extract conditions from:
#
# > betterColnames
# [1] "Sample_03_GE4_LS_180314100836" "Sample_04_G3_SW_180314112243"  "Sample_05_G2_MB"               "Sample_06_GE1_AP"
# [5] "Sample_08_G4_CM"               "Sample_09_GE3_CW"              "Sample_10_GE2_ZF"              "Sample_11_G1_AP"
# [9] "Sample_12_GE1_JTK"

myConditions <- c("GlycEtOH", "Gluc", "Gluc", "GlycEtOH", "Gluc", "GlycEtOH", "GlycEtOH", "Gluc", "GlycEtOH")
table(myConditions)

# rows to filter out (e.g. Decoy hits, proteins with only 1 or 2 peptides)
minPeptide <- 2
bool_minPepOK <- dat$Razor...unique.peptides > minPeptide

# how many do we take
sum(bool_minPepOK)

# build up the generic matrix
qmat <- dat[bool_minPepOK, id_x]

# add better colnames and proteins as rownames
rownames(qmat) <- myproteins[bool_minPepOK]
colnames(qmat) <- betterColnames

# look at the  structure of the object
str(qmat)

# look at it
par(mfrow = c(1,1)) # only 1 plot per sheet
image(as.matrix(t(qmat))) # --> not too much visible why??

# maybe here it becomes obvious
hist(as.matrix(qmat), breaks = 100)



#
#
#   transformation, handling NAs, normalization, descriptive statistics
#
#


# missing values
# handle zeros as NAs
qmat[qmat == 0] <- NA

# log 2 transformation
lg2mat <- log2(qmat)

# is there an effect
image(as.matrix(t(lg2mat)))
hist(as.matrix(lg2mat), breaks = 100)


#calculate an average expression per protein
lg2_meanExpressionPerProtein <- log2(rowMeans(qmat, na.rm = TRUE))

# be aware taking average of the qmat and then log2 is NOT the same as averaging the log2Values
# plot(log2(rowMeans(qmat,na.rm = TRUE)) ~ rowMeans(lg2mat, na.rm = TRUE))

bool_hasNA <- rowSums(is.na(lg2mat)) > 0
sum(bool_hasNA)

# where are our NAs
par(mfrow = c(2,1)) # 2 plots as 2 rows and 1 column
hist(lg2_meanExpressionPerProtein, breaks = 100, xlim = c(20,35))
hist(lg2_meanExpressionPerProtein[bool_hasNA], breaks = 100, xlim = c(20,35))
# conclusion ->  NAs are mainly found in the low intensity range (close to limit of detection = LOD)

# filter out NAs to work with completely filled matrix
bool_noNA <- rowSums(is.na(lg2mat)) == 0
sum(bool_noNA)
hist(lg2_meanExpressionPerProtein[bool_noNA], breaks = 100, xlim = c(20,35))
hist(lg2_meanExpressionPerProtein[bool_hasNA], breaks = 100, xlim = c(20,35))

#
# clean log2 matrix
#
clean_lg2Mat <- lg2mat[bool_noNA, ]

# lets look at the different samples
par(mfrow = c(1,1))
limma::plotDensities(clean_lg2Mat, legend = "bottomright")
# conclusion -> obviously not all samples are having identical means/medians


# Normalization (subtract median values for each column)
# find("colMedians")
#  [1] "package:matrixStats"

sampleMedians <- matrixStats::colMedians(as.matrix(clean_lg2Mat))
norm_lg2Mat <- clean_lg2Mat - sampleMedians

limma::plotDensities(norm_lg2Mat, legend = "topleft")

#
# Descriptive Statistics (QC)
#

myCor <- cor(norm_lg2Mat)
image(myCor) # visualize correlations

# use heatmap.2 to get a much nicer picture of your correlations
heatmap.2(myCor, dendrogram = "both", trace = "none", margins = c(15,15))
# conclusion -> at least we see that GE and G are well separated

heatmap.2(as.matrix(norm_lg2Mat), dendrogram = "both",trace = "none", margins = c(15,15))
# conclusion -> clustering is not as nice as for the correlations

# multi dimensional scaling
plotMDS(as.matrix(norm_lg2Mat))



# mean expression per group
GEgroup_idx <- grep(x = colnames(norm_lg2Mat), pattern = "_GE")
Ggroup_idx <- setdiff(1:ncol(norm_lg2Mat),GEgroup_idx)
colnames(norm_lg2Mat)[GEgroup_idx]
colnames(norm_lg2Mat)[Ggroup_idx]

GE_proteinMeans <- rowMeans(norm_lg2Mat[GEgroup_idx])
G_proteinMeans <- rowMeans(norm_lg2Mat[Ggroup_idx])

lg2_foldChanges_GvsGE <- G_proteinMeans - GE_proteinMeans

# Bland-Altmann (MA-plot)

plot(lg2_foldChanges_GvsGE ~ lg2_meanExpressionPerProtein[bool_noNA])
abline(h = c(-1,1), col = "red")
abline(h = 0, col = "blue")


# now let's test what proteins are differentially expressed using t-test on
pValueVector <- vector(length = nrow(norm_lg2Mat))

i <- 20
t.test(x = norm_lg2Mat[i, GEgroup_idx], y = norm_lg2Mat[i, Ggroup_idx], paired = FALSE)
res <- t.test(x = norm_lg2Mat[i, GEgroup_idx], y = norm_lg2Mat[i, Ggroup_idx], paired = FALSE)
res$p.value

# usin for loop
for (i in 1:length(pValueVector)) {
  res <- t.test(x = norm_lg2Mat[i, GEgroup_idx], y = norm_lg2Mat[i, Ggroup_idx], paired = FALSE)
  pValueVector[i] <- res$p.value
}

# apply
pV <- sapply(1:nrow(norm_lg2Mat), function(x) t.test(norm_lg2Mat[x,GEgroup_idx], norm_lg2Mat[x,Ggroup_idx])$p.value)

pV <- lapply(1:nrow(norm_lg2Mat), function(x) t.test(norm_lg2Mat[x,GEgroup_idx], norm_lg2Mat[x,Ggroup_idx]))
sapply(pV, function(x){x$p.value})
sapply(pV, function(x){x$statistic})
sapply(pV, function(x){x$estimate})
sapply(pV, function(x){x$stderr})

# do we get the same
plot(-1*log10(pV) , -1*log10(pValueVector))

# volcano plot
plot(-1*log10(pValueVector) ~ lg2_foldChanges_GvsGE)






