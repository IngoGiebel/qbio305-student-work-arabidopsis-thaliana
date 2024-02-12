########################################################################
##                            GROUP 5                                 ##
########################################################################

library(PopGenome)
library(ggplot2)
library(vcfR)
library (readr)
library(dplyr)
library(tidyr)
library(tibble)
library(BiocManager)
library(KEGGREST)
library(org.At.tair.db)
library(Rgraphviz)
library(topGO)
library(biomaRt)
library(ggplot2)
library(AnnotationDbi)
library(clusterProfiler)
library(scales)
library(readr)
library(GenomicRanges)
library(dplyr)
library(adegenet)




setwd("C:/Users/alizf/Documents/Semester3/Quantitative Population & Genetics/PART2/Arabidopsis project/Data Analysis")
at.VCF <- read.vcfR("group_5_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz")



########################################################################
##              Principle Component Analysis (PCA)                    ##
########################################################################

# Convert VCF to genind object
genind_vcf <- vcfR2genind(at.VCF)

# Scale genind object for PCA
genind_vcf_scaled = scaleGen(genind_vcf, NA.method = "mean")

# Perform PCA
pca <- dudi.pca(genind_vcf_scaled, cent = TRUE, scale = FALSE, scannf = FALSE, nf = 10)

pca_axis <- pca$li

# You can also directly provide a list of names
# first twenty are from sweden, last twenty from spain
pop <- c("SWE" , "SWE" , "SWE" , "SWE" , "SWE" , "SWE" , "SWE" , "SWE" , "SWE" , "SWE" , "SWE" , "SWE" , "SWE" , "SWE" , "SWE" , "SWE" , "SWE" , "SWE" , "SWE" , "SWE" , "ESP" , "ESP" , "ESP" , "ESP" , "ESP" , "ESP" , "ESP" , "ESP" , "ESP" , "ESP", "ESP" , "ESP" , "ESP" , "ESP" , "ESP" , "ESP" , "ESP" , "ESP" , "ESP" , "ESP")
pca_axis$Population <- pop

# remake data.frame
pca_2 <- as_tibble(data.frame(pca_axis, pop))

# first convert to percentage variance explained
pca_eigenval <- pca$eig[1:10]
pve <- data.frame(PC = 1:10, pve = pca_eigenval/sum(pca_eigenval)*100)


## Plot PCA

# plot pca PC1 and PC2
ggplot(pca_2, aes(Axis1, Axis2, col = pop)) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("red", "blue")) +
  coord_equal() +
  theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  ggtitle("Arabidopsis thaliana Accessions from Spain (ESP) and Sweden (SWE)") +
  theme(plot.title = element_text(hjust = 0.5, size=15),
        plot.margin = margin(20, 20, 20, 20))







# Start and stop positions of chromosomes
"
  CHROM    Start      End
1     1   443538 30404145
2     2   316267 12415019
3     3 14946344 23400776
4     4   931554 18536293
5     5   160635 26959553
"

At_Chr1 <- readVCF("group_5_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz", numcols=89, tid="1", frompos = 443538, topos = 30404145, include.unknown = TRUE)
chromname <- "1"


# Number of biallelix abd polyallelic sites
At_Chr1@n.biallelic.sites + At_Chr1@n.polyallelic.sites

# To check starting position and last position of genome class
At_Chr1@region.names

#####Define populations in dataset####

population_info <- read_delim("pop_group5.txt", delim = "\t")
populations <- split(population_info$sample, population_info$pop)

At_Chr1 <- set.populations(At_Chr1, populations, diploid = T)
At_Chr1@populations

####Setting up sliding windows###

# total number of sites
At_Chr1@n.sites
chr1 <- At_Chr1@n.sites + 1

# set window size and window jump
window_size <- 100
window_jump <- 50

# use seq to find the start points of each window
window_start_Chr1 <- seq(from = 1, to = chr1, by = window_jump)
# add the size of the window to each start point
window_stop_Chr1 <- window_start_Chr1 + window_size

# remove windows from the start and stop vectors
window_start_Chr1 <- window_start_Chr1[which(window_stop_Chr1 < chr1)]
window_stop_Chr1 <- window_stop_Chr1[which(window_stop_Chr1 < chr1)]

# save as a data.frame
windows_Chr1 <- data.frame(start = window_start_Chr1, stop = window_stop_Chr1,
                           mid = window_start_Chr1 + (window_stop_Chr1-window_start_Chr1)/2)

# sliding window dataset
At_sw_Chr1 <- sliding.window.transform(At_Chr1, width = 100, jump = 50, type = 2)


####### Calculating sliding window estimates of nucleotide diversity and differentiation #####

# Pi
At_sw_Chr1 <- diversity.stats(At_sw_Chr1, pi = TRUE)

# FST
At_sw_Chr1 <- F_ST.stats(At_sw_Chr1, mode = "nucleotide")



#### calculate neutrality statistics ####
At_sw_Chr1 <- neutrality.stats(At_sw_Chr1)


####Extracting statistics for visualization####

# Pi
nd_Chr1 <- At_sw_Chr1@nuc.diversity.within/100

# FST
fst_Chr1 <- t(At_sw_Chr1@nuc.F_ST.pairwise)

# dxy
dxy_Chr1 <- get.diversity(At_sw_Chr1, between = T)[[2]]/100

# Tajma's D
td_Chr1 <- At_sw_Chr1@Tajima.D/100


# make population name vector
pops <- c("SWE","ESP")
colnames(nd_Chr1) <- paste0(pops, "_pi")

# get column names
x <- colnames(fst_Chr1)

# Loop through each population and replace the corresponding population name in the column names
for (i in 1:length(pops)) {
  pattern <- paste0("pop", i)
  x <- sub(pattern, pops[i], x)
}

# replace forward slash
x <- sub("/", "_", x)

paste0(x, "_fst")
paste0(x, "_dxy")

colnames(fst_Chr1) <- paste0(x, "_fst")
colnames(dxy_Chr1) <- paste0(x, "_dxy")
colnames(td_Chr1) <- paste0(pops, "_td")

# Create a tibble
At_data_Chr1 <- as_tibble(data.frame(windows_Chr1, td_Chr1, nd_Chr1, fst_Chr1, dxy_Chr1))





###############################################################################
### B O X P L O T S  ###
###############################################################################



###
# nucleotide diversity (Pi) Boxplot
###

# select nucleotide diversity data and calculate means
At_data_Chr1 %>% select(contains("pi")) %>% summarise_all(mean)
pi_g_Chr1 <- At_data_Chr1 %>% select(contains("pi")) %>% gather(key = "populations", value = "pi")
pi_g_Chr1$log_pi <- log10(pi_g_Chr1$pi)


# Wilcoxon test
wilcox_test_pi_Chr1 <- wilcox.test(log_pi ~ populations, data = pi_g_Chr1)
print(wilcox_test_pi_Chr1)

comparison_data_Chr1 <- pi_g_Chr1 %>% filter(populations %in% c("SWE_pi", "ESP_pi"))

# Perform Wilcoxon rank-sum test
wilcox_test_pi_Chr1 <- wilcox.test(log_pi ~ populations, data = comparison_data_Chr1)
print(wilcox_test_pi_Chr1)


# Boxplot Nucleotide diversity
a_pi_Chr1 <- ggplot(pi_g_Chr1, aes(populations, log_pi, fill = populations)) +
  geom_boxplot(color = "black") +  # Border color
  #scale_fill_manual(values = c("red", "blue","orange","magenta")) +  # Box fill colors
  theme_light() +
  xlab(NULL) +
  ylab("Log10(pi)") +
  ggtitle(paste("Chromosome", chromname)) +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, size = 55),
        axis.title.x = (element_text(size=43)),
        axis.title.y = (element_text(size=43)),
        axis.text.y = (element_text(size=43)),
        axis.text.x = (element_text(size=43))
        ) +
  annotate("text", x = 1.5, y = max(pi_g_Chr1$log_pi)-0.2, label = paste("p =", format.pval(wilcox_test_pi_Chr1$p.value, digits = 3)), size=28)

a_pi_Chr1


###
# TajimaD (td) Boxplot
###

At_data_Chr1 %>% select(contains("td_Chr1")) %>% summarise_all(mean)
td_g_Chr1 <- At_data_Chr1 %>% select(contains("td")) %>% gather(key = "populations", value = "td")

# Remove rows with missing values
td_g2_Chr1 <- na.omit(td_g_Chr1)

# Log transformation
td_g2_Chr1$log_td <- log10(td_g2_Chr1$td)

# Wilcoxon test
wilcox_test_td_Chr1 <- wilcox.test(log_td ~ populations, data = td_g2_Chr1)

comparison_data_Chr1 <- td_g2_Chr1 %>% filter(populations %in% c("SWE_td", "ESP_td"))

# Perform Wilcoxon rank-sum test
wilcox_test_td_Chr1 <- wilcox.test(log_td ~ populations, data = comparison_data_Chr1)
print(wilcox_test_td_Chr1)


# Boxplot of Tajimas D
a2_td_g2_Chr1 <- ggplot(td_g2_Chr1, aes(populations, log_td, fill = populations)) +
  geom_boxplot(color = "black") +  # Border color
  theme_light() +
  xlab("Log10(td)") +
  ylab("Log10(Tajima's D)") +
  ggtitle(paste("Tajima's D of", chromname)) +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, size = 55),
        axis.title.x = (element_text(size=43)),
        axis.title.y = (element_text(size=43)),
        axis.text.y = (element_text(size=43)),
        axis.text.x = (element_text(size=43))
  ) +
  annotate("text", x = 1.5, y = -1, label = paste("Wilcox test p =", format.pval(wilcox_test_td_Chr1$p.value, digits = 3)), size=10)


a2_td_g2_Chr1


#############################################
##             D E N S I T Y               ##
#############################################

###
# Tajima's D density
###

subset_ESP_td <- log(At_data_Chr1$ESP_td[!is.na(At_data_Chr1$ESP_td) & !is.nan(log(At_data_Chr1$ESP_td))])
subset_SWE_td <- log(At_data_Chr1$SWE_td[!is.na(At_data_Chr1$SWE_td) & !is.nan(log(At_data_Chr1$SWE_td))])

tajimas_d_g <- At_data_Chr1 |>
  dplyr::select(tidyselect::contains("td")) |>
  tidyr::gather(key = "population", value = "td", na.rm = TRUE)

tajimas_d_g$log_td <- log10(tajimas_d_g$td)


# Kruskal-Wallis test
td_data <- list(
  ESP = subset_ESP_td,
  SWE = subset_SWE_td
)

kruskal_td_dist <- kruskal.test(td_data)
print(kruskal_td_dist)


ggplot2::ggplot(
  tajimas_d_g,
  ggplot2::aes(x = log_td)) +
  ggplot2::geom_density(
    ggplot2::aes(fill = population),
    alpha = 0.5,
    color = "black"
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Log10(Tajima's D)",
    title = paste("Chromosome", chromname),
    fill = "Population"
    ) +
  theme_classic() +
  ggplot2::theme(
    plot.title = element_text(hjust = 0.5, size = 55),
    axis.text = element_text(size = 30),
    axis.title.y = (element_text(size=30)),
    legend.title = element_text(size=30),
    legend.text = element_text(size=20),
    legend.key.height = unit(2, "lines"),
  ) +
  xlim(-3.8, -1.3) +
  annotate("text", x = -3, y = 1.8, label = paste("p =", format.pval(kruskal_td_dist$p.value, digits = 3)), size=28)



################################################
##           F A C E T P L O T S              ##
################################################

## FST and dxy

# Select data of interest
hs_fst_dxy_Chr1 <- At_data_Chr1 %>%
  select(mid, SWE_ESP_dxy, SWE_ESP_fst)

# Mutate to set FST and dXY values smaller than zero to zero
hs2_fst_dxy_Chr1 <- hs_fst_dxy_Chr1 %>%
  select(mid, SWE_ESP_dxy, SWE_ESP_fst) %>%
  mutate(across(c(SWE_ESP_fst, SWE_ESP_dxy),
                ~ ifelse(. < 0, 0, .)))

hs <- At_data_Chr1 %>%
  select(mid, ESP_pi, SWE_pi, SWE_ESP_fst, SWE_ESP_dxy) %>%
  mutate(across(c(SWE_ESP_fst, SWE_ESP_dxy), ~ ifelse(. < 0, 0, .)))

# Use gather to rearrange everything
hs_fst_dxy_g_Chr1 <- gather(hs2_fst_dxy_Chr1, -mid, key = "stat", value = "value")

# Reorder the levels of the stat factor
hs_fst_dxy_g_Chr1$stat <- factor(hs_fst_dxy_g_Chr1$stat, levels = c("SWE_ESP_fst","SWE_ESP_dxy"))

# Take the logarithm of the value variable
hs_fst_dxy_g_Chr1$log_value <- log10(hs_fst_dxy_g_Chr1$value)

# FACETPLOT
a_fst_dxy_g_Chr1 <- ggplot(hs_fst_dxy_g_Chr1, aes(mid / 10^6, log_value, colour = stat)) + geom_line() +
  facet_wrap(~stat, scales = "free_y", ncol = 1) +
  xlab("Position (Mb)") +
  theme_light() +
  ggtitle(paste("Chromosome", chromname)) +
  theme(legend.position = "none",
        plot.title = element_text(size=25),
        strip.text = element_text(size = 15))


# Show the plot
a_fst_dxy_g_Chr1






##################################
#  Plotting FST along Chromosome #
#  Visualizing FST outliers      #
##################################

# Filter out zero and negative values
positive_data <- At_data_Chr1[At_data_Chr1$SWE_ESP_fst > 0, ]

# Set threshold for 95% and 99%
threshold_95 <- quantile(positive_data$SWE_ESP_fst, 0.95, na.rm = TRUE)
threshold_99 <- quantile(positive_data$SWE_ESP_fst, 0.99, na.rm = TRUE)

# HISTOGRAM
ggplot(positive_data, aes(x = SWE_ESP_fst)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black",
                 alpha = 0.7) +
  labs(title = paste("Chromosome", chromname), x = "SWE_ESP_fst", y
       = "Frequency") +
  geom_vline(xintercept = threshold_95, colour = "orange", linetype =
               "dashed", size = 1)+
  geom_vline(xintercept = threshold_99, colour = "red", linetype =
               "dashed", size = 1) +
  theme_light() +
  theme(plot.title = element_text(size = 25),
        axis.text = element_text(size = 10),
        axis.title.y = (element_text(size=12)),
        axis.title.x = (element_text(size=12)))


###################################
##      Fst Scatterplot          ##
###################################

# Add outlier columns based on thresholds
positive_data$outlier_95 <- ifelse(positive_data$SWE_ESP_fst >
                                     threshold_95, "Outlier", "Non-outlier")
positive_data$outlier_99 <- ifelse(positive_data$SWE_ESP_fst >
                                     threshold_99, "Outlier", "Non-outlier")

# Estimate mean FST value
mean_fst_Chr1 <- mean(positive_data$SWE_ESP_fst)



# Create a plot with marked outliers, FST mean line, and customized legend
ggplot(positive_data, aes(mid, SWE_ESP_fst)) +
  geom_point(aes(colour = ifelse(outlier_99 == "Outlier", "Outlier 99%", ifelse(outlier_95 == "Outlier", "Outlier 95%", "Non-outlier"))), size = 3) +
  scale_colour_manual(values = c("Non-outlier" = "black", "Outlier 95%" = "orange", "Outlier 99%" = "red")) +
  geom_hline(yintercept = mean_fst_Chr1, colour = "blue") +
  xlab("Position on Chromosome (Mb)") + ylab(expression(italic(F)[ST])) +
  theme_light() +
  ggtitle(paste("Locus Specific Estimate of FST on Chromosome", chromname)) +
  guides(color = guide_legend(override.aes = list(size = 3, shape = 16), title = NULL)) +
  theme(plot.title = element_text(size=25),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title.y = (element_text(size=15)),
        axis.title.x = (element_text(size=15)))



###
# Filter and print coordinates of outliers
###

# Filter rows where the 'outlier_95' column has the value "Outlier"
out_95 <- filter(positive_data, outlier_95 == "Outlier")

# Filter rows where the 'outlier_99' column has the value "Outlier"
out_99 <- filter(positive_data, outlier_99 == "Outlier")

# Print the first 20 rows of the 'out_95' data frame
print(out_95, n = 20)

# Print the first 20 rows of the 'out_99' data frame
print(out_99, n = 30)
###################################

out_99 %>% select(start, stop, mid, SWE_ESP_fst, outlier_99)
out_95 %>% select(start, stop, mid, SWE_ESP_fst, outlier_95)

#write.table(out_99, file = "out_99_Chromosome4.csv", sep = ",", row.names = FALSE)
#write.table(out_95, file = "out_95_Chromosome4.csv", sep = ",", row.names = FALSE)





################################################################################
### GENE ENRICHTMENT
################################################################################


############################
# Prepare input file for   #
# Gene Enrichment Analysis #
############################

######
# Running stacks.sh in the terminal!
######

############
# Rscript to process stacks output file "*.p.fst_SWE-ESP.tsv"
############

# read --fstats output file from stacks
fstats <- read.table("group_5_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.p.fst_SWE-ESP.tsv", header=FALSE)

# Set the column names manually and replace spaces with underscores
colnames(fstats) <- gsub(" ", "_", c(
  "Locus ID", "Pop 1 ID", "Pop 2 ID", "Chr", "BP", "Column", "Fishers P",
  "Odds Ratio", "CI Low", "CI High", "LOD", "AMOVA Fst", "Smoothed AMOVA Fst",
  "Smoothed AMOVA Fst P-value"
))

# Create a subset of the dataframe with selected columns
fstats_fst <- data.frame(
  Chr = fstats$Chr,
  BP = fstats$BP,
  Fishers_P = fstats$Fishers_P,
  AMOVA_Fst = fstats$AMOVA_Fst
)


#############
# Read annotation file
#############
annot<-read.delim("At_defense_only.gff", header = FALSE)
head(annot)
names(annot)

print(head(annot, n=1))

# Set the column names manually and replace spaces with underscores
colnames(annot) <- c(
  "gene_id", "chr", "tair_version", "type", "start", "end", "empty1",
  "strand", "empty2", "gene_id2", "gene_id3", "type_name", "short_description", "curation"
)


# Create a subset of the dataframe with selected columns
annot_subset <- data.frame(
  gene_id = annot$gene_id,
  chr = annot$chr,
  start = annot$start,
  end = annot$end
)

print(head(annot_subset, n=5))

#remove "Chr" from each entry in the annot_subset$chr column to make it similar as in fstats_fst$Chr:
annot_subset$chr <- sub("Chr", "", annot_subset$chr)


# Convert data frames to GRanges objects
gr_fst <- GRanges(seqnames = fstats_fst$Chr,
                  ranges = IRanges(start = fstats_fst$BP, end = fstats_fst$BP),
                  fst = fstats_fst$AMOVA_Fst,
                  p_value = fstats_fst$Fishers_P)

length(fstats_fst$Chr)
print(head(gr_fst, n=1))

gr_annotation <- GRanges(seqnames = annot_subset$chr,
                         ranges = IRanges(start = annot_subset$start, end = annot_subset$end),
                         gene_id = annot_subset$gene_id)

# Find overlaps between fst_data and annotation_data
ovl <- findOverlaps(gr_fst, gr_annotation, type = "any", select = "all", ignore.strand = TRUE)

# Extract Gene IDs based on overlaps
overlapping_genes <- as.data.frame(gr_annotation[subjectHits(ovl)])
overlapping_fst <- as.data.frame(gr_fst[queryHits(ovl)])

# Print or further process the results
print(head(overlapping_genes,n=5))
print(head(overlapping_fst,n=5))

# Rename 'seqnames' to 'chrom'
colnames(overlapping_genes)[colnames(overlapping_genes) == "seqnames"] <- "chrom"
colnames(overlapping_fst)[colnames(overlapping_fst) == "seqnames"] <- "chrom"

# Remove the 'strand' column
overlapping_genes <- overlapping_genes[, !(names(overlapping_genes) %in% "strand")]
overlapping_fst <- overlapping_fst[, !(names(overlapping_fst) %in% "strand")]

# Print or further process the results
print(head(overlapping_genes,n=5))
print(head(overlapping_fst,n=5))

# Write both overlapping_genes and overlapping_fst as text files
write.table(overlapping_genes, file = "overlapping_genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(overlapping_fst, file = "overlapping_fst.txt", sep = "\t", row.names = FALSE, quote = FALSE)


####################
# Run this part in terminal
###################

# Copy these file to Cheops and run bedtools intersect to get loci which are within
# genomic boundaries of defense related genes
$ module load bedtools/
  $ bedtools intersect -wa -wb -a overlapping_genes.txt -b overlapping_fst.txt > overlapping_genes_fst.txt

##################
# Process "overlapping_genes_fst.txt"
# in R
#################

#Read annotation file

ovl_genes_fst<-read.delim("overlapping_genes_fst.txt", header = FALSE)
head(ovl_genes_fst)
names(ovl_genes_fst)

print(head(ovl_genes_fst, n=1))

# Set the column names
colnames(ovl_genes_fst) <- c(
  "chrom", "gene_start", "gene_end", "gene_width", "gene_id","chrom", "snp_start", "snp_end", "snp_width",
  "fst", "p_value")

print(head(ovl_genes_fst, n=5))

# subset ovl_genes_fst
# Create a subset of the "ovl_genes_fst" dataframe with selected columns
ovl_genes_fst_sub <- data.frame(
  gene_id = ovl_genes_fst$gene_id,
  fst = ovl_genes_fst$fst,
  p_value = ovl_genes_fst$p_value
)

print(head(ovl_genes_fst_sub, n=5))

# Drop Integer After Decimal Point in gene_id:
ovl_genes_fst_sub$gene_id <- sub("\\.\\d+", "", ovl_genes_fst_sub$gene_id)

# Order by fst (in descending order) and then by p_value (in ascending order) within each gene_id:
sorted_data <- ovl_genes_fst_sub %>%
  arrange(gene_id, desc(fst), p_value)

# Keep Unique gene_id with Highest fst but Lowest p_value:
unique_data <- sorted_data %>%
  group_by(gene_id) %>%
  filter(row_number() == 1)

print(head(unique_data, n=10))

########################

#### GO OUTPUT FILE
# Write both overlapping_genes and overlapping_fst as text files
write.table(unique_data, file = "GO_input.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Read the input file (gene expression values and significance)
#gene_list <- read.csv("output.csv")

gene_list <- read.table("GO_input.txt", header=TRUE) # output.csv
univ <- gene_list[, 3]
names(univ) <- gene_list[, 1]
univ <- univ[!is.na(univ)]

# Define a function to return TRUE/FALSE for p-values < 0.05
selection <- function(allScore) { return(allScore < 0.05) }

# Prepare topGO data
# Change the ontology to MF or BP or CC to view different ontology enrichments.
tGOdata <- new("topGOdata", description = "Simple session", ontology = "MF", geneSel = selection, allGenes = univ, nodeSize = 3, mapping = "org.At.tair.db", annot = annFUN.org)

# Run Gene Enrichment Analysis
results.ks <- runTest(tGOdata, algorithm = "elim", statistic = "ks")

# Generate a table of enriched GO terms
goEnrichment <- GenTable(tGOdata, KS = results.ks, orderBy = "KS", topNodes = 50)
goEnrichment$KS <- as.numeric(goEnrichment$KS)
goEnrichment <- goEnrichment[goEnrichment$KS < 0.05, ]
goEnrichment <- goEnrichment[, c("GO.ID", "Term", "KS")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep = ", ")

class(goEnrichment$KS)
# If the class is not numeric, you'll need to convert it to numeric:
goEnrichment$KS <- as.numeric(goEnrichment$KS)

# check for NAs
any(is.na(goEnrichment$KS))
head(goEnrichment[is.na(goEnrichment$KS), ])

# remove NAs
goEnrichment <- goEnrichment[!is.na(goEnrichment$KS), ]

# Save the enriched GO terms table to a CSV file
write.csv(goEnrichment, "sample_project.csv")

# Read the exported file
#defense <- read.csv("export_file.csv")

# Select the top N terms
#ntop <- 10
#ggdata <- defense[1:ntop, ]

# Select the top N terms
ntop <- 10
ggdata <- goEnrichment[1:ntop, ]

ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term))  # Fix order

# Plotting using ggplot2
gg1 <- ggplot(ggdata,
              aes(x = Term, y = -log10(KS), size = -log10(KS), fill = -log10(KS))) +
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5, 12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
  xlab('') + ylab('Enrichment score') +
  labs(
    title = 'Pathways Enriched in Genes with High Fst Values - MF',
    subtitle = 'Top 10 Terms Ordered by Kolmogorov-Smirnov p-Value',
    caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001'
  ) +
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "longdash", "solid"),colour = c("black", "black", "black"),size = c(0.5, 1.5, 3)) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, vjust = 1),
    axis.text.x = element_text(angle = 0, size = 12, hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12), axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'), axis.line = element_line(colour = 'black'),
    legend.key = element_blank(),legend.key.size = unit(1, "cm"),legend.text = element_text(size = 16),
    title = element_text(size = 12)) + coord_flip()

# Display and save the plot
print(gg1)

ggsave("GO_results_defense.png", plot = gg1, width = 11, height = 10, dpi = 300)

