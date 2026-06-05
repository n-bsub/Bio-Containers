# Set personal R library path (for HPC use)

# Get the input folder from the environment (set in your sbatch with --env INPUT_DIR=...)
input_dir <- Sys.getenv("INPUT_DIR", unset = "")
if (!nzchar(input_dir) || !dir.exists(input_dir)) {
  stop("INPUT_DIR is not set or does not exist: ", input_dir)
}

# (optional) tiny helper so paths are shorter to write
f <- function(...) file.path(input_dir, ...)

# Load packages (avoid tidyverse metapackage to reduce dependency issues)
library(readxl)
library(writexl)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(data.table)
library(glue)
library(ggplot2)
library(DESeq2)
library(vegan)
library(RColorBrewer)
library(WGCNA)
library(edgeR)
library(limma)
library(matrixStats)

enableWGCNAThreads()

# Note: Avoid loading UniProt.ws due to conflict with dplyr::select


Blast_GO_wide = readRDS(f("Blast_GO_wide_CHEMM.RDS"))

# If using all genes and not just fungi, uncomment below
#Blast_GO_wide_fungi = Blast_GO_wide

# If using all genes and not just fungi, comment below
Blast_GO_wide_fungi = Blast_GO_wide %>%
  filter(grepl('fungi', Org.XP))

CPM = read.table(f("RSEM.gene.counts.matrix"), sep = "\t", 
                 stringsAsFactors = FALSE, header = TRUE, row.names = NULL,
                 check.names = FALSE) 
colnames(CPM)[1]="#gene_id"
TMM = read.table(f("RSEM.gene.TMM.EXPR.matrix"), sep = "\t", 
                 stringsAsFactors = FALSE, header = TRUE, row.names = NULL,
                 check.names = FALSE) 
colnames(TMM)[1]="#gene_id"

#Left-join 
CPM_fungi = left_join(Blast_GO_wide_fungi, CPM, by = "#gene_id")
TMM_fungi = left_join(Blast_GO_wide_fungi, TMM, by = "#gene_id")

# Check to see if Blast_GO_wide_fungi's Trinity gene ids are a subset of those in CPM and TMM
all(unique(Blast_GO_wide_fungi$`#gene_id`) %in% unique(CPM_fungi$`#gene_id`))
all(unique(Blast_GO_wide_fungi$`#gene_id`) %in% unique(TMM_fungi$`#gene_id`))

#Order #gene_id column and remove duplicates
CPM_fungi = CPM_fungi[order(CPM_fungi$'#gene_id'), ]
CPM_fungi_U = CPM_fungi[!duplicated(CPM_fungi$`#gene_id`),]

TMM_fungi = TMM_fungi[order(TMM_fungi$'#gene_id'), ]
TMM_fungi_U = TMM_fungi[!duplicated(TMM_fungi$`#gene_id`),]

# Dust only so exclude Carpet samples from CPM_fungi_U
CPM_fungi_U = CPM_fungi_U[,-(17:25)]

# Prepare data for samples-RH or metadata information

CPM_fungi_meta = CPM_fungi_U[,-(1:7)]
rownames(CPM_fungi_meta) = CPM_fungi_U[,1] 

CPM_fungi_meta1 = CPM_fungi_meta %>%
  t() %>%
  data.frame() %>%
  mutate(ERH = case_when(str_detect(rownames(.), "50") ~ "50%",
                         str_detect(rownames(.), "85") ~ "85%",
                         str_detect(rownames(.), "95") ~ "95%")) %>%
  relocate(ERH, .before = everything()) 

# Create samples-RH or metadata file
Metadata = CPM_fungi_meta1 %>%
  select(c(ERH)) %>%
  data.frame() %>%
  mutate(Samples = rownames(.)) %>%
  relocate(Samples, .before = everything())

# For CPM, we want 12 or more samples to have expression values of 1 and over
# Change 12 or more samples to 6 or more sample since we are doing Dust only here
# Testing: Changing from 6 to 3
CPM_fungi_U_fil = CPM_fungi_U %>%
  filter(rowSums(across(8:16, ~ . >=1)) >= 3)

# Prepare voom input from the counts matrix
counts_matrix <- CPM_fungi_U_fil[, 8:ncol(CPM_fungi_U_fil)]  # raw counts only
rownames(counts_matrix) <- CPM_fungi_U_fil[, 1]              # gene IDs as rownames

# Match metadata rownames to column names of counts
Metadata <- Metadata[match(colnames(counts_matrix), Metadata$Samples), ]
rownames(Metadata) <- Metadata$Samples

# Build design matrix (use ~1 for unsupervised WGCNA, or ~ERH + Site for supervised limma)
design <- model.matrix(~1, data = Metadata)

# Create DGEList object
dge <- DGEList(counts = counts_matrix)
dge <- calcNormFactors(dge, method = "TMM")

# Apply voom transformation
voom_out <- voom(dge, design = design, plot = TRUE)
voom_data <- t(voom_out$E)  # transpose to samples x genes


voom_final2 = t(voom_data)
#Check metadata rows columns with Voom data
all(rownames(Metadata) %in% colnames(voom_final2))
all(rownames(Metadata) == colnames(voom_final2))


# After voom and collapsing by Sprot.Blast.XP
library(matrixStats)

voom_final2 = voom_final2 %>%
  t()
#norm_data <- voom_final2[,colSums(voom_final2) >= 0] # use 2 or higher as filter
norm_data <- voom_final2

## For testing pupose
#norm_data = t(voom_df[,-19])


# check that there is not too much missing data (flagging OTUs and samples with too many missing values)
# should return "[1] TRUE"
gsg = goodSamplesGenes(norm_data, verbose = 3);
gsg$allOK

# Plot voom
pdf("boxplot_6.pdf")
boxplot(t(norm_data), main='voom-normalized log2-CPM (filtered)', las=2)
dev.off()

# clustering SAMPLES (not OTUs) by their similarities in OTU abundance, 
clustTree1 = hclust(dist(norm_data), method = "average")

# Plot the tree: Open a graphic output window of size 12 by 9 inches
sizeGrWindow(12,9)
#pdf(file = "Gene_module_clusters_loc.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))

pdf("samples_cluster.pdf")
plot(clustTree1, main = "Gene module clusters", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()


# we apply a series of power transformations on the edges of the network till we fit the scale free topology model of the network. 
powers = c(c(1:10), seq(from = 12, to=50, by=2))
# Call the network topology analysis function, this takes 30ish seconds
sft = pickSoftThreshold(norm_data, powerVector = powers, 
                        verbose = 5, 
                        networkType = "signed hybrid")
# Plot the results:
sizeGrWindow(9, 5)
#pdf("Powers.pdf",width=12,height=8)
par(mfrow = c(1,2));
cex1 = 0.9

pdf("sft_plot.pdf")

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",
     main = paste("Scale independence"));
abline(h=0.8, col = "blue")
abline(h=0.7, col = "blue")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

power_norm_data<-9


# Now we define the modules from the network
# we break the data into multiple blocks that are then rejoined (blockwiseModules)
# if we want to try with one block we can attempt a maxBlockSize = ## possibly up to the size of the entire dataset for better modules, 
# but this will be extremely slow for large datasets. You also have control over the minmodulesize.
# The Topological Overlap Measure (TOM) makes networks less sensitive to spurious connections or to connections missing due to random noise. 
# The central idea of TOM is to count the direct connection strengths as well as connection strengths "mediated" by shared neighbors.
# PAM = Partitioning Around Medoids. When pamRespectsDendro = FALSE, dendrogram is ignored and only dissimilarity information is
# used to assign any previously unassigned objects to clusters detected in earlier steps
net = blockwiseModules(norm_data, power = power_norm_data,
                       corType="pearson",  #maxPOutliers = 0.10,
                       maxBlockSize = 50000,
                       networkType = "signed hybrid", 
                       TOMType = "signed", minModuleSize = 100,
                       reassignThreshold = 0, 
                       #deepSplit = 2, 
                       mergeCutHeight = 0.1,
                       numericLabels = FALSE, 
                       pamRespectsDendro = FALSE,
                       saveTOMs = FALSE, 
                       #replaceMissingAdjacencies = TRUE,
                       #saveTOMFileBase = "MB", 
                       #minKMEtoStay = 0.4, 
                       verbose = 3)

#plotDendroAndColors(
# net$dendrograms[[1]],
#net$colors,
#"Module colors",
#dendroLabels = FALSE,
#hang = 0.03,
#addGuide = TRUE, guideHang = 0.05
#)



for (block in 1:length(net$dendrograms)) {
  plotDendroAndColors(
    net$dendrograms[[block]],
    labels2colors(net$colors[net$blockGenes[[block]]]),
    groupLabels = "Module colors",
    main = paste("Gene dendrogram and module colors in block", block),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )
}



table(net$colors)

#Module-gene mapping 
module_colors <- net$colors
gene_names <- colnames(norm_data)
module_gene_table <- data.frame(Gene = gene_names, Module = module_colors)
module_gene_table = left_join(module_gene_table, Blast_GO_wide_fungi[,-2], by = c("Gene" = "#gene_id"))

#Relocate and write table
module_gene_table = module_gene_table %>%
  relocate(Module, .after = everything()) %>%
  distinct(Gene, Module, .keep_all = TRUE)

write.table(module_gene_table, "gene_module_mapping.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Split modules

# Get unique module names
modules <- unique(module_gene_table$Module)

# Loop through each module and write to a file
for (mod in modules) {
  df_mod <- module_gene_table %>%
    filter(Module == mod) %>%
    select(Gene, Module)
  
  # Write to a file named after the module (e.g., blue_module.txt)
  write.table(df_mod,
              file = file.path("modules", paste0(mod, "_module.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

#for halla with voom data

voom_halla = t(voom_data)
voom_halla_rownames = cbind(Gene = rownames(voom_halla), voom_halla)
voom_halla_rownames_df <- as.data.frame(voom_halla_rownames)

Voom_Sprot_map = unique(CPM_fungi_U_fil[,c(1,3)])
colnames(Voom_Sprot_map)[1] = "Gene"

voom_halla_sprot = left_join(voom_halla_rownames_df, Voom_Sprot_map, by = "Gene")
voom_halla_sprot = voom_halla_sprot %>%
  relocate(Sprot.Blast.XP, .after = Gene)
voom_halla_sprot = voom_halla_sprot[,-1]

write.table(voom_halla_sprot, "voom_halla_sprot.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(colnames(voom_data), "wgcna_input_gene_list_background.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# 1. Calculate module eigengenes
MEs <- moduleEigengenes(norm_data, colors = module_colors)$eigengenes

# 2. Save module eigengenes
saveRDS(MEs, file = "MEs.RDS")

# 3. Save metadata (RH, Site, etc.)
saveRDS(Metadata, file = "Metadata.RDS")

# 4. Save norm_data
saveRDS(norm_data, "norm_data.RDS")

# 5. Save module assignments (gene-module mapping)
saveRDS(module_colors, file = "module_colors.RDS")

head(names(module_colors))

# ---------------- Part 2: Module–Trait correlation (robust alignment) ----------------

# datExpr comes from your pipeline; ensure it's numeric
datExpr <- as.matrix(norm_data)        # WGCNA expects samples x genes (rows = samples)
storage.mode(datExpr) <- "double"

# Load metadata from INPUT_DIR and standardize
datTraits <- data.table::fread(f("Metadata_wgcna.txt"), data.table = FALSE)

# Clean header names (strip potential BOM) and whitespace
names(datTraits) <- trimws(sub("^\ufeff", "", names(datTraits)))

# Pick the sample ID column: accept either "Samples" or "Sample"
id_col <- if ("Samples" %in% names(datTraits)) "Samples" else
  if ("Sample"  %in% names(datTraits)) "Sample"  else
    stop("Metadata must contain a 'Samples' or 'Sample' column. Found: ",
         paste(names(datTraits), collapse = ", "))

# Build rownames safely (trim, check empties/dups)
ids <- trimws(as.character(datTraits[[id_col]]))
if (anyNA(ids) || any(ids == "")) {
  bad <- which(is.na(ids) | ids == "")
  stop("Metadata has empty/NA sample IDs in rows: ", paste(head(bad, 10), collapse = ", "))
}
if (anyDuplicated(ids)) {
  dups <- unique(ids[duplicated(ids)])
  stop("Metadata has duplicated sample IDs (first few): ", paste(head(dups, 10), collapse = ", "))
}
rownames(datTraits) <- ids
datTraits[[id_col]] <- NULL

# Decide whether samples are rows or columns in datExpr; transpose if needed
rn <- rownames(datExpr); if (is.null(rn)) rn <- character(0)
cn <- colnames(datExpr); if (is.null(cn)) cn <- character(0)
over_r <- sum(rownames(datTraits) %in% rn)
over_c <- sum(rownames(datTraits) %in% cn)

if (over_r == 0 && over_c == 0) {
  cat("Example datExpr rownames: ", paste(head(rn, 5), collapse = ", "), "\n")
  cat("Example datExpr colnames: ", paste(head(cn, 5), collapse = ", "), "\n")
  cat("Example trait IDs: ", paste(head(rownames(datTraits), 5), collapse = ", "), "\n")
  stop("No overlap between sample IDs in metadata and datExpr row/col names.")
}

if (over_c > over_r) {
  # Samples are in columns; transpose so samples become rows
  datExpr <- t(datExpr)
}

# Final trim on rownames
rownames(datExpr) <- trimws(rownames(datExpr))

# Align by intersection (keeps only shared samples and puts them in the same order)
common <- intersect(rownames(datExpr), rownames(datTraits))
if (length(common) < 2) {
  cat("Matched samples: ", length(common), "\n")
  cat("In traits not in expr (first 10): ",
      paste(head(setdiff(rownames(datTraits), rownames(datExpr)), 10), collapse = ", "), "\n")
  cat("In expr not in traits (first 10): ",
      paste(head(setdiff(rownames(datExpr), rownames(datTraits)), 10), collapse = ", "), "\n")
  stop("Too few overlapping samples between expression and metadata.")
}

datExpr   <- datExpr[common, , drop = FALSE]
datTraits <- datTraits[common, , drop = FALSE]
stopifnot(identical(rownames(datExpr), rownames(datTraits)))

# Sanity: module_colors must match # of genes (columns of datExpr)
if (length(module_colors) != ncol(datExpr)) {
  stop("Length of module_colors (", length(module_colors),
       ") does not equal number of genes/columns in datExpr (", ncol(datExpr), ").")
}

# Compute MEs and correlate with traits
MEs0 <- moduleEigengenes(datExpr, colors = module_colors)$eigengenes
MEs  <- orderMEs(MEs0)

moduleTraitCor    <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

# Create text matrix for labeling
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# Plot heatmap
pdf("Module_Trait_Correlation_Heatmap_Dust.pdf", width = 10, height = 8)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = "Module–Trait Relationships")
dev.off()

# Save correlations and p-values to CSV
write.csv(moduleTraitCor, "Module_Trait_Correlations.csv")
write.csv(moduleTraitPvalue, "Module_Trait_Pvalues.csv")

