library(readr)
library(stringr)
library(writexl)

# Set directory containing your enrichment result files
results_dir <- "Enrichment_results"

# Find all *.GOseq.enriched files
enriched_files <- list.files(results_dir, pattern = "\\.GOseq\\.enriched$", full.names = TRUE)

# Create named list to store dataframes by module name
goseq_list <- list()

for (file in enriched_files) {
  # Get module name from filename, e.g., "blue" from "blue_module.txt.GOseq.enriched"
  filename <- basename(file)
  module_name <- str_remove(filename, "_module\\.txt\\.GOseq\\.enriched")
  
  # Read the enrichment result file
  df <- read_tsv(file, col_types = cols(.default = "c"))
  
  # Store in list under module name
  goseq_list[[module_name]] <- df
}

# Truncate all long strings to avoid Excel limit
MAX_CHAR <- 32000
goseq_list <- lapply(goseq_list, function(df) {
  df[] <- lapply(df, function(col) {
    if (is.character(col)) {
      substr(col, 1, MAX_CHAR)
    } else {
      col
    }
  })
  return(df)
})

# Write to a single Excel workbook
write_xlsx(goseq_list, path = "Combined_GOseq_Enrichment.xlsx")
