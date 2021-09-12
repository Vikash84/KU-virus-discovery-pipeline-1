args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else {
  require(data.table)
  require(dplyr)
  require(taxonomizr)
  data <- as.data.frame(fread(args[1]))
}

if (nrow(data)==0) {
  quit()
} else{
  uniq_taxid <- unique(data$"TAXID")
}

all_taxons <- as.data.frame(getTaxonomy(uniq_taxid, args[3]))
all_taxons <- tibble::rownames_to_column(all_taxons, "TAXID")
all_taxons$TAXID <- as.integer(all_taxons$TAXID)

taxon_mapped_data <- inner_join(data, all_taxons, by = "TAXID")

write.table(taxon_mapped_data, file = args[2], row.names=FALSE, col.names = TRUE, sep="\t", quote = FALSE)
