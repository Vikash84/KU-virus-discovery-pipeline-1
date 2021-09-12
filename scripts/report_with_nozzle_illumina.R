# A script to generate the pipeline result report using Nozzle API.
# 
# Author: Kijin Kim (skkujin@gmail.com)
###############################################################

suppressPackageStartupMessages(library("argparser"))
require( "Nozzle.R1" );
require("dplyr");
require("data.table");
require("DT");

# --- Parsing functions ---

refMapFileToTable <- function(text_file) {
  ret <- read.table(text_file, header = TRUE, sep = "", strip.white = TRUE)
  ret <- select(ret, -c("N_COVERED_BASES", "AVG_BASEQ", "AVG_MAPQ"))
  return(ret)
}

AssembleSummaryFileToTable <- function(text_file) {
  ret <- read.table(text_file, header = FALSE, sep = ":", fill = TRUE, strip.white = TRUE)
  colnames(ret) <- c("Feature", "Value")
  return(ret)
}

blastFileToTable <- function(text_file) {
  ret <- read.table(text_file, header = TRUE, sep = "\t", strip.white = TRUE, quote="")
  ret$BITSCORE <- as.integer(ret$BITSCORE)
  return(ret)
}

uniq_ref_species <- function(blastTable) {
  ret <- blastTable[order(-blastTable$"BITSCORE"),]
  ret$species <- as.character(ret$species)
  ret <- ret[!duplicated(ret[ , "species"]), ]
  cnt <- as.data.frame(table(blastTable$"species"))
  colnames(cnt) <- c("species", "NUM_CONTIGS")

  ret <- inner_join(ret, cnt, by = "species")
  ret <- ret[c("NUM_CONTIGS", "species", "PER_IDENT", "EVALUE", "BITSCORE")]
  ret <- ret[order(-ret$"BITSCORE"),]
  colnames(ret) <- c("NUM_CONTIGS", "REF_SPECIES", "BEST_PER_IDENT", "BEST_EVALUE", "BEST_BITSCORE")
  return(ret)
}

# create parser object
p <- arg_parser("Report generating program")
  
# --- Script parameter parsing ---

p <- add_argument(p, "--prefix", help="Sample prefix")

p <- add_argument(p, "--multiqc_html", help="Multiqc html file path")

p <- add_argument(p, "--mapping_summary", help="Reference mapping summary file path")

p <- add_argument(p, "--classification_result_folder", help="The directory path where heatmap files exist")

p <- add_argument(p, "--assembly_summary", help="Contigs summary file path")

p <- add_argument(p, "--assembly_length_histogram", help="Contigs length histogram file path")

p <- add_argument(p, "--blastn_table", help="Blastn result file path")

p <- add_argument(p, "--megablast_table", help="Megablast result file path")

p <- add_argument(p, "--blastx_table", help="Blastx result file path")

p <- add_argument(p, "--blast_table_output_dir", help="Ouput directory path where blast html files will be stored")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parse_args(p)

prefix <- args$prefix


# qc data


if(is.na(args$multiqc_html)){
  qc_link <- paste(prefix, "/qc/multiqc_report.html", sep="")
} else{
  qc_link <- args$multiqc_html
}

# reference mapping data

if(is.na(args$mapping_summary)){
  ref_map_summary_file <- paste(prefix, "/mapping/", prefix, ".filtered_reference_mapping_collection.txt", sep="")
} else{
  ref_map_summary_file <- args$mapping_summary
}

# classification heatmap image files

if(is.na(args$classification_result_folder)){
  base_dir_path <- paste(prefix, "/classification", sep="")
} else{
  base_dir_path <- args$classification_result_folder
}

# de novo assembly data

if(is.na(args$assembly_summary)){
  assemble_summary_file <- paste(prefix, "/assembly/", prefix, ".contig_summary.txt", sep="")
  assemble_contig_length_histogram_file <- paste(prefix, "/assembly/", prefix, ".contigs_length_histogram.png", sep="")
} else{
  assemble_summary_file <- args$assembly_summary
}

# blast data

if(is.na(args$blastn_table)){
  blast_blastn_file <- paste(prefix, "/analysis/", prefix, ".blastn.txt", sep="")
} else{
  blast_blastn_file <- args$blastn_table
}

if(is.na(args$megablast_table)){
  blast_megablast_file <- paste(prefix, "/analysis/", prefix, ".megablast.txt", sep="")
} else{
  blast_megablast_file <- args$megablast_table
}

if(is.na(args$blastx_table)){
  blast_blastx_file <- paste(prefix, "/analysis/", prefix, ".blastx.txt", sep="")
} else{
  blast_blastx_file <- args$blastx_table
}




# --- Input Data ---


# reference mapping data
ref_map_summary_table <- refMapFileToTable(ref_map_summary_file)
#ref_map_virus_list <- 

# classification heatmap image files
classification_order_heatmap_file <- paste(base_dir_path, "/", prefix, "_order.svg", sep = "")
classification_family_heatmap_file <- paste(base_dir_path, "/", prefix, "_family.svg", sep = "")
classification_genus_heatmap_file <- paste(base_dir_path, "/", prefix, "_genus.svg", sep = "")
classification_species_heatmap_file <- paste(base_dir_path, "/", prefix, "_species.svg", sep = "")

# de novo assembly data
assemble_summary_table <- AssembleSummaryFileToTable(assemble_summary_file)

# blast data

blast_blastn_table <- blastFileToTable(blast_blastn_file)
blast_uniq_blastn_table <- uniq_ref_species(blast_blastn_table)
blast_megablast_table <- blastFileToTable(blast_megablast_file)
blast_uniq_megablast_table <- uniq_ref_species(blast_megablast_table)
blast_blastx_table <- blastFileToTable(blast_blastx_file)
blast_uniq_blastx_table <- uniq_ref_species(blast_blastx_table)

table_header_order <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "TAXID", "REF_ID", "REF_TITLE", "QUERY_ID", "REF_ID", "EVALUE", "BITSCORE", "PER_IDENT", "QUERY_LEN", "ALN_LEN", "MISMATCH", "GAPOPEN", "QSTART", "QEND", "REFSTART", "REFEND")

blast_blastn_table$TAXID <- as.factor(blast_blastn_table$TAXID)
blast_blastn_table <- blast_blastn_table[, table_header_order]
blastn_datatable <- datatable(blast_blastn_table, filter = 'top', options = list(pageLength = 20, autoWidth = TRUE))
blastn_html_link <- paste(prefix, ".full_blastn_table.html", sep = "")
saveWidget(blastn_datatable, blastn_html_link)

blast_megablast_table$TAXID <- as.factor(blast_megablast_table$TAXID)
blast_megablast_table <- blast_megablast_table[, table_header_order]
megablast_datatable <- datatable(blast_megablast_table, filter = 'top', options = list(pageLength = 20, autoWidth = TRUE))
megablast_html_link <- paste(prefix, ".full_megablast_table.html", sep = "")
saveWidget(megablast_datatable, megablast_html_link)

blast_blastx_table$TAXID <- as.factor(blast_blastx_table$TAXID)
blast_blastx_table <- blast_blastx_table[, table_header_order]
blastx_datatable <- datatable(blast_blastx_table, filter = 'top', options = list(pageLength = 20, autoWidth = TRUE))
blastx_html_link <- paste(prefix, ".full_blastx_table.html", sep = "")
saveWidget(blastx_datatable, blastx_html_link)

#===================================================================================================
# Report
#===================================================================================================

report <- newReport( paste(prefix, "KU Pipeline Analysis Report" ) );
report <- setReportSubTitle( report, "A report that showcases results generated from KU pipeline" );

# --- References ---

# create some references (even though they are included at the end of the report the need to be created
# first so they can be referenced in other elements)
#simpleCitation <- newCitation( authors="Nils Gehlenborg", title="Nozzle.R1 Package", year="2013", url="https://github.com/parklab/Nozzle" );
#webCitation <- newCitation( title="The Cancer Genome Atlas Website", url="http://tcga.cancer.gov/" );
#fullCitation <- newCitation( authors="Nils Gehlenborg", title="Nozzle: a report generation toolkit for data analysis pipelines", publication="Bioinformatics", issue="29", pages="1089-1091", year="2013", url="http://bioinformatics.oxfordjournals.org/content/29/8/1089" );

#report <- addToReferences( report, simpleCitation, webCitation, fullCitation );

# --- Results ---

# create objects

multiqc_link <- newParagraph( asLink( "MultiQC report", url=qc_link ));

ref_map_table <- newTable( ref_map_summary_table,
				"Reference mapping summary. Column follows bamcov format.");

classification_order_heatmap <- newFigure(classification_order_heatmap_file, "")
classification_family_heatmap <- newFigure(classification_family_heatmap_file, "")
classification_genus_heatmap <- newFigure(classification_genus_heatmap_file, "")
classification_species_heatmap <- newFigure(classification_species_heatmap_file, "")

assemble_table <- newTable( assemble_summary_table,
        "Summary statistics of assembled contigs.");

assemble_contig_length_histogram_figure <- newFigure( assemble_contig_length_histogram_file,
                                                      "")
						
blast_table_1 <- newTable( blast_uniq_blastn_table,
				"blastn result", file=blastn_html_link);

for ( i in 1:dim( blast_uniq_blastn_table )[1] )
{
  subtable_with_ref <- subset(blast_blastn_table, species == blast_uniq_blastn_table[i,"REF_SPECIES"])
  subtable_with_ref <- subtable_with_ref[c("QUERY_ID", "REF_ID", "PER_IDENT", "QUERY_LEN", "ALN_LEN",	"REFSTART", "REFEND", "EVALUE", "BITSCORE")]
  subtable_with_ref$REF_ID <- as.character(subtable_with_ref$REF_ID)
  result1 <- addTo( newResult( ""),
                    addTo( newSection( "Contigs assigned to ", blast_uniq_blastn_table[i,2] ), newTable(subtable_with_ref) ) );

  blast_table_1 <- addTo( blast_table_1, result1, row=i, column=1 );
}

blast_table_2 <- newTable( blast_uniq_megablast_table,
                           "megablast result", file=megablast_html_link);

for ( i in 1:dim( blast_uniq_megablast_table )[1] )
{
  subtable_with_ref <- subset(blast_megablast_table, species == blast_uniq_megablast_table[i,"REF_SPECIES"])
  subtable_with_ref <- subtable_with_ref[c("QUERY_ID", "REF_ID", "PER_IDENT", "QUERY_LEN", "ALN_LEN",	"REFSTART", "REFEND", "EVALUE", "BITSCORE")]
  subtable_with_ref$REF_ID <- as.character(subtable_with_ref$REF_ID)
  result1 <- addTo( newResult( ""),
                    addTo( newSection( "Contigs assigned to ", blast_uniq_megablast_table[i,2] ), newTable(subtable_with_ref) ) );
  
  blast_table_2 <- addTo( blast_table_2, result1, row=i, column=1 );
}

blast_table_3 <- newTable( blast_uniq_blastx_table,
                           "blastx result", file=blastx_html_link);
	
for ( i in 1:dim( blast_uniq_blastx_table )[1] )
{
  subtable_with_ref <- subset(blast_blastx_table, species == blast_uniq_blastx_table[i,"REF_SPECIES"])
  subtable_with_ref <- subtable_with_ref[c("QUERY_ID", "REF_ID", "PER_IDENT", "QUERY_LEN", "ALN_LEN",	"REFSTART", "REFEND", "EVALUE", "BITSCORE")]
  subtable_with_ref$REF_ID <- as.character(subtable_with_ref$REF_ID)
  result1 <- addTo( newResult( ""),
                    addTo( newSection( "Contigs assigned to ", blast_uniq_blastx_table[i,2] ), newTable(subtable_with_ref) ) );
  
  blast_table_3 <- addTo( blast_table_3, result1, row=i, column=1 );
}

report <- addToResults( report,
				addTo( newSubSection( "QC" ), multiqc_link ),
				addTo( newSubSection( "Reference Mapping" ), ref_map_table ), 
				addTo( newSubSection( "Classification heatmap" ), classification_order_heatmap, classification_family_heatmap, classification_genus_heatmap, classification_species_heatmap ),
				addTo( newSubSection( "Assembly Summary" ), assemble_table, assemble_contig_length_histogram_figure),
				addTo( newSubSection( "Blastn Result" ), blast_table_1 ), 
				addTo( newSubSection( "Megablast Result" ), blast_table_2 ), 
				addTo( newSubSection( "Blastx Result" ), blast_table_3 )); 
				
# set report maintainer information
report <- setMaintainerName( report, "Kijin Kim" );
report <- setMaintainerEmail( report, "skkujin@gmail.com" );
report <- setMaintainerAffiliation( report, "Saarland University" );

# set the copyright notice for this report
report <- setCopyright( report, owner="Kijin Kim", year=2021, statement="All rights reserved.", url="https://github.com/KijinKims/Long-read-Seq-Virome-Pipeline" ); 

# set contact information for error reports
report <- setContactInformation( report, email="skkujin@gmail.com", subject="Contact", message="Please email me if you have any question.", label="Contact" );

# --- HTML and RData file generation for report ---

writeReport( report, filename=paste(prefix,"_analysis_report", sep = "") )
