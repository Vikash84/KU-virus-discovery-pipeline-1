# A script to generate the pipeline result report using R markdown.
# 
# Author: Kijin Kim (skkujin@gmail.com)
###############################################################

suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("rmarkdown"))
suppressPackageStartupMessages(library("knitr"))

# create parser object
p <- arg_parser("Nanopore analysis report generating program")
  
# --- Script parameter parsing ---

p <- add_argument(p, "--prefix", help="Sample prefix")

p <- add_argument(p, "--nanoplot_html", help="NanoPlot html file path")

p <- add_argument(p, "--mapping_summary", help="Reference mapping summary file path")

p <- add_argument(p, "--classification_result_folder", help="The directory path where heatmap files exist")

p <- add_argument(p, "--classification_kraken_html", help="Kraken result html file path")

p <- add_argument(p, "--classification_kaiju_html", help="Kaiju result html file path")

p <- add_argument(p, "--contigs_len", help="Contigs length file path")

p <- add_argument(p, "--blastn_table", help="Blastn result file path")

p <- add_argument(p, "--megablast_table", help="Megablast result file path")

p <- add_argument(p, "--blastx_table", help="Blastx result file path")

p <- add_argument(p, "--blast_table_output_dir", help="Ouput directory path where blast html files will be stored")

# get command line options,
# if options not found on command line then set defaults automatically generated with prefix, 
args <- parse_args(p)
prefix <- args$prefix

# qc data
if(is.na(args$nanoplot_html)){
  qc_link <- paste("qc/", prefix, "_NanoPlot-report.html", sep="")
} else{
  qc_link <- args$nanoplot_html
}

# reference mapping data

if(is.na(args$mapping_summary)){
  ref_map_summary_file <- paste(prefix, "/mapping/", prefix, ".reference_mapping.txt", sep="")
} else{
  ref_map_summary_file <- args$mapping_summary
}

# classification heatmap image files

if(is.na(args$classification_result_folder)){
  classification_base_dir_path <- "classification"
} else{
  classification_base_dir_path <- args$classification_result_folder
}

if(is.na(args$classification_kraken_html)){
  classification_kraken_link <- paste(classification_base_dir_path, "/", prefix, ".kraken.html", sep = "")
} else{
  classification_kraken_link <- args$classification_kraken_html
}

if(is.na(args$classification_kaiju_html)){
  classification_kaiju_link <- paste(classification_base_dir_path, "/", prefix, ".kaiju.html", sep = "")
} else{
  classification_kaiju_link <- args$classification_kaiju_html
}

# de novo assembly data
if(is.na(args$assembly_summary)){
  contigs_len_file <- paste(prefix, "/assembly/", prefix, ".contigs.len", sep="")
} else{
  contigs_len_file <- args$contigs_len
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


rmd_template_file <- paste(dirname(rstudioapi::getSourceEditorContext()$path), "/6_report_template_nanopore.Rmd", sep="")
render(rmd_template_file,
      output_file = paste(prefix,"/",prefix,"_analysis_report", sep = "")
      params=list(prefix=prefix,
                  qc_link_=qc_link, 
                  ref_map_summary_file=ref_map_summary_file, 
                  classification_base_dir_path=classification_base_dir_path, 
                  classification_kraken_link=classification_kraken_link, 
                  classification_kaiju_link=classification_kaiju_link,
                  contigs_len_file=contigs_len_file,
                  blast_blastn_file=blast_blastn_file,
                  blast_megablast_file=blast_megablast_file,
                  blast_blastx_file=blast_blastx_file))