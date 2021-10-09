args = commandArgs(trailingOnly=TRUE)

require(MetaComp)

assignment <- load_edge_assignment(args[1], type = "kraken")

plot_merged_assignment(assignment, "class", row_limit = 10, min_row_abundance = 1.0e-10, plot_title = args[2], filename = paste(args[2], "_", "class.svg", sep = ""))
plot_merged_assignment(assignment, "order", row_limit = 10, min_row_abundance = 1.0e-10, plot_title = args[2], filename = paste(args[2], "_", "order.svg", sep = ""))
plot_merged_assignment(assignment, "family", row_limit = 10, min_row_abundance = 1.0e-10, plot_title = args[2], filename = paste(args[2], "_", "family.svg", sep = ""))
plot_merged_assignment(assignment, "genus", row_limit = 10, min_row_abundance = 1.0e-10, plot_title = args[2], filename = paste(args[2], "_", "genus.svg", sep = ""))
plot_merged_assignment(assignment, "species", row_limit = 10, min_row_abundance = 1.0e-10, plot_title = args[2], filename = paste(args[2], "_", "species.svg", sep = ""))

