args = commandArgs(trailingOnly=TRUE)

require(MetaComp)

assignment <- load_edge_assignment(args[1], type = "kraken")

plot_edge_assignment(assignment, "family", "Project", "column", "family.svg")
plot_edge_assignment(assignment, "phylum", "Project", "column", "phylum.svg")
plot_edge_assignment(assignment, "species", "Project", "column", "species.svg")