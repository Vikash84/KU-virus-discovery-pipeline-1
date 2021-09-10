args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  require(data.table)
  require(dplyr)
  data <- as.data.frame(fread(args[1]))
}

if (nrow(data)==0) {
  quit()
} else{

filtered_data <- data %>%
    filter(score > 0.95, pvalue < 0.05) 
}

if (nrow(filtered_data)>0){
    write.table(filtered_data, file = args[1], row.names=FALSE, col.names = TRUE, sep="\t", quote = FALSE)
} else {
    file.remove(args[1])
}
