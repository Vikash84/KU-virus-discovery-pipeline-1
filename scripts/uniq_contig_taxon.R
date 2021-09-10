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
  uniq_data <- data[row.names(unique(data[,c("V1", "V3")])),]
}

uniq_data <- uniq_data %>% group_by(V3) %>% 
    mutate(mx = min(V12)) %>% 
    arrange(mx, V3) %>% 
    select(-mx)

write.table(uniq_data, file = args[1], row.names=FALSE, col.names = FALSE, sep="\t", quote = FALSE)
