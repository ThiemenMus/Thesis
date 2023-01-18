args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) < 1) {
  stop("At least one argument must be supplied (NSLOTS).n", call.=FALSE)
}

#Parameters##################### 
n_cores <- as.integer(args[1])
n_genes <- c(100, 400, 1200)
n_cells <- c(60, 200, 500, 1000)
sim_style <- c('clean', 'dirty')
################################

library(bc3net)
library(data.table)
library(parallel)
library(MASS)

mcinf <- function(i) {

  expr <- fread(file = sprintf('%s/%s', path_to_matrices, filenames[i]))
  
  # bc3net
  expr_matrix <- matrix(unlist(expr), nrow=nrow(expr)) #bc3 uses rows as genes, columns as plants; 
  #asks for matrix as input
  rownames(expr_matrix) <- c(1:nrow(expr)) #genes
  colnames(expr_matrix) <- c(1:ncol(expr)) #plants
  bc3net_grn <- bc3net(expr_matrix) 
  #expr must contain the names of the genes as rownames or an error is thrown, same for GENIE3
               
  print(i)
  return(bc3net_grn)
}

for (n in n_genes) {
  for (i in n_cells) {
    for (j in sim_style) {
  
      path_to_matrices <- sprintf('sergio_output/%sg/%sc/%s', n, i, j)
      filenames <- dir(path_to_matrices)
      output_name <- sprintf('%s%sg%s_bc3', substr(j, 1, 1), n/100, as.integer(i/100))
      
      inf_results <- mclapply(1:length(filenames), mcinf, mc.cores=n_cores)
      assign(output_name, inf_results)
      save(list=output_name, file = sprintf('inference_output/%s.RData', output_name))
    }
  }
}
