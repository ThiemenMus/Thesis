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

library(GENIE3)
library(data.table)
library(tidyverse)
library(parallel)
library(MASS)

mcinf <- function(i) {
  
  print(sprintf('genie cores: %s', parallel::detectCores()))
  
  expr <- fread(file = sprintf('%s/%s', path_to_matrices, filenames[i]))
  
  # GENIE3
  expr_matrix <- matrix(unlist(expr), nrow=nrow(expr)) #genie uses rows as genes, columns as plants; ~bc3
  rownames(expr_matrix) <- c(1:nrow(expr)) #genes
  colnames(expr_matrix) <- c(1:ncol(expr)) #plants
  adj_mat <- GENIE3(
    expr_matrix,
    regulators = NULL,
    targets = NULL,
    treeMethod = "RF",
    K = "sqrt",
    nTrees = 1000,
    nCores = 1, #otherwise there are 8*8 = 64 cores being used?
    returnMatrix = TRUE,
    verbose = FALSE
  )
  genie_grn <- getLinkList(adj_mat) %>%
        as_tibble() %>%
        rename(reg = regulatoryGene, target = targetGene)
               
  print(i)
  return(genie_grn)
}

for (n in n_genes) {
  for (i in n_cells) {
    for (j in sim_style) {
  
      path_to_matrices <- sprintf('sergio_output/%sg/%sc/%s', n, i, j)
      filenames <- dir(path_to_matrices)
      output_name <- sprintf('%s%sg%s_genie', substr(j, 1, 1), n/100, as.integer(i/100))
      
      inf_results <- mclapply(1:length(filenames), mcinf, mc.cores=n_cores)
      assign(output_name, inf_results)
      save(list=output_name, file = sprintf('inference_output/%s.RData', output_name))
    }
  }
}
