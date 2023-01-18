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

library(pcalg)
library(data.table)
library(parallel)
library(MASS)

mcinf <- function(i) {

  expr <- fread(file = sprintf('%s/%s', path_to_matrices, filenames[i]))
  
  # PC
  expr_t <- t(expr) #PC uses rows as plants, columns as genes
  suffStat <- list(C = cor(expr_t), n = nrow(expr_t))
  varNames <- as.character(1:ncol(expr_t))
  pc_grn <- pc(suffStat, indepTest = gaussCItest, labels = varNames,
               alpha = 0.01, numCores = 1) #prevent 8*8 = 64 cores being used?
               
  print(i)
  return(pc_grn)
}

for (n in n_genes) {
  for (i in n_cells) {
    for (j in sim_style) {
  
      path_to_matrices <- sprintf('sergio_output/%sg/%sc/%s', n, i, j)
      filenames <- dir(path_to_matrices)
      output_name <- sprintf('%s%sg%s_pc', substr(j, 1, 1), n/100, as.integer(i/100))
      
      inf_results <- mclapply(1:length(filenames), mcinf, mc.cores=n_cores)
      assign(output_name, inf_results)
      save(list=output_name, file = sprintf('inference_output/%s.RData', output_name))
    }
  }
}



