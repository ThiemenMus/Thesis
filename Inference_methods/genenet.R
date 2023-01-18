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

library(GeneNet)
library(data.table)
library(parallel)
library(MASS)

mcinf <- function(i) {

  expr <- fread(file = sprintf('%s/%s', path_to_matrices, filenames[i]))

  # GeneNet
  expr_t <- t(expr) #GeneNet uses rows as plants, columns as genes ~PC
  expr_pcor <- ggm.estimate.pcor(expr_t) #estimate partial correlation matrix from expression matrix
  expr_edges <- network.test.edges(expr_pcor) #convert pcor to edge list
  genenet_grn <- extract.network(expr_edges, cutoff.ggm=0.95)
               
  print(i)
  return(genenet_grn)
}

for (n in n_genes) {
  for (i in n_cells) {
    for (j in sim_style) {
  
      path_to_matrices <- sprintf('sergio_output/%sg/%sc/%s', n, i, j)
      filenames <- dir(path_to_matrices)
      output_name <- sprintf('%s%sg%s_gn', substr(j, 1, 1), n/100, as.integer(i/100))
      
      inf_results <- mclapply(1:length(filenames), mcinf, mc.cores=n_cores)
      assign(output_name, inf_results)
      save(list=output_name, file = sprintf('inference_output/%s.RData', output_name))
    }
  }
}
