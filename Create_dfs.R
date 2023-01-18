library(tidyverse)
library(data.table)
library(gridExtra)
library(pcalg)
library(GeneNet)
library(bc3net)
library(GENIE3)

genes <- c("100g", "400g", "1200g")
obs <- c("60c", "200c", "500c", "1000c")

#######################################################################################################################
### Read in true GRN file from SERGIO ###
#########################################

read_grn <- function(interaction_file, n_genes) {
  
  grn_txt <- read.delim(interaction_file, header=FALSE)
  # 'input_file_taregts' file containing GRN structure and its parameters. for legend, cfr. https://github.com/PayamDiba/SERGIO
  
  
  ### Extracting edges from interaction file ###
  
  grn_edges <- data.frame(matrix(ncol=2,nrow=0, dimnames=list(NULL, c("reg", "target"))))
  count <- 0 
  #count number of regulators in GRN as affirmation that number of edges is correct
  
  for (i in seq(grn_txt[,1])) {
    #iterate over each target gene in the GRN
    gene_row <- as.integer(unlist(strsplit(grn_txt[i,], ','))) 
    #extract info of target gene i
    n_regs <- gene_row[2]
    count <- count + gene_row[2]
    
    for (j in seq(n_regs)) {  
      # iterate over each regulator of the target gene,
      # adding both to the final dataframe as an edge
      df <- data.frame(gene_row[2 + j] + 1, gene_row[1] + 1)
      # '+ 1' to adjust for the 0-index naming of the genes
      grn_edges <- rbind(grn_edges, setNames(df, names(grn_edges)))
    }
  }
  control1 <- count == nrow(grn_edges)
  #control whether all edges have been added to the data frame
  
  ### Control for cyclic edges, which should not be present.
  cyclic <- 0
  for (i in seq(nrow(grn_edges))) {
    edge_to_find <- data.frame(reg=grn_edges[[i, 2]], target=grn_edges[[i, 1]])
    if (nrow(merge(edge_to_find, grn_edges)) > 0) {
      count <- count + 1
    }
  }
  control2 <- cyclic == 0
  
  ### Control for presence of n genes in the GRN
  both <- list()
  for (i in grn_edges["reg"]) {both <- append(both, i)}
  for (i in grn_edges["target"]) {both <- append(both, i)}
  
  control3 <- length(unique(both)) == n_genes
  
  cat('all edges extracted:', control1, '\n') 
  cat('no cyclic edges:', control2, '\n')
  cat(sprintf('%s genes present in GRN:', n_genes), control3)
  
  grn_edges <- as_tibble(grn_edges) %>%
    mutate(reg = as.character(reg), target = as.character(target))
  return(as_tibble(grn_edges))

}

### Reading folder of R objects ###
###################################

### load folder of RData ###
load_rdata <- function(folder) {
  filenames <- dir(folder)
  for (i in filenames) {
    load(sprintf('%s/%s', folder, i), envir=.GlobalEnv)
  }
}

### load folder of rds objects ###
load_rds <- function(folder) {
  
  filenames <- dir(folder)
  for (i in filenames) {
    results <- str_replace(i, ".rds", "")
    sub <- readRDS(sprintf('%s/%s', folder, i))
    assign(results, sub)
  }
}

### inference results to GRN ###
################################

### PC result to GRN
pc_to_grn <- function(inf_out) {
  
  grn <- tribble(~reg, ~target, ~p, ~directed)
  p_vals <- inf_out@pMax
  adj <- summary(inf_out)
  
  for (i in 1:nrow(adj)) {
    targets_of_i <- which(adj[i,-(1:i)] == 1) + i # + 1 to offset the lowered indices because of subsetting
    regs_of_i <- which(adj[-(1:i),i] == 1) + i
    sect <- intersect(targets_of_i, regs_of_i)
    
    if (length(sect) > 0) {
      targets_of_i <- targets_of_i[!(targets_of_i %in% sect)]
      regs_of_i <- regs_of_i[!(regs_of_i %in% sect)]
    }
    
    for (target in targets_of_i) {
      grn <- add_row(grn, reg=i, target=target, p=p_vals[i, target], directed=TRUE)
    }
    for (reg in regs_of_i) {
      grn <- add_row(grn, reg=reg, target=i, p=p_vals[reg, i], directed=TRUE)
    }
    for (node in sect) {
      grn <- add_row(grn, reg=i, target=node, p=p_vals[i, node], directed=FALSE)
    }
  }
  grn <- mutate(grn, method='pc', reg = as.character(reg), target = as.character(target))
  return(grn)
}

### bc3 result to GRN
bc3_to_grn <- function(inf_out) {
  
  grn <- as_edgelist(inf_out)
  
  if (length(grn) > 0) {
    grn <- as_tibble(grn) %>%
    mutate(weight=E(inf_out)$weight, directed=FALSE, method='bc3') %>%
    arrange(desc(weight)) %>%
    rename(reg = V1, target = V2)
    return(grn)
  } else {
    return(NA)
  }
}

### gn result to GRN
gn_to_grn <- function(inf_out) {
  
  if (nrow(inf_out) == 0) {
    return(NA)
  }
  
  grn <- as_tibble(inf_out) %>%
    dplyr::select(node1, node2, pval) %>%
    rename(reg = node1, target = node2) %>% #does genenet return directed edges? Not by default!
    mutate(directed = FALSE, method='gn', reg = as.character(reg), target = as.character(target))
  
  return(grn)
}

### genie result to GRN: only keeping the highest weight edge of both options.
genie_to_grn <- function(inf_out, size=FALSE) {
  if (size) {
    inf_out <- inf_out[1:size,]
  }
  grn <- inf_out %>%
    mutate(edge_key = case_when(reg < target ~paste(reg, target, sep='_'),
                                reg > target ~paste(target, reg, sep='_'))) %>%
    arrange(edge_key) %>%
    distinct(edge_key, .keep_all=TRUE) %>%
    dplyr::select(reg, target, weight) %>%
    arrange(desc(weight))%>%
    mutate(directed = TRUE, method='genie', reg=as.character(reg), target=as.character(target))
  
  return(grn)
}
    
### inference results to GRNs ###
inf_to_grns <- function(inf_results, size=FALSE) {
  
  results_class <- class(inf_results[[1]])
  if (length(results_class) > 1) {
    results_class <- results_class[1]
  }
                          
  n_sims <- length(inf_results)
  grn_list <- vector("list", n_sims)
  
  ### PC
  if (results_class == 'pcAlgo') {
    print('PC')
    for (i in 1:n_sims) {
      grn_list[[i]] <- pc_to_grn(inf_results[[i]])
      if (i %in% seq(0, n_sims, 10)) {
        print(paste(i, "% has been processed"))
      }
    }
  }
  ### bc3
  else if (results_class == 'igraph') {
    print('bc3')
    for (i in 1:n_sims) {
      grn_list[[i]] <- bc3_to_grn(inf_results[[i]])
      if (i %in% seq(0, n_sims, 10)) {
        print(paste(i, "% has been processed"))
      }
    }
  }
  ### genenet
  else if (results_class == 'data.frame') {
    print('genenet')
    for (i in 1:n_sims) {
      grn_list[[i]] <- gn_to_grn(inf_results[[i]])
      if (i %in% seq(0, n_sims, 10)) {
        print(paste(i, "% has been processed"))
      }
    }
  }
  ### Genie : delete directed edges with lowest weight
  else if (is(inf_results[[1]], "tbl_df")) {
    print('genie')
    for (i in 1:n_sims) {
      grn_list[[i]] <- genie_to_grn(inf_results[[i]], size=size)
      if (i %in% seq(0, n_sims, 10)) {
        print(paste(i, "% has been processed"))
      }
    }
  }
  else {
    print('No valid object class')
  }
  return(grn_list)
}

### Summarizing statistics of interest per list of GRNs ###
##################################################

### Precision and recall ###

summarise_grn <- function(fake_grn, true_grn, true_grn_inverse=FALSE) {
  
  if (is(true_grn_inverse, 'logical')) {
    true_grn_inverse <- mutate(true_grn, reg2 = target, target2 = reg) %>%
      select(reg2, target2) %>%
      rename(reg=reg2, target=target2)
  }
  
  if (is(fake_grn, 'logical')) {
    
    stats <- tribble(~prec, ~rec, ~TP, ~PP, 0, 0, 0, 0)
    
  } else {
    
    TP <- nrow(inner_join(fake_grn, true_grn)) + nrow(inner_join(fake_grn, true_grn_inverse))
    PP <- nrow(fake_grn) # Proposed positives = TP + FP
    prec <- TP / PP
    rec <- TP / nrow(true_grn) # nrow(true_grn) = TP + FN (all positives)
    
    stats <- tribble(~prec, ~rec, ~TP, ~PP, prec, rec, TP, PP)
  }
  
  return(stats) 
  
}

### list of GRNs to statistics ###

grns_to_stats <- function(grn_list, true_grn) {
  
  stats <- tribble(~precision, ~recall, ~TP, ~PP)
  true_grn_inverse <- mutate(true_grn, reg2 = target, target2 = reg) %>%
    select(reg2, target2) %>%
    rename(reg=reg2, target=target2)
  
  if (is(grn_list[[1]], 'logical')) {
    i <- 2
    method <- NA
    while(is.na(method) & i <= length(grn_list)) {
      if (is(grn_list[[i]], 'logical')) {
        i <- i + 1
      } else {
        method <- grn_list[[i]]$method[1]
      }
    }
    if (i > length(grn_list)) {
      print('No results found')
      return(NA)
    }
  } else {
    method <- grn_list[[1]]$method[1]
  }
  
  
  for (i in 1:length(grn_list)) {
    
    fake_grn <- grn_list[[i]]
    stats <- rbind(stats, summarise_grn(fake_grn, true_grn, true_grn_inverse=true_grn_inverse))
    
    cat(sprintf('%s/%s\r', i, length(grn_list)))
  }
  
  stats <- mutate(stats, method = method)
  return(stats)
}


### ranked precision of list of GRNs ###
########################################
grns_to_rankprec <- function(grn_list, true_grn) {
  grn_rank <- tribble(~reg, ~target, ~weight, ~directed, ~method, ~rank)
  ### add rank and concatenate all grns
  for (i in 1:length(grn_list)) {
    if (!(is(grn_list[[i]], 'logical'))) {
      grn_list[[i]] <- mutate(grn_list[[i]], rank = 1:nrow(grn_list[[i]]))
      grn_rank <- rbind(grn_rank, grn_list[[i]])
    }
  }
  grn_rank <- arrange(grn_rank, rank)
  
  ### add 1 for TP 0 for FP
  true_regs <- true_grn$reg; fake_regs <- grn_rank$reg
  true_targets <- true_grn$target; fake_targets <- grn_rank$target
  TPs <- c()
  for (i in 1:nrow(true_grn)) {
    curr_reg <- true_regs[i]
    curr_target <- true_targets[i]
    TPs <- c(TPs, which(fake_regs == curr_reg & fake_targets == curr_target))
    TPs <- c(TPs, which(fake_regs == curr_target & fake_targets == curr_reg))
  }
  
  grn_rank <- mutate(grn_rank, TP = 0)
  for (i in TPs) {
    grn_rank[i, 'TP'] = 1
  }
  grn_rankprec <- group_by(grn_rank, rank) %>%
    summarise(rank_precision=sum(TP)/100, rank_precision_adj=sum(TP)/(sum(rank) / mean(rank)), count = sum(rank) / mean(rank))
  
  return(list(grn_rank, grn_rankprec))
}

##################################################################################################################################################################################################

true_grn1 <- read_grn("Interaction100g.txt", 100)
true_grn4 <- read_grn("Interaction400g.txt", 400)
true_grn12 <- read_grn("Interaction1200g.txt", 1200)

for (g in genes) {
  for (o in obs) {
    
    print(sprintf('%s Genes, %s Observations', g, o))
    path <- sprintf('inference_output/%s/%s', g, o)
    
    filenames <- dir(path)
    filenames <- filenames[grepl('rds', filenames)]
    print(filenames)
    
    if (g == "100g") {
      true_grn <- true_grn1 
    } else if (g == "400g") {
      true_grn <- true_grn4
    } else if (g == "1200g") {
      true_grn <- true_grn12
    }
    
    for (f in filenames) {
      f_object <- readRDS(sprintf('%s/%s', path, f))
      
      if (grepl('genie', f)) {                         
      #to avoid huge tibbles
        f_grns <- inf_to_grns(f_object, size=nrow(true_grn))
      } else {
        f_grns <- inf_to_grns(f_object)
      }
      
      f_stats <- grns_to_stats(f_grns, true_grn)
      f_rankedprec <- grns_to_rankprec(f_grns, true_grn)
      f_merge <- list(f_stats, f_rankedprec)
      
      output_name <- sprintf('%s_stats', unlist(strsplit(f, ".rds")))
      saveRDS(f_merge, file=sprintf('%s/%s.rds', path, output_name))
      
    }
  }
}
      
      
      
      
      
