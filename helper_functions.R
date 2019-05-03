library(dplyr)
library(readr)
library(phyloseq)
library(reshape2)
library(ggplot2)
library(glmnet)
library(vegan)



loadPhyloseq <- function(dir){
  ps <- readRDS(file.path(dir, "data/ps_w_tree.RDS"))
  #ps_dds <- readRDS(file.path(dir, "data/ps_dds.RDS"))
  return(ps)
}

cleanSampleData <- function(ps){
  df <- sample_data(ps)
  df <- select(df, -c("ForwardFastqFile" , "ReverseFastqFile", "Description", "mdi_median_1","pdi_median_1"))
  df <- na.omit(df)
  sample_data(ps) <- df
  return(ps)
}

getCommonOtus <- function(ps, numSamples){
  otu <- otu_table(ps)
  if(nrow(otu) < ncol(otu)){
    # Flip so that we have a lot of rows rather than a lot of columns
    #Taxa by Samples
    otu <- t(otu)
  }
  keep <- apply(otu, 1, function(taxa_dist) return(sum(taxa_dist > 0) > numSamples))# Keep taxa present in > 25% of samples
  ps_filt <- prune_taxa(keep, ps)
 # ps_filt <- prune_samples(colSums(otu_table(ps_filt)) > 0, ps_filt) #Drop samples that are empty
  
  names <- apply(tax_table(ps_filt), 1, function(phy_names){
    return(paste(phy_names[[5]], phy_names[[6]], sep = "_"))
  })
  
  taxa_names(ps_filt) <- make.unique(names)
  
  return(ps_filt)
}

writeFasta <- function(seqs, file_name){
  headers <- paste(">seq" , seq(1, length(seqs)), sep = "") 
  #Sequences must all be the same length
  seqs <- sapply(seqs, function(x) return(strtrim(x, 440)))
  fasta = paste(headers, seqs, sep = "\n")
  write(fasta, file_name)
}



standardize <- function(df){
  #Standardize
  df_stand <- as.data.frame(apply(df, 2, function(var){
    if(range(var)[1] == 0 & range(var)[2] == 1){
      return(var)
    }
    else{
      return( (var - mean(var, na.rm = T)) / sd(var, na.rm = T))
    }
  }))
  return(df_stand)
}

getBasicInput <- function(ps, y){
  merc_vars <- c("maternal_hair_thg_median", "child_hair_thg_median", "child_hair_thg")
  df <- sample_data(ps)
  df <- select(df, -c(merc_vars))
  other <- getOther(y)
  df <- select(df, -c(other))
  return(standardize(df))
}

getMercInput <- function(ps, y){
  df <- sample_data(ps)
  other <- getOther(y)
  df <- select(df, -c(other))
  return(standardize(df))
}

getOtuInput <- function(ps, y){
  merc_vars <- c("maternal_hair_thg_median", "child_hair_thg_median", "child_hair_thg")
  df <- sample_data(ps)
  df <- select(df, -c(merc_vars))
  
  otu <- otu_table(ps)
  if(taxa_are_rows(ps)){
    otu <- t(otu) #We will concatenate with df so samples should be rows
  }
  input <- cbind(df, otu)
  colnames(input) <- c(colnames(df), colnames(otu))
  other <- getOther(y)
  input <- select(input, -c(other))
  
  return(standardize(input))
}

getAllInput <- function(ps, y){
  df <- sample_data(ps)
  otu <- otu_table(ps)
  if(taxa_are_rows(ps)){
    otu <- t(otu) #We will concatenate with df so samples should be rows
  }
  input <- cbind(df, otu)
  colnames(input) <- c(colnames(df), colnames(otu))
  other <- getOther(y)
  input <- select(input, -c(other))
  return(standardize(input))
}


getOther <- function(y){
  other = ""
  if(y == "bayley36_mdi_adj"){
    other = "bayley36_pdi_adj"
  }
  if(y == "bayley36_pdi_adj"){
    other = "bayley36_mdi_adj"
  }
  return(other)
}

getFormula <- function(x, y){
  #x=Variables to input  #y=target
  vars <- paste(x, collapse = "+")
  ind <- paste(y, " ~")
  form <- formula(paste(ind, vars, sep = ""))
  return(form)
}


regression <- function(y, input, method, reg_strength){
  #input will have all input variables and the target
  x <- colnames(input)[colnames(input) != y]
  
  if(method == "cvglmnet"){
    target <- input[,y]
    input <- select(input, -y)
    fit <- cv.glmnet(as.matrix(input), target, family = "gaussian", alpha = 0.95)
    lambda = fit$lambda[length(fit$lambda)]
    train_preds = predict(fit, as.matrix(input), s = lambda)
    print(train_preds)
    fit$mse = sum((train_preds - target)^2)
    
    #Save lambda of choice in model
    #Deal with the regularization strength issue
    if(reg_strength %in% c("lambda.min", "lambda.1se")){
      fit$lambda_chose <- reg_strength
    }
    else{
      #Full reg strength is 1, no reg strength is 0
      lambda = fit$lambda[round(length(fit$lambda) * (1-reg_strength))]
      fit$lambda_chose <- lambda
    }
  }
  if(method == "basic"){
    form <- getFormula(x, y)
    fit <- lm(form, input)
    fit$rval <- summary(fit)$adj.r.squared
  }
  
  return(fit)
}



regression_randomNoise <- function(y, df, otu, method, num_vars, permutations = 10, reg_strength){
  #other <- getOther(y)
  #df_use <- select(df, -c(other))
  df_use <- df
  
  set.seed(0)
  rsquared_cum <- c()
  ranges = apply(otu, 2, function(taxa) return(list(range(taxa))))
  
  for(perm in seq(permutations)){
    numRands <- 0
    #add random noise column as control
    #Range of random noise column will match range of a randomly selected otu
    for(i in seq(1, num_vars)){
      tmp <- paste("random", i, sep = "")
      range = sample(ranges, 1)
      range = range[[1]][[1]]
      df_use[ , tmp] <- runif(nrow(df_use), min = range[1], max = range[2])
      numRands = numRands + 1
    }
    
    df_stand <- standardize(df_use)
    
    model <- regression(y = y, input = df_stand, method = method, reg_strength)
    rsquared <- getRvalues(method = method, model = model, reg_strength = reg_strength)
    rsquared_cum <- c(rsquared_cum, rsquared)
  }
  
  return(mean(rsquared_cum))
}


plotVariationExplained <- function(df, title = "", ylab = 'var_exp'){
  tmp <- df
  tmp_melt <- melt(tmp)
  colnames(tmp_melt) <- c("Model", "value")
  tmp_melt$value <- round(tmp_melt$value, 4)
  p <- ggplot(tmp_melt, aes(x = Model, y = value, colour = Model, fill = Model)) 
  p <- p + geom_bar(stat = "identity")
  p <- p + geom_text(aes(label=value), vjust = -0.3, size = 10)
  if(ylab == 'var_exp'){
    p <- p + ylab("Variance explained (R^2)")
  }else if(ylab == 'dev.ratio'){
    p <- p + ylab("Fraction of null deviance explained")
  }
  
  p <- p + ggtitle(title)
  print(p)
}

getSigCoefs <- function(model, method){
  if(method == "basic"){
    coefs <- summary(model)$coefficients
    return(coefs[coefs[,4]  < .07, , drop = F])
  }
  if(method == "cvglmnet"){
    return(coef(model, s = model$lambda_chose)[as.numeric(coef(model, s= model$lambda_chose)) != 0, ])
  }

}


getRvalues <- function(method, model, reg_strength){
  
  if(method=="glmnet"){
    num_lamb <- 30
    return(model$dev.ratio[num_lamb])
  }
  if(method == "basic"){
    return(summary(model)$adj.r.squared)
  }

  if(method == "cvglmnet"){
    
    #Deal with the regularization strength issue
    if(reg_strength %in% c("lambda.min", "lambda.1se")){
      return(model$glmnet.fit$dev.ratio[which(model$lambda == model[[reg_strength]])])
    }
    #If lambda.min and lambda.1se are not allowing us to fit a model at all, we lower the regularization manually
    else{
      #Full reg strength is 1, no reg strength is 0
      lambda = model$lambda[round(length(model$lambda) * (1-reg_strength))]
      dev_exp <- model$glmnet.fit$dev.ratio[which(model$lambda == lambda)]
      print(dev_exp)
      return(dev_exp) #% deviance explained
    }

  }
}

runLinearRegression <- function(y, ps, method, permutations, reg_strength, title = ""){
  
  #Get data together
  basic_input <- getBasicInput(ps, y = y)
  merc_input <- getMercInput(ps, y = y)
  otu_input <- getOtuInput(ps, y = y)
  all_input <- getAllInput(ps, y = y)
  
  #Run regression on each dataset
  set.seed(1)
  basic <- regression(y, basic_input, method, reg_strength)
  basic_merc <- regression(y, merc_input, method, reg_strength)
  basic_otu <- regression(y, otu_input, method, reg_strength)
  basic_all <- regression(y, all_input, method, reg_strength)
  
  #Make predictions using variables of random noise
  df <- sample_data(ps)
  otu_mat <- otu_table(ps)
  if(nrow(otu_mat) != nsamples(ps)){
    otu_mat <- t(otu_mat)
  }
  
  #Deal with regularization parameter

  random1 = regression_randomNoise(y, basic_input, otu_mat, method, num_vars = 3,  permutations = permutations, reg_strength)
  
  random2 = regression_randomNoise(y, basic_input, otu_mat, method, num_vars = ntaxa(ps), permutations = permutations, reg_strength = reg_strength)
  
  random3 = regression_randomNoise(y, basic_input, otu_mat, method, num_vars = ntaxa(ps) + 3, permutations, reg_strength = reg_strength)
  
  #Plot results
  ylab = "var_exp"
  if(method == "cvglmnet"){
    ylab= "dev.ratio"
  }
  
  
  plotVariationExplained(data.frame(
    basic = getRvalues(method, basic, reg_strength),
    basic_merc = getRvalues(method, basic_merc, reg_strength), 
    basic_micro = getRvalues(method, basic_otu, reg_strength),
    basic_micro_merc = getRvalues(method, basic_all, reg_strength),
    rand1 = random1, rand2 = random2, rand3 = random3),
    title = title, ylab = ylab)

  
  print("Basic Model: ")
  sig_vars_basic <- getSigCoefs(basic, method)
  print(sig_vars_basic)
  
  print("Basic + Merc Model: ")
  sig_vars_merc <- getSigCoefs(basic_merc, method)
  print(sig_vars_merc)
  
  print("Basic + Micro Model: ")
  sig_vars_micro <- getSigCoefs(basic_otu, method)
  print(sig_vars_micro)
  
  print("Basic + Micro + Merc Model: ")
  sig_vars_all <- getSigCoefs(basic_all, method)
  print(sig_vars_all)
  return(list(basic, basic_merc, basic_otu, basic_all))
}

runPermanova <- function(ps, method){
  otu <- otu_table(ps)
  if(taxa_are_rows(ps)){
    otu <- t(otu) #want samples as rows
  }
  dist_bray <- vegdist(otu, method = method)
  vars <- paste(sample_variables(ps), collapse = "+")
  form <-  formula(paste("dist_bray ~", vars))
  perm = adonis2(form, data = as.data.frame(sample_data(ps)), method = "bray", perm = 9999)
}

recordSignificance <- function(ps, model, method, column_name){
  sig_vars_micro <- getSigCoefs(model, method = method)
  
  tmp_taxa <- as.data.frame(tax_table(ps))
  
  if(!(column_name %in% colnames(tmp_taxa))){
    tmp_taxa[ , column_name] <- as.numeric(rep(0, ntaxa(ps)))
  }
  sig_taxa <- names(sig_vars_micro)[names(sig_vars_micro) %in% taxa_names(ps)]
  tmp_taxa[ , column_name] <- as.numeric(as.character(tmp_taxa[ , column_name]))
  tmp_taxa[sig_taxa, column_name] <- ifelse(sig_vars_micro[sig_taxa] > 0, "pos", "neg")
  tax_table(ps) <- tax_table(as.matrix(tmp_taxa))
  return(ps)
}


plotCCA <- function(ps, dist, color, shape, type ){
  vars <- paste(colnames(sample_data(ps)), collapse = "+") #Coordinate system
  ind <- paste("~")
  form <- formula(paste(ind, vars, sep = ""))
  #Place samples based on distances between their microbiomes
  cca <- ordinate(physeq = ps, method = "CCA", distance = dist, formula = form)
  cap_plot <- plot_ordination(ps, 
                              ordination = cca,
                              type = type, 
                              color= color, 
                              shape = shape, 
                              axes = c(1,2)) + geom_point(size=5)
 
  
  arrowmat <- vegan::scores(cca, display = "bp")
  # Add labels, make a data.frame
  arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
  # Define the arrow aesthetic mapping
  arrow_map <- aes(xend = CCA1 * 4, 
                   yend = CCA2 * 4, 
                   x = 0, 
                   y = 0, 
                   shape = NULL, 
                   color = NULL, 
                   label = labels)
  
  label_map <- aes(x = 4.5 * CCA1, 
                   y = 4.5 * CCA2, 
                   shape = NULL, 
                   color = NULL, 
                   label = labels)
  
  arrowhead = arrow(length = unit(0.02, "npc"))
  
  # Make a new graphic
  cap_plot <- cap_plot + 
    geom_segment(
      mapping = arrow_map, 
      size = .5, 
      data = arrowdf, 
      color = "gray", 
      arrow = arrowhead
    ) + 
    geom_text_repel(
      mapping = label_map, 
      size = 4,  
      data = arrowdf, 
      show.legend = FALSE
    )
  return(cap_plot)
}

library(ggrepel)
