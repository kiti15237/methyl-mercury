library(dplyr)
library(readr)
library(phyloseq)
library(reshape2)
library(ggplot2)
library(glmnet)




loadPhyloseq <- function(dir){
  ps <- readRDS(file.path(dir, "data/ps.RDS"))
  ps_dds <- readRDS(file.path(dir, "data/ps_dds.RDS"))
  return(list(ps, ps_dds))
}

cleanSampleData <- function(ps){
  df <- sample_data(ps_dds)
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
  ps_filt <- prune_samples(colSums(otu_table(ps_filt)) > 0, ps_filt) #Drop samples that are empty
  
  names <- apply(tax_table(ps_filt), 1, function(phy_names){
    return(paste(tax_table(ps_filt)[1,], collapse = "_"))
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


regression <- function(y, input, method ){
  #input will have all input variables and the target
  x <- colnames(input)[colnames(input) != y]
  
  if(method == "cvglmnet"){
    target <- input[,y]
    input <- select(input, -y)
    fit <- cv.glmnet(as.matrix(input), target,family = "gaussian")
    train_preds = predict(fit, as.matrix(input), s = fit$lambda[length(fit$lambda)])
    fit$mse = sum((train_preds - target)^2)
    print(fit$mse)
  }
  if(method == "basic"){
    form <- getFormula(x, y)
    fit <- lm(form, input)
    fit$rval <- summary(fit)$adj.r.squared
  }
  if(method == "glmnet"){
    input <- input[,x]
    target <- input[,y]
    fit <- glmnet(x = as.matrix(input), target, family = "gaussian")
    num_lamb <- length(fit$lambda)
    rsquared <- fit$dev.ratio[num_lamb]
    print(rsquared)
  }
  
  return(fit)
}



regression_randomNoise <- function(y, df, otu, method, num_vars, permutations = 10){
  other <- getOther(y)
  df_use <- select(df, -c(other))
  
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
    
    x <- colnames(df_use)
    x <- x[x != y]
    form <- getFormula(x, y)
    df_stand <- as.data.frame(apply(df_use, 2, function(var) return( (var - mean(var, na.rm = T)) / sd(var, na.rm = T))))
    
    if(method == "glm"){
      fit <- glm(form, df_stand, family = "gaussian")
      print(fit$aic)
    }
    if(method == "basic"){
      fit <- lm(form, df_stand)
      rsquared <- summary(fit)$adj.r.squared
    }
    if(method == "glmnet"){
      input <- df_stand[,x]
      target <- df_stand[,y]
      fit <- glmnet(x = as.matrix(input), target, family = "gaussian")
      num_lamb <- length(fit$lambda)
      rsquared <- fit$dev.ratio[num_lamb]
      print(rsquared)
    }
    if(method == "cvglmnet"){
      input <- df_stand[,x]
      target <- df_stand[,y]
      fit <- cv.glmnet(as.matrix(input), target,family = "gaussian")
      train_preds = predict(fit, as.matrix(input), s = fit$lambda[length(fit$lambda) - 5])
      fit$mse = sum((train_preds - target)^2)
      rsquared <- fit$mse
    }
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
  }else if(ylab == 'mse'){
    p <- p + ylab("Mean Squared Error")
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
    return(coefficients(model))
  }

}


getRvalues <- function(method, model, s='lambda.min'){
  if(method=="glmnet"){
    num_lamb <- 30
    return(model$dev.ratio[num_lamb])
  }
  if(method == "basic"){
    return(summary(model)$adj.r.squared)
  }
  #if(method == "cvglmnet"){
  #  return(model$mse)
  #}
  if(method == "cvglmnet"){
    return(model$glmnet.fit$dev.ratio[which(model$lambda == model[[s]])])
  }
}

runLinearRegression <- function(y, ps, method, permutations, title = ""){
  
  #Get data together
  basic_input <- getBasicInput(ps, y = y)
  merc_input <- getMercInput(ps, y = y)
  otu_input <- getOtuInput(ps, y = y)
  all_input <- getAllInput(ps, y = y)
  
  #Run regression on each dataset
  set.seed(1)
  basic <- regression(y, basic_input, method)
  basic_merc <- regression(y, merc_input, method)
  basic_otu <- regression(y, otu_input, method)
  basic_all <- regression(y, all_input, method)
  
  #Make predictions using variables of random noise
  df <- sample_data(ps)
  otu_mat <- otu_table(ps)
  if(nrow(otu_mat) != nsamples(ps)){
    otu_mat <- t(otu_mat)
  }
  random1 = regression_randomNoise(y, df, otu_mat, method, num_vars = 3,  permutations = permutations)
  
  random2 = regression_randomNoise(y, df, otu_mat, method, num_vars = ntaxa(ps), permutations = permutations)
  
  random3 = regression_randomNoise(y, df, otu_mat, method, num_vars = ntaxa(ps) + 3, permutations)
  
  #Plot results
  ylab = "var_exp"
  if(method == "cvglmnet"){
    ylab= "mse"
  }
  plotVariationExplained(data.frame(
    basic = getRvalues(method, basic),
    basic_merc = getRvalues(method, basic_merc), 
    basic_micro = getRvalues(method, basic_otu),
    basic_micro_merc = getRvalues(method, basic_all),
    rand1 = random1, rand1 = random2, rand3 = random3),
    title = title, ylab = ylab)

  
  print("Basic Model: ")
  print(getSigCoefs(basic, method))
  
  print("Basic + Merc Model: ")
  print(getSigCoefs(basic_merc, method))
  
  print("Basic + Micro Model: ")
  print(getSigCoefs(basic_otu, method))
  
  print("Basic + Micro + Merc Model: ")
  print(getSigCoefs(basic_all, method))
  return(list(basic, basic_merc, basic_otu, basic_all))
}
