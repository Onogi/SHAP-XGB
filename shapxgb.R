shapxgb <- function(PhenotypeFile, GenotypeFile, Early_stopping_rounds = 3, Nrounds = 1000){
  
  #This is a script for conducting the SHAP-XGB algorithm for QTL mapping
  #The reference of SHAP-XGB is:
  #Ishibashi and Onogi (2024) biorxiv (https://doi.org/10.1101/2024.01.15.575690).
  
  #The followings are the formats
  ##Phenotype file
  #         TraitA  TraitB  TraitC
  #Line1    20.4    180     NA
  #Line2    12.0    145     10.3
  #...      ...     ...     ...
  #Line100  11.9    138     12.4
  
  ##=>The first column is the line ID
  ##=>The first row is the trait name
  ##=>NA is allowed. 
  ##=>Only quantitative values are allowed
  
  ##Genotype file
  #         Marker1 Marker2 Marker3
  #Line1    1       0       0
  #Line2    1       1       0
  #...      ...     ...     ...
  #Line100  -1      0       -1
  
  ##=>The first column is the line ID
  ##=>The first row is the marker name
  ##=>NA is allowed (imputed with averages)
  ##=>Genotypes are coded in the additive manner
  ###In the above example, 1:homo, 0:hetero, -1:homo
  ###But other coding such as [2:homo, 1:hetero: 0:homo] are also acceptable
  
  #Required packages (can be downloaded from CRAN at Sep. 6, 2024)
  ##xgboost (1.7.5.1)
  ##SHAPforxgboost (0.1.1)
  #Numbers are the versions used in the reference
  #These were run with R version 4.2.2.
  
  #Outputs
  ##Directories are created for each trait
  ##Each directory has a nested directory named as LocalImportance
  ##LocalImportance includes
  ###shap_values_average.csv: instance- (line/genotype-) wise importance for main effects 
  ###=>This includes N x P values
  ###shap_int_average_i.csv: instance-wise importance for interactions
  ###=>i takes a value from 1 to N. Each file includes P x (P + 1)/2 values.
  ##These values are the averages across bootstrap samples
  
  ##Global SHAP values are also output at each trait directory
  ###Shap_values_globalImportance.csv: global importance for main effects ("confound" effects)
  ###=>this includes P values
  ###Shap_int_globalimportance.csv: global importance for interaction effects
  ###=>this includes P * (P + 1)/2 values.
  
  #Read files
  ##Phenotypes
  Y_Ori <- as.matrix(read.csv(PhenotypeFile, header = TRUE, row.names = 1))
  ##Genotypes
  X_Ori <- as.matrix(read.csv(GenotypeFile, header = TRUE, row.names = 1))
  
  #Number of markers
  P <- ncol(X_Ori)
  #Number of traits
  E <- ncol(Y_Ori)
  
  #To extract and name interactions
  upper_tri <- upper.tri(matrix(0, P, P), diag = TRUE)
  intlabel <- matrix(colnames(X_Ori), P, P)
  temp <- matrix(colnames(X_Ori), P, P, byrow = TRUE)
  intlabel <- paste(intlabel[upper_tri], temp[upper_tri], sep = "-")
  rm(temp)
  
  #Missing genotypes are imputed with marker averages.
  for(i in 1:P){
    X_Ori[, i][is.na(X_Ori[, i])] <- mean(X_Ori[, i], na.rm = TRUE)
  }
  
  #Load the packages
  library(xgboost)
  library(SHAPforxgboost)
  
  #Number of bootstrap samples
  Nboot <- 100
  
  #Apply the SHAP-XGB algorithm to each trait
  for (target in 1:E) {
    
    #Create directory for the target trait
    folder_name <-paste0("SHAP-XGB.For.", colnames(Y_Ori)[target])
    dir.create(folder_name)
    dir.create(file.path(folder_name, "LocalImportance"))
    
    #Use lines with phenotypic values
    Y <- Y_Ori[complete.cases(Y_Ori[, target]), , drop = FALSE]
    X <- X_Ori[complete.cases(Y_Ori[, target]), , drop = FALSE]
    
    #Number of lines
    N <- nrow(X)
    
    #SHAP values
    shap_values_average <- matrix(0, nrow = N, ncol = P)
    
    #SHAP interaction values
    shap_int_average <- as.list(numeric(N))
    for(i in 1:N) shap_int_average[[i]] <- numeric(P * (P + 1)/2)
    
    #Bootstrap
    for (b in 1:Nboot){
      
      cat("trait", target, "bootstrap", b, "\n")
      
      #Divide data randomly
      Train.pop <- sort(sample(1:N, N, replace = TRUE))
      Val.pop <- c(1:N)[-unique(Train.pop)]
      
      #Train model
      Train <- xgb.DMatrix(data = X[Train.pop, ], label = Y[Train.pop, target])
      Val <- xgb.DMatrix(data = X[Val.pop, ], label = Y[Val.pop, target])
      model <- xgb.train(data = Train,
                         nrounds = Nrounds,
                         watchlist = list(train = Train, eval = Val),
                         early_stopping_rounds = Early_stopping_rounds,
                         verbose = 0)
      
      #Calculate SHAP values
      shap_values <- shap.values(xgb_model = model, X_train = Train)
      shap_values_average <- shap_values_average + as.matrix(shap_values$shap_score)
      
      shap_int <- shap.prep.interaction(xgb_mod = model, X_train = Train)
      shap_int <- shap_int[ , -c(P + 1), -c(P + 1)]#removing bias
      shap_int <- 2 * shap_int#doubled because only the upper triangle is used
      for (i in 1:N) {
        diag(shap_int[i, , ]) <- diag(shap_int[i, , ])/2#halve the diagonal elements
        shap_int_average[[i]] <- shap_int_average[[i]] + shap_int[i, , ][upper_tri]
      }
    }#b
    
    #Average bootstrap samples
    ##Main effects
    shap_values_average <- shap_values_average/Nboot
    rownames(shap_values_average) <- rownames(X)
    colnames(shap_values_average) <- colnames(X)
    write.csv(shap_values_average, 
              file = file.path(folder_name, "LocalImportance/shap_values_average.csv"))
    ##Interaction
    for (i in 1:N) {
      shap_int_average[[i]] <- shap_int_average[[i]]/Nboot
      names(shap_int_average[[i]]) <- intlabel
      write.csv(shap_int_average[[i]], 
                file = file.path(folder_name, "LocalImportance", paste0("shap_int_average_", i, ".csv")))
    }
    
    #Calculate global importance
    ##Main
    write.csv(colMeans(abs(shap_values_average)), 
              file = file.path(folder_name, "shap_values_globalImportance.csv"))
    
    ##Interaction
    shap_int_global <- numeric(P * (P + 1)/2)
    for(i in 1:N) shap_int_global <- shap_int_global + abs(shap_int_average[[i]])
    shap_int_global <- shap_int_global/N
    write.csv(shap_int_global, 
              file = file.path(folder_name, "shap_int_globalimportance.csv"))
  }#target
}