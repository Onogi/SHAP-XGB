#This is a script for conducting SHAP-XGB.
#The reference of SHAP-XGB is:
#Ishibashi and Onogi (2024) biorxiv (https://doi.org/10.1101/2024.01.15.575690).

#When applying SHAP-XGB to users' data, replace the following file names with their own csv file names.
PhenotypeFile <- "DaysToFlowering.csv"
GenotypeFile <- "Genotype.csv"

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
###Shap_values_average.csv: instance- (line/genotype-) wise importance for main effects
###=>This includes N x P values
###Shap_int_average_i.csv: instance-wise importance for interactions
###=>i takes a value from 1 to P. Each file includes N x P values.
##These values are the averages across bootstrap samples

##Global SHAP values are also output at each trait directory
###Shap_values_globalImportance.csv: global importance for main effects ("compouond" effects)
###=>this includes P x 1 values
###Shap_int_globalimportance.csv: global importance for interaction effects
###=>this consists of a P x P matrix.


#Read files#####################################################################
#Phenotypes
Y_Ori <- as.matrix(read.csv(PhenotypeFile, header = TRUE, row.names = 1))

#Genotypes
X_Ori <- as.matrix(read.csv(GenotypeFile, header = TRUE, row.names = 1))

#Number of markers
P <- ncol(X_Ori)

#Number of traits
E <- ncol(Y_Ori)

#Missing genotypes are imputed with marker averages.
for(i in 1:P){
  X_Ori[, i][is.na(X_Ori[, i])] <- mean(X_Ori[, i], na.rm = TRUE)
}


#Applying SHAP-XGB##############################################################
library(xgboost)
library(SHAPforxgboost)

#Number of bootstrap samples
Nboot <- 100

#Apply SHAP-XGB to each trait
for (target in 1:E) {
  
  #Create directory for the target trait
  folder_name <-paste0("SHAP-XGB.at.", colnames(Y_Ori)[target])
  dir.create(folder_name)
  dir.create(file.path(folder_name, "LocalImportance"))
  
  #Use lines with phenotypic values
  Y <- Y_Ori[complete.cases(Y_Ori[, target]), ]
  X <- X_Ori[complete.cases(Y_Ori[, target]), ]
  
  #Number of lines
  N <- nrow(X)
  
  #Shap values for each marker and line in each bootstrap sample
  shap_values_results <- matrix(0, nrow = N * P, ncol = Nboot)
  colnames(shap_values_results) <- 1:Nboot
  
  #Shap interaction values for each marker combination and line in each bootstrap sample
  ##The object becomes quite large when P is high.
  shap_int_results <- matrix(0, nrow = P * N * P, ncol = Nboot)
  colnames(shap_int_results) <- 1:Nboot
  
  #Bootstrap
  for (i in 1:Nboot){
    
    cat(target, i, "\n")
    
    #Divide data randomly
    Train.pop <- sort(sample(1:N, N, replace = TRUE))
    Val.pop <- c(1:N)[-unique(Train.pop)]
    
    #Train model
    Train <- xgb.DMatrix(data = X[Train.pop, ], label = Y[Train.pop, target])
    Val <- xgb.DMatrix(data = X[Val.pop, ], label = Y[Val.pop, target])
    model <- xgb.train(data = Train,
                       nrounds = 1000,
                       watchlist = list(train = Train, eval = Val),
                       early_stopping_rounds = 3,
                       verbose = 0)
    
    #Calculate Shap values
    shap_values <- shap.values(xgb_model = model, X_train = Train)
    shap_int <- shap.prep.interaction(xgb_mod = model, X_train = Train)
    
    #Save values
    shap_values_results[ ,i] <- unlist(shap_values$shap_score)
    shap_int_results[, i] <- as.vector(shap_int[ , -P, -P])
  }#i
  
  #Average bootstrap samples
  ##Main effects
  shap_values_average <- rowMeans(shap_values_results)
  write.csv(shap_values_average, file = file.path(folder_name, "LocalImportance/Shap_values_average.csv"))
  ##Interaction (divide the results into P chunks)
  for (i in 1:P) {
    start_row <- (i - 1) * (N * P) + 1
    end_row <- i * (N * P)
    chunk <- shap_int_results[start_row:end_row, ]
    shap_int_average <- rowMeans(chunk)
    write.csv(shap_int_average, file = file.path(folder_name, "LocalImportance", paste0("Shap_int_average_", i, ".csv")))
  }
  
}#target

#remove large objects and clean up
rm(shap_values_results, shap_int_results, shap_values_average, chunk, shap_int_average)
gc();gc();gc();gc()


#Calculate global importance####################################################
#Main
SHAP <- matrix(NA, nrow = P, ncol = E)
colnames(SHAP) <- colnames(Y_Ori)
for(target in 1:E){

  N <- sum(!is.na(Y_Ori[, target]))
  Env <- colnames(Y_Ori)[target]
  
  result <- unlist(read.csv(paste0("SHAP-XGB.at.", Env, "/LocalImportance/Shap_values_average.csv"), 
                            header = TRUE, row.names = 1))
  result <- matrix(result, nrow = N, ncol = P)
  SHAP[, target] <- apply(abs(result), 2, mean)
  write.csv(SHAP[, target], file = file.path(folder_name, "Shap_values_globalImportance.csv"))
}

#Interaction
SHAP.int <- as.list(numeric(E))
names(SHAP.int) <- colnames(Y_Ori)
for(target in 1:E){

  N <- sum(!is.na(Y_Ori[, target]))
  Env <- colnames(Y_Ori)[target]
  SHAP.int[[target]] <- matrix(NA, P, P)
  
  for(j in 1:(P - 1)){
    
    result <- as.matrix(read.csv(paste0("SHAP-XGB.at.", Env, "/LocalImportance/Shap_int_average_", j, ".csv"), 
                                 header = TRUE, row.names = 1))
    #Two notes:
    ##the number of columns of result is 1.
    ##j needs not to be P because interactions with Pth marker are already included in the files until P.
    for(k in (j+1):P){
      #SHAP interaction values are doubled
      SHAP.int[[target]][k, j] <- mean(abs(result[((k - 1) * N + 1):(k * N), ])) * 2
    }
  }
  write.csv(SHAP.int[[target]], file = file.path(folder_name, "Shap_int_globalimportance.csv")) 
}

