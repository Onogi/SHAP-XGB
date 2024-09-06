#Here, rice heading date data is used as an example.
#Onogi A et al. (2016) Theor Appl Genet 129:805–817 (https://doi.org/10.1007/s00122-016-2667-5).
#The data is redistributed here.

#Drawing figures for the main and interaction effects are also demonstrated.


#Read files#####################################################################
#Phenotypes
Y_Ori <- as.matrix(read.csv("DaysToFlowering.csv", header = TRUE, row.names = 1))
#Environment names of Y_Ori contain errors and are fixed here
colnames(Y_Ori)[7] <- "Tsukuba2008E"
colnames(Y_Ori)[8] <- "Tsukuba2008L"

#Genotypes
X_Ori <- as.matrix(read.csv("Genotype.csv", header = TRUE, row.names = 1))

#Number of markers
P <- ncol(X_Ori)

#Number of traits
E <- ncol(Y_Ori)

#Missing genotypes are imputed with marker averages.
for(i in 1:P){
  X_Ori[, i][is.na(X_Ori[, i])] <- mean(X_Ori[, i], na.rm = TRUE)
}

#Map of markers
Map <- read.csv("Map.csv", header = TRUE, row.names = 1)

#Markers nearest to Major QTLs
MajorQTL <- c(92, 60, 112, 108)
names(MajorQTL) <- c("HD1", "HD6", "DTH8", "HD2")


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


#Plot main and interaction global importance####################################
#Colors for different chromosomes
Col <- NULL
for(chr in 1:12) Col <- c(Col, rep(chr%%2 + 1, sum(Map$Chr == chr)))
Col[Col == 1] <- 4

#Main
tiff("SHAPvaluesForMaineffects.tiff", unit = "cm", res = 400, height = 18, width = 24)
par(mfcol = c(3, 3))
par(mar = c(1, 2.5, 1, 0.3))
par(mgp = c(1.5, 0.5, 0))
for(i in c(5, 1, 7, 8, 9, 2, 4, 6, 3)){#change the order of environments along with latitudes.
  plot(SHAP[, i], 
       col = Col, 
       main = colnames(Y_Ori)[i], 
       pch = 1, 
       xlab = "", 
       ylab = "Global Shap values", 
       xaxt = "n")
  for(j in MajorQTL) abline(v = j, lty=3, col="gray70")
}
dev.off()


#Heat map for interactions
tiff("SHAPvaluesForInteractions.tiff", unit = "cm", res = 400, height = 18, width = 18)
par(mfrow = c(3, 3))
par(mar=c(0.5,0.5,0.5,0.5))
for(i in c(5, 1, 7, 8, 9, 2, 4, 6, 3)){
  image(SHAP.int[[i]], xaxt ="n", yaxt = "n")
  text(0.4, 0.8, colnames(Y_Ori)[i])
  for(j in MajorQTL) abline(v = j/162, lty=3, col="gray85", lwd=0.5)
  for(j in MajorQTL) abline(h = j/162, lty=3, col="gray85", lwd=0.5)
}
dev.off()


#Create plots summarizing interactions##########################################
#This is Figure 7 in the reference.
#Markers within 20cM from major QTLs are regarded as the signals of the QTLs
#Interactions between markers that are 20cM apart from each other are regarded as the main effects
a <- matrix(rownames(Map), nrow = P, ncol = P, byrow = TRUE)
b <- matrix(rownames(Map), nrow = P, ncol = P, byrow = FALSE)
intlabel <- paste(a[lower.tri(a, diag = FALSE)], b[lower.tri(b, diag = FALSE)], sep="_")

#Give colors and pch to interactions between major QTLs
Neighbor <- as.list(numeric(length(MajorQTL)))
names(Neighbor) <- names(MajorQTL)
for(i in 1:length(MajorQTL)){
  Neighbor[[i]] <- paste0("M", which(Map$Chr == Map$Chr[MajorQTL[i]] & 
                                       Map$Dist < Map$Dist[MajorQTL[i]] + 20 & 
                                       Map$Dist > Map$Dist[MajorQTL[i]] - 20))
}

ColorInt <- matrix("gray80", ncol = P, nrow = P)
colnames(ColorInt) <- rownames(ColorInt) <- rownames(Map)
#Colors are:
#HD1 - HD6  : 2
#HD1 - DTH8 : 3
#HD1 - HD2  : 4
#HD6 - DTH8 : 5
#HD6 - HD2  : 6
#DTH8 - HD2 : 7
PchInt <- matrix(1, ncol = P, nrow = P)
colnames(PchInt) <- rownames(PchInt) <- rownames(Map)
#Pch are:
#HD1 - HD6  : 2
#HD1 - DTH8 : 3
#HD1 - HD2  : 4
#HD6 - DTH8 : 5
#HD6 - HD2  : 6
#DTH8 - HD2 : 8
for(m in 1:(length(MajorQTL) - 1)){
  Q1 <- names(MajorQTL)[m]
  N1 <- Neighbor[[m]]
  for(k in (m + 1):length(MajorQTL)){
    Q2 <- names(MajorQTL)[k]
    N2 <- Neighbor[[k]]
    
    if(Q1 == "HD1"){
      if(Q2 == "HD6") {colint <- 2; pchint <- 2}
      if(Q2 == "DTH8") {colint <- 3; pchint <- 3}
      if(Q2 == "HD2") {colint <- 4; pchint <- 4}
    }else{
      if(Q1 == "HD6"){
        if(Q2 == "DTH8") {colint <- 5; pchint <- 5}
        if(Q2 == "HD2") {colint <- 6; pchint <- 6}
      }else{
        if(Q1 == "DTH8"){
          if(Q2 == "HD2") {colint <- 7 ; pchint <- 8}     
        }else{
          cat("error", m, k, "\n")
        }
      }
    }
    
    ColorInt[N1, N2] <- colint
    ColorInt[N2, N1] <- colint
    PchInt[N1, N2] <- pchint
    PchInt[N2, N1] <- pchint
  }
}
ColorInt <- ColorInt[lower.tri(ColorInt, diag = FALSE)]
PchInt <- PchInt[lower.tri(PchInt, diag = FALSE)]

#function to jitter x axis accoding to data density
jitterxaxis <- function(v, p, width){
  #count neighbors
  n <- length(v)
  count <- numeric(n)
  for(i in 1:n){
    count[i] <- sum(v >= v[i] - width & v <= v[i] + width)
  }
  
  #jitter x axis according to count
  x_jitter <- numeric(n)
  for(i in 1:n){
    x_jitter[i] <- jitter(0, amount = log(count[i]) * p)
  }
  x_jitter
}

SHAP.int.RemoveMain <- as.list(numeric(E))
for(i in 1:E){
  
  value <- SHAP.int[[i]]
  
  #remove contamination of main effects
  for(m in 1:161){
    for(k in (m+1):162){
      if(Map$Chr[m] == Map$Chr[k]){
        if(abs(Map$Dist[m] - Map$Dist[k]) < 20){
          value[k, m] <- 0
        }
      }
    }
  }
  
  value <- value[lower.tri(value, diag=FALSE)]
  x_jitter <- jitterxaxis(value, p = 4e-2, width = 0.01)
  SHAP.int.RemoveMain[[i]] <- cbind(value, x_jitter)
  colnames(SHAP.int.RemoveMain[[i]]) <- c("SHAP", "Xaxis")
}

#plot as a tiff image
tiff("SummaryOfInteraction.tiff", unit="cm", res=400, width=17.5, height=8)
par(mar=c(3, 3, 1, 0.5))
par(mgp=c(2, 1, 0))
par(cex=0.5)
i <- 5
plot(SHAP.int.RemoveMain[[i]][, "Xaxis"] + 1, 
     SHAP.int.RemoveMain[[i]][, "SHAP"], 
     col = ColorInt, pch = PchInt,
     xlim = c(0.5, 9.5), ylim = c(0, 0.2),
     xlab = "", ylab = "Global SHAP interaction score", main = "", xaxt = "n")
for(i in 2:E){
  w <- c(5, 1, 7, 8, 9, 2, 4, 6, 3)[i]
  points(SHAP.int.RemoveMain[[w]][, "Xaxis"] + i, 
         SHAP.int.RemoveMain[[w]][, "SHAP"], 
         col = ColorInt, pch = PchInt)
}
axis(1, at = 1:E, labels = colnames(Y_Ori)[c(5, 1, 7, 8, 9, 2, 4, 6, 3)], las = 1)
legend(0.5, 0.2, legend = c("Hd1-Hd6","Hd1-DTH8","Hd1-Hd2","Hd6-DTH8","Hd6-Hd2","DTH8-Hd2"),
       col = 2:7, pch = c(2:6, 8), horiz = TRUE)
dev.off()


#Dependency plots###############################################################
#Dependency plots are drawn using my.shap.plot.dependence.R
#This function was developed from shap.plot.dependence in SHAPforxgboost
#such that regression lines and 95% CI can be drawn.

my.shap.plot.dependence <- function (data_long, x, y = NULL, color_feature = NULL, data_int = NULL, 
                                     dilute = FALSE, smooth = TRUE, size0 = NULL, add_hist = FALSE, 
                                     add_stat_cor = FALSE, alpha = NULL, jitter_height = 0, jitter_width = 0,
                                     Title = NULL, AddLine = TRUE,
                                     ...) 
{
  library(ggplot2)
  library(SHAPforxgboost)
  
  if (is.null(y)) 
    y <- x
  data0 <- data_long[variable == y, .(variable, value)]
  data0$x_feature <- data_long[variable == x, rfvalue]
  if (!is.null(color_feature) && color_feature == "auto") {
    color_feature <- strongest_interaction(X0 = data0, Xlong = data_long)
  }
  if (!is.null(color_feature)) {
    data0$color_value <- data_long[variable == color_feature, 
                                   rfvalue]
  }
  if (!is.null(data_int)) 
    data0$int_value <- data_int[, x, y]
  nrow_X <- nrow(data0)
  if (is.null(dilute)) 
    dilute = FALSE
  if (dilute != 0) {
    dilute <- ceiling(min(nrow(data0)/10, abs(as.numeric(dilute))))
    set.seed(1234)
    data0 <- data0[sample(nrow(data0), min(nrow(data0)/dilute, 
                                           nrow(data0)/2))]
  }
  if (x == "dayint") {
    data0[, `:=`(x_feature, as.Date(data0[, x_feature], format = "%Y-%m-%d", 
                                    origin = "1970-01-01"))]
  }
  if (is.null(size0)) {
    size0 <- if (nrow(data0) < 1000L) 
      1
    else 0.4
  }
  if (is.null(alpha)) {
    alpha <- if (nrow(data0) < 1000L) 
      1
    else 0.6
  }
  plot1 <- ggplot(data = data0, 
                  aes(x = x_feature, y = if (is.null(data_int)) value else int_value, 
                      color = if (!is.null(color_feature)) color_value else NULL)) + 
    geom_jitter(size = size0, 
                width = jitter_width,
                height = jitter_height,
                alpha = alpha, ...) + 
    labs(y = if (is.null(data_int)) {
      paste0("SHAP value for ", y)
    } else {
      paste0("SHAP interaction values for\n", x, " and ", y)
    }, 
    x = x, 
    color = if (!is.null(color_feature)) paste0(color_feature, "\n", "(Feature value)") else NULL) + 
    scale_color_gradient(low = "#FFCC33", 
                         high = "#6600CC",
                         guide = guide_colorbar(barwidth = 10, barheight = 0.3)) + 
    theme_bw() + 
    theme(legend.position = "bottom",
          plot.margin = unit(c(0.5, 2, 0.5, 0.5), "lines"),
          legend.title = element_text(size = 10), 
          legend.text = element_text(size = 8))
  if (smooth) {
    plot1 <- plot1 + geom_smooth(method = "loess", color = "red", 
                                 size = 0.4, se = FALSE)
  }
  #plot1 <- plot.label(plot1, show_feature = x)
  if (add_stat_cor) {
    plot1 <- plot1 + ggpubr::stat_cor(method = "pearson")
  }
  if (add_hist) {
    plot1 <- ggExtra::ggMarginal(plot1, type = "histogram", 
                                 bins = 50, size = 10, color = "white")
  }
  if (AddLine){
    data1 = data.frame(x_feature = data0[data0$color_value == -1, x_feature],
                       int_value = data0[data0$color_value == -1, int_value])
    plot1 <- plot1 + geom_smooth(data = data1, method = "lm", color = "#FFCC33", size=0.5)
    data1 = data.frame(x_feature = data0[data0$color_value == 1, x_feature],
                       int_value = data0[data0$color_value == 1, int_value])
    plot1 <- plot1 + geom_smooth(data = data1, method = "lm", color = "#6600CC", size=0.5)
  }
  plot1 <- plot1 + ggtitle(Title)
  plot1
}

#Create required objects
#Because the function requires the outputs of xgboost,
#run xgboost once and replace the values with the averages obtained with bootstrap.
Train <- xgb.DMatrix(data = as.matrix(X_Ori), label = Y_Ori[,1])
Val <- xgb.DMatrix(data = as.matrix(X_Ori), label = Y_Ori[,1])
model <- xgb.train(data = Train,
                   nrounds = 1000,
                   watchlist = list(train = Train, eval = Val),
                   early_stopping_rounds = 3,
                   verbose = 0)

#Extract the shap values
shap_long <- shap.prep(xgb_model = model, X_train = as.matrix(X_Ori))
shap_int <- shap.prep.interaction(xgb_mod = model, X_train = as.matrix(X_Ori))
#These objects are used as templates

#read bootstrap averages of local importance and replace values with them.
Shap_long <- as.list(numeric(E))
names(Shap_long) <- colnames(Y_Ori)
Shap_int <- as.list(numeric(E))
names(Shap_int) <- colnames(Y_Ori)
M <- unique(as.character(shap_long$variable))
for(env in colnames(Y_Ori)){
  
  fill <- !is.na(Y_Ori[, env])
  
  #main effects
  v <- matrix(unlist(read.csv(paste0("SHAP-XGB.at.", env, "/LocalImportance/Shap_values_average.csv"),
                              row.names=1)), ncol = P)
  colnames(v) <- paste0("M", 1:P)
  for(m in M){
    shap_long$value[shap_long$variable == m][fill] <- v[, m]
    shap_long$value[shap_long$variable == m][!fill] <- NA
  }
  Shap_long[[env]] <- shap_long
  
  #interaction effect
  for(i in 1:162){
    v <- matrix(unlist(read.csv(paste0("SHAP-XGB.at.", env, "/LocalImportance/Shap_int_average_", i, ".csv"),
                                row.names=1)), ncol=162)
    colnames(v) <- paste0("M", 1:162)
    for(m in 1:162){
      shap_int[fill, i, m] <- v[, m] * 2
      shap_int[!fill, i, m] <- NA
    }
    Shap_int[[env]] <- shap_int
  }
}

#See top interactions at Tsukuba2009 and Ishigaki2008 as examples
##Tsukuba2009
target <- E
v <- order(abs(SHAP.int[[target]]), decreasing = TRUE)[1:10]
pos <- cbind(v %% P, v %/% P + 1)
pos <- cbind(pos, diag(SHAP.int[[target]][pos[, 1], pos[, 2]]))
for(j in 1:nrow(pos)) cat("#", pos[j, 1], pos[j, 2], pos[j, 3], "\n")
##Following interactions were included
# 112 92 0.1329156  =>DTH8-HD1
#  92 60 0.1015208  =>HD1-HD6
# 108 92 0.06825149 =>HD2-HD1

##Ishigaki2008
target <- 4
v <- order(abs(SHAP.int[[target]]), decreasing = TRUE)[1:10]
pos <- cbind(v %% P, v %/% P + 1)
pos <- cbind(pos, diag(SHAP.int[[target]][pos[, 1], pos[, 2]]))
for(j in 1:nrow(pos)) cat("#", pos[j, 1], pos[j, 2], pos[j, 3], "\n")
##Following interactions were included
# 111 108 0.06290753 =>DTH8-HD2

#Create dependency plots for these interactions
Wj <- 0.08
##Tsukuba2009 
#### 112 92 (DTH8-HD1)
data_long <- Shap_long[[E]]
levels(data_long$variable)[levels(data_long$variable) == "M112"] <- "DTH8"
levels(data_long$variable)[levels(data_long$variable) == "M92"] <- "Hd1"
data_int <- Shap_int[[E]]
dimnames(data_int)[[2]][c(112, 92)] <- c("DTH8", "Hd1")
dimnames(data_int)[[3]][c(112, 92)] <- c("DTH8", "Hd1")
plot1 <- my.shap.plot.dependence(data_long = data_long,
                                 data_int = data_int,
                                 x = "DTH8",
                                 y = "Hd1",
                                 color_feature = "Hd1",
                                 smooth = FALSE,
                                 size0 = 0.5,
                                 alpha = 1,
                                 jitter_width = Wj,
                                 Title = "Tsukuba2009")

##Tsukuba2009 
#### 92 60 (HD1-HD6)
data_long <- Shap_long[[E]]
levels(data_long$variable)[levels(data_long$variable) == "M92"] <- "Hd1"
levels(data_long$variable)[levels(data_long$variable) == "M60"] <- "Hd6"
data_int <- Shap_int[[E]]
dimnames(data_int)[[2]][c(92, 60)] <- c("Hd1", "Hd6")
dimnames(data_int)[[3]][c(92, 60)] <- c("Hd1", "Hd6")
plot2 <- my.shap.plot.dependence(data_long = data_long,
                                 data_int = data_int,
                                 x = "Hd6",
                                 y = "Hd1",
                                 color_feature = "Hd1",
                                 smooth = FALSE,
                                 size0 = 0.5,
                                 alpha = 1,
                                 jitter_width = Wj,
                                 Title = "Tsukuba2009")

##Tsukuba2009 
#### 108 92 (HD2-HD1)
data_long <- Shap_long[[E]]
levels(data_long$variable)[levels(data_long$variable) == "M108"] <- "Hd2"
levels(data_long$variable)[levels(data_long$variable) == "M92"] <- "Hd1"
data_int <- Shap_int[[E]]
dimnames(data_int)[[2]][c(108, 92)] <- c("Hd2", "Hd1")
dimnames(data_int)[[3]][c(108, 92)] <- c("Hd2", "Hd1")
plot3 <- my.shap.plot.dependence(data_long = data_long,
                                 data_int = data_int,
                                 x = "Hd2",
                                 y = "Hd1",
                                 color_feature = "Hd1",
                                 smooth = FALSE,
                                 size0 = 0.5,
                                 alpha = 1,
                                 jitter_width = Wj,
                                 Title = "Tsukuba2009")

##Ishigaki2008 
#### 111 108 (DTH8-HD2)
data_long <- Shap_long[[4]]
levels(data_long$variable)[levels(data_long$variable) == "M111"] <- "DTH8"
levels(data_long$variable)[levels(data_long$variable) == "M108"] <- "Hd2"
data_int <- Shap_int[[4]]
dimnames(data_int)[[2]][c(111, 108)] <- c("DTH8", "Hd2")
dimnames(data_int)[[3]][c(111, 108)] <- c("DTH8", "Hd2")
plot4 <- my.shap.plot.dependence(data_long = data_long,
                                 data_int = data_int,
                                 x = "DTH8",
                                 y = "Hd2",
                                 color_feature = "Hd2",
                                 smooth = FALSE,
                                 size0 = 0.5,
                                 alpha = 1,
                                 jitter_width = Wj,
                                 Title = "Ishigaki2008")

#combine ggplot objects using library patchwork
#This is Figure 8 in the reference.
library(patchwork)
tiff("DependencyPlot.tiff", unit="cm", res=400, width=18.5, height=18)
(plot1 | plot2)/(plot3 | plot4)
dev.off()

