
# Here is my SDM script template1 for various projects

# The first thing to do is to doanload and load required packages

library(sp) #required for reading and mapping coordinates
library(raster) # required for reading and working with raster data
library(rgdal) # required for reading and working with shapefiles
library(rgeos) #calculating gdifference
library(tidyverse) # requited for manipulating data

# Here are the different modeling algorithm packages to load
library(dismo) #working with SDMs
library(mgcv) #To run GAM
library (rJava)# equired for running Maxent
library(ROCR) # required for model cross validations i.e., AUC and TSS
install.packages("rJava")
install.packages("ROCR")


#.......................SDM Data Analysis..........................

# the first thing to do is to load our species occurrence data as a csv file. Here are analysing chimpanzee data from a PA in central africa
occur <- read.csv("chimps.csv", header = T)
head(occur)

# Next is to convert our csv data to spatial points objects and with the right ptojection
coords <-SpatialPoints(occur[,2:3], proj4string = CRS("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"))

# Alternatively, we could read our species occurrence data directly as a shapefile, if the datasets are available in shapefile format

#the next thing is to poad and stack predictor datasets as a list from the data file path
#make sure all rasters have thesame spatial extents (rows and columns values), resolutions, and coordinate reference systems
files <- list.files(path = "F:/SDM_template/", pattern = ".tif$", full.names=T)
files
predictors <- raster::stack(files)
plot(predictors)

# Our next step is to load the study area polygon for mapping our species occurrence and stacked raster data
# ensure the study area is projected to same spatial extent as your raster data
studyarea <- readOGR("study_area.shp")
plot(studyarea)


# My next step now is to create pseudo absence or background points
# Before I do that, I will first of all create buffers around the presence locations in order to ensure the presence points do not overlap with absence points created
buffer <- gBuffer(coords, width=1000, byid=TRUE)
buffer1 <- polygons(buffer) # this ensures the buffered points are inserted into a polygon object 

# the next step now is to eliminate all buffered areas from the study area 
buff_elim <- gDifference(studyarea, buffer1) 

# I can now Sample my random or background points across the elimated buffered areas so as to avoid overlaps with presence locations
set.seed(0) # set seed to create the same set every time
random_points <- spsample(buff_elim, 1000, type='random')

# It is a good idea to always Visualize the presence and pseudoabsences points in your study area before proceeding with further analysis
plot(predictors, 1, add=TRUE) # Here I am testing the points plot on the first raster layer from stacked predictors
points(random_points, col='red')
points(coords, col='green')
plot(studyarea, border= "blue", add = TRUE)


#.To perform my SDM analysis, I will extract presence and absence point locations from my predictors, and bind these datasets as a dataframe

# I will start by extracting raster values for presence points 
pres_vals <- raster::extract(predictors, coords)
pres_vals

# Now for absence or random points
abs_vals <- raster::extract(predictors, random_points)
abs_vals


# Next is to append and bind the extracted presence data to the presence and absence points as data frames 
xy_coords <- data.frame(coords) #Converts spatial points to data
pres_vals1 <- cbind(xy_coords, pres_vals)

# employ a similar approach for absence data
random_coords <- as.data.frame(random_points)
abs_vals1 <- cbind(random_coords, abs_vals)
abs_vals1

# The next step is to Create a new column with the value 1 for presence records and 0 for absence records
occurrence <- c(rep(1, nrow(pres_vals1)), rep(0, nrow(abs_vals1)))

# Next is to append and bind the 0 and 1 values to the extracted values for all points as a dataframe
sdm_data <- data.frame(cbind(occurrence, rbind(pres_vals1, abs_vals1)))
sdm_data

# it could be that the sdm_data contains not applicable values. Is is therefore important to check and omit such datasets in order to avoid model prediction bias 
NA_available <- is.na(sdm_data)
summary(NA_available) # results generated will have false and true, with false meaning there are no NA and True meaning there are NA. In my case there some few NA

# my next step is to delete all NA
sdm_data1 <- na.omit(sdm_data)
sdm_data1

# Next is to verify if all NA have been deleted
NA_available1 <- is.na(sdm_data1)
summary(NA_available1) # results 
# all NA have been deleted in my case, implying we can proceed with the modeling analysis

# Before training models, it is a good idea to always check for predictor correlations
cor<- sdm_data1 %>% select(ecoguard_patrols1, elevation1, human_pressure1) # these are the predictors in my dataset
cor(cor) #verify if values <0.7?


#...............SDM Modeling................................----

# I am going to test three deifferent models, including glm, gam, and MaxEnt

# I will start with gam model 

# the first step in the modeling process is to Separate the presence values from the absences values
pres <- sdm_data1[sdm_data1$occurrence == 1,]#presence
abs <- sdm_data1[sdm_data1$occurrence == 0,]#absence

# the next step is to create testing and training datasets
# I am creating a 5 k-fold groups (i.e., splitting the data by 20% into 5 groups)
k <- 5
grouppres <- kfold(pres, k) # this argument generates a group for presence data
groupabs <- kfold(abs, k) # this argument generates a group for absence or background data
total<-c(grouppres, groupabs)
total
sdm_data1$kfold<-total
sdm_data1$kfold

#Next is to train my data for SDM
train <- sdm_data1[!sdm_data1$kfold == 1, ]  # we can put any number between 1 to 5 considering we have a five seperate data folds
test <- sdm_data1[sdm_data1$kfold == 1, ]    # similar input here

# Next is to run the gam model
# lests try for testing group 1
m1 <- gam(occurrence ~  s(elevation1) + s(human_pressure1) + s(ecoguard_patrols1), 
          data = train,  family= binomial("logit"))
mgcv::summary.gam(m1) #gets the model summary

# Next is to get the predicted model values for the test data 
predicted_values_gam <- as.data.frame(predict(m1, test, type="response"))
predicted_values_gam
Values <- as.data.frame(test$occurrence) #get the ones and zeros into a data frame
values


# Next is to validate the gam model by getting the ROC-AUC values
pred_roc_gam <- ROCR::prediction(predicted_values_gam, test$occurrence) 
pred_roc_gam
perf_roc_gam <- performance(pred_roc_gam,"tpr","fpr")  # this argument extracts the true and false postive values of the model
perf_roc_gam
plot(perf_roc_gam,col="grey",lwd=2,type="l") # this argument is to Visualize the ROC in the form of a curve
perf_auc_gam <- performance(pred_roc_gam,"auc") # this argument generates the  AUC
perf_auc_gam
AUC_gam <- as.numeric(slot(perf_auc_gam,"y.values")) # creates Values needed for the ensemble
AUC_gam

# the next thing is to get the TSS value
# TSS is point at which specificity + sensitivity is maximum
perf <- performance(pred_roc_gam, "sens", "spec")
perf
TSS_GAM<-(max(perf@x.values[[1]] + perf@y.values[[1]]))-1
TSS_GAM

# Next is to create suitability or SDM maps and save results
gam_suitability <- predict(predictors, m1, progress='text', type="response")
plot(gam_suitability)
writeRaster(gam_suitability, filename="chimps_gam_prediction.tif", format="GTiff", overwrite=TRUE) # saves my predicted map as a tif file



#...............MaxEnt Analysis..........................----

#The MaxEnt model worls with the package "dismo", which has Maxent enbedded as a function
# the process also requires java, hence, the need to first install the rjava package installed 
# Now, install MaxEnt and copy the maxent.jar file and paste in the R dismo folder,
# example: place it in D:/TempR/R/win-library/3.5/dismo/java/maxent.jar

.libPaths()#verifies path for installed packages

# use a random value to test the maxent model (here, I am using the value 1)
test1=1  
train_p <- pres[grouppres != test1, c("x","y")]  #Training presence points
train_a <- abs[groupabs != test1, c("x","y")]  #Training absence points

# running the MaxEnt model
Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_271') 
library(rJava)

max1 <- maxent(predictors, p=train_p, a=train_a, args=c('responsecurves=true', 'writeplotdata=true'))
max1 #check rsults results
plot(max1) #plots the percentage variable contributions

# Next is to get predicted model values for the test data 
predicted_values_max <- as.data.frame(predict(max1, test, type="response")) #gets model predicted values
Values<-as.data.frame(test$occurrence) #get binary response
Values1<-as.numeric(test$occurrence)#use numeric argument if dataframe does not work
predicted_values_max1 <- as.numeric(predict(max1, test, type="response"))#use numeric argument if dataframe does not work

#Next is o conduct the MaxEnt model validation by getting the ROC-AUC 
pred.rocr.max <- ROCR::prediction(predicted_values_max, Values) 
perf.max <- performance(pred.rocr.max,"tpr","fpr")  #Extracting the true and false positives values
plot(perf.max,col="grey",lwd=2,type="l") # Visualize the extracted values in a curve
perf.auc.max <- performance(pred.rocr.max,"auc") #getting the AUC valuess
perf.auc.max
AUC_MAX <- as.numeric(slot(perf.auc.max,"y.values")) # getting the Value needed for model ensemble
AUC_MAX

# Next is to get the TSS values
perf2 <- performance(pred.rocr.max, "sens", "spec")
TSS_MAX<-(max(perf2@x.values[[1]] + perf2@y.values[[1]]))-1
TSS_MAX

# Final step is to create a maxent suitability or SDM map and export results as a tif file
pred_max = predict(max1, predictors, progress='text', timer=TRUE)
plot(pred_max)
writeRaster(pred_max, filename="chimps_suitability_MAXENT.tif", format="GTiff", overwrite=TRUE) #saves the suitability map!



#...............GLM Analysis..........................----

m2 <- glm(occurrence ~ elevation1 + human_pressure1 + ecoguard_patrols1, 
          data = train,  family= binomial("logit"))
summary(m2)

# Next is to get the predicted model values for the test data 
predicted_values_glm <- as.data.frame(predict(m2, test, type="response"))
predicted_values_glm
Values <- as.data.frame(test$occurrence) #get the ones and zeros into a data frame
values

# Next is to validate the glm model by getting the ROC-AUC values
pred_roc_glm <- ROCR::prediction(predicted_values_glm, test$occurrence) 
pred_roc_glm
perf_roc_glm <- performance(pred_roc_glm,"tpr","fpr")  # this argument extracts the true and false postive values of the model
perf_roc_glm
plot(perf_roc_glm,col="grey",lwd=2,type="l") # this argument is to Visualize the ROC in the form of a curve
perf_auc_glm <- performance(pred_roc_glm,"auc") # this argument generates the  AUC
perf_auc_glm
AUC_glm <- as.numeric(slot(perf_auc_glm,"y.values")) # creates Values needed for the ensemble
AUC_glm

# the next thing is to get the TSS value
perf_glm <- performance(pred_roc_glm, "sens", "spec")
perf
TSS_GLM<-(max(perf_glm@x.values[[1]] + perf_glm@y.values[[1]]))-1
TSS_GLM

# Next is to create GLM suitability or SDM maps and save results
glm_suitability <- predict(predictors, m2, progress='text', type="response")
plot(glm_suitability)
writeRaster(glm_suitability, filename="chimps_glm_prediction.tif", format="GTiff", overwrite=TRUE) # saves my predicted map as a tif file



#..........................Ensemble modelling.................----
# Because all models performed relatively well, I will now create a weighted model ensemble for all three SDMs
AUC_sum <-  AUC_gam + AUC_glm + AUC_MAX   # argument sums the individual AUC values predicted in each model
weight_AUC_gam <- (AUC_gam/AUC_sum)           # this argument Calculates the weighted AUC's for each model 
weight_AUC_glm <- (AUC_glm/AUC_sum)
weight_AUC_max <- (AUC_MAX/AUC_sum)

# Next is to create the model ensemble
model_esmble <- ((gam_suitability * weight_AUC_gam) + (pred_max * weight_AUC_max) 
                 + (glm_suitability * weight_AUC_glm)) 
plot(model_esmble)
writeRaster(model_esmble, filename="chimps_ensemble.tif", format="GTiff", overwrite=TRUE) #saves my ensemble model predictions



#Alternative ensemble calculation

suitability_avg <- ((gam_suitability) + (pred_max) + (glm_suitability))/3 
plot(suitability_avg)

AUC_avr <-  (AUC_gam + AUC_glm + AUC_MAX)/3
AUC_avr


# Next is to get ROC-AUC values for the weighted model
# We start by extracting the raster values for the presences and absences data
data_pts <- sdm_data1[,c("x","y")]   # gets the coordinates
esm_pred_vals <- raster::extract(model_esmble, data_pts) #Extract ensemble values for each location
esm_pred_vals <- as.data.frame(esm_pred_vals)
pres_abs_vals <- sdm_data1[,"occurrence"]   #this defines the "actual" presence values
pres_abs_vals <- as.data.frame(pres_abs_vals)
predict<-cbind(esm_pred_vals, pres_abs_vals)

#...Get the ROC-AUC value ----
pred.rocr.esm <- ROCR::prediction(predict$esm_pred_vals,predict$pres_abs_vals)
perf.esm<- performance(pred.rocr.esm,"tpr","fpr")
plot(perf.esm,col="grey",lwd=2,type="l")
perf.auc.esm<- performance(pred.rocr.esm,"auc")
ESM_AUC <- as.numeric(slot(perf.auc.esm,"y.values"))
ESM_AUC

# Get the TSS value 
perf <- performance(pred.rocr.esm, "sens", "spec")
TSS_all<-(max(perf@x.values[[1]] + perf@y.values[[1]]))-1
TSS_all

