######### For this Scenario  with high batch effects, cluster label and spatial coordinates arebased on  
######### three DLFPC datasets (ID: 151507, 151669 and 151673) from each of three donors (Visium platform)

setwd("./Simulation/simulated_datasets/")
load('simulated_datasets_Scenario2_high.rds')

### In the rds data file, posList is a list including spaital coordinates for all samples,
### yList is a list including true cluster labels for all samples; 
### XtList is a list including the count matrix for all samples
str(posList)
str(yList)


library(PRECAST)
XList <- lapply(XtList, function(x) log(1+x))

XList1 <- lapply(XList, scale, scale=FALSE) 
platform_use <- "Visium"; maxIter_use <- 30

### We run the simulation experiments using the low-level main function ICM.EM();
### User can use the high level function PRECAST() by following the tutorial at https://feiyoung.github.io/PRECAST/.

## Use fixed K
hq <- 15; hK <- 7
resList <- ICM.EM(XList1,q=hq, K=hK, posList= posList,
                  platform = platform_use,maxIter = maxIter_use,
                  Sigma_equal =F, coreNum=length(hK))
reslist <- selectModel(resList) ## select a best model
str(resList)
str(reslist)
## use ?ICM.EM to see the explanation of results in resList object.

## Choose K using MBIC
hq <- 15; hK <- 2:9
resList <- ICM.EM(XList1,q=hq, K=hK, posList= posList,
                  platform = platform_use,maxIter = maxIter_use,
                  Sigma_equal =F, coreNum=length(hK))
reslist <- selectModel(resList) ## select a best model
str(reslist)

