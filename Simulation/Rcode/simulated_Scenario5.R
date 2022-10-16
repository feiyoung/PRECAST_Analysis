###### setwd("F:\\Research paper\\IntegrateDRcluster\\AnalysisCode\\PRECAST_Analysis")

######### For this Scenario  with high batch effects, cluster label and spatial coordinates are based on  
######### three DLFPC datasets (ID: 151507, 151669 and 151673) from each of three donors (Visium platform)

setwd("./Simulation/simulated_datasets/")
for(r in 1:3){
  load(paste0("simulated_datasets_Scenario5_sample", r, ".rds"))
}
### In the rds data file, posList is a list including spaital coordinates for all samples,
### yList is a list including true cluster labels for all samples; 
### X1, X2, and X3 are the Logcount matrices for three samples.
str(posList)
str(yList)

library(PRECAST)
XList <- list(X1, X2, X3)
XList1 <- lapply(XList, scale, scale=FALSE) 
platform_use <- "Visium"; maxIter_use <- 30

### We run the simulation experiments using the low-level main function ICM.EM();
### User can use the high-level function PRECAST() by following the tutorial at https://feiyoung.github.io/PRECAST/.

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

