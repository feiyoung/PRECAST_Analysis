######### For this Scenario  with High batch effects, cluster label and spatial coordinates are generated from a Potts model.
######### 65*65, 60*60, 60*60 spots corresponding to rectangular lattices  from  a $K$-state ($K$=7)  Potts model
######### Thus, the platform is ST.

setwd("./Simulation/simulated_datasets/")
load('simulated_datasets_Scenario1_high.rds')

library(PRECAST)
XList <- lapply(XtList, function(x) log(1+x))

XList1 <- lapply(XList, scale, scale=FALSE) 
platform_use <- "ST"; maxIter_use <- 30

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
