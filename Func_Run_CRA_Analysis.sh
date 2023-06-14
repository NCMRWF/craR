#!/bin/bash
module load gnu/R/3.4.3
##################################################
#for i in {1..31..1}
#do
#echo $i

date=20210517
RUNPATH=/home/greeshma/TestArea/CRA_Test/CRA_TEST/RMarkdown/CRA_Analysis_Package_Mod/CRA_Package_V1
ActINPUT_OBS=obs_${date}.nc
ActINPUT_MODEL=model_${date}.nc
OUT_PATH=/home/greeshma/TestArea/CRA_Test/CRA_TEST/RMarkdown/CRA_Analysis_Package_Mod/CRA_Package_V1
Threshold=40
case=CTL
OUT_PREFIX=${date}_Thr${Threshold}_${case}
Ngrids=6 #### No of grids to be shifted
##################################################

mkdir -p $OUT_PATH
cd $RUNPATH
cat > Calc_CRA_${date}_${case}_Th${Threshold}.R << EOF
library(raster)
library("SpatialVx")
###
rotate <- function(x) t(apply(x, 2, rev))	
antirotate <- function(x) apply(t(x),2,rev)
################## INPUTS ##############################
ActOrf <- raster("${ActINPUT_OBS}")
ActMrf <- raster("${ActINPUT_MODEL}")
################## Identify the Matched Objects ########
p<-raster::as.matrix(ActOrf)
q<-raster::as.matrix(ActMrf)

p_rot<-rotate(p)
q_rot<-rotate(q)

hold <- make.SpatialVx(p_rot,q_rot)
print(hold)
look1 <- FeatureFinder(hold, thresh = ${Threshold},min.size = c(30,30))  ### have to test different min.size and Threshold

look2 <- minboundmatch( look1, type = "multiple", mindist =5) ### have to test mindist
look <- MergeForce( look2 )
#plot(look)

SummFeats=summary(look)
ObsSummary<-SummFeats\$X
ForecastSummary<-SummFeats\$Y
num_featO=dim(ObsSummary)
num_featO<-num_featO[1]
print(paste0('No. of identified features in Obs= ',num_featO[1]))
num_featM=dim(ForecastSummary)
num_featM<-num_featM[1]
print(paste0('No. of identified features in Forecast= ',num_featM[1]))

### get the number of matched features ###
MatchObj<-print(look\$matches)
Num_Match<-MatchObj[,1]
TotalNumFeat=length(Num_Match)
MatchFeatModel<-MatchObj[,1]
MatchFeatObs<-MatchObj[,2]

############### Detaching Packages which stop running ´shift´ function ##########
detach(package:SpatialVx)
detach(package:turboEM)
detach(package:quantreg)
detach(package:SparseM)
detach(package:numDeriv)
detach(package:doParallel)
detach(package:iterators)
detach(package:foreach)
detach(package:smatr)
detach(package:smoothie)
detach(package:fields)
detach(package:maps)
detach(package:spam)
detach(package:dotCall64)
detach(package:spatstat)
detach(package:rpart)
detach(package:nlme)
detach(package:spatstat.data)
library(fields)
##
for (ii in 1:TotalNumFeat){
        nam=paste("Area_ObF",ii,sep="")
        assign(nam,SummFeats\$X[ii,3])
        nam2=paste("Area_ModelF",ii,sep="")
        assign(nam2,SummFeats\$Y[ii,3])
        ind_Obs<-(look\$X.labeled)
        ind_Model<-(look\$Y.labeled)
        ind_Obs[ind_Obs!=MatchFeatObs[ii]]<-NA
        ind_Model[ind_Model!=MatchFeatModel[ii]]<-NA
        image.plot(ind_Obs)
        image.plot(ind_Model)

	tempObs=as.matrix(ActOrf)
	tempObs<-rotate(tempObs)
	new_Obs<-(ind_Obs*tempObs)/MatchFeatObs[ii]
#	image.plot(new_Obs)
	newrastO<-ActOrf
	new_Obs<-antirotate(new_Obs)
	values(newrastO)=new_Obs
	plot(newrastO)

	tempModel=as.matrix(ActMrf)
	tempModel<-rotate(tempModel)
	new_Model<-(ind_Model*tempModel)/MatchFeatModel[ii]
	newrastM<-ActMrf
	new_Model<-antirotate(new_Model)
	values(newrastM)=new_Model
	plot(newrastM)

orf<-newrastO
mrf<-newrastM

######################
Obj_Obs=paste('Obs_${date}','_Object',ii,'.nc',sep="")
Obj_Model=paste('Model_${date}','_Object',ii,'.nc',sep="")
ObObjnc <- writeRaster(orf, filename=Obj_Obs, format="CDF", overwrite=TRUE)
MoObjnc <- writeRaster(mrf, filename=Obj_Model, format="CDF", overwrite=TRUE)
print("-----------------------------------------------------------------")
print(paste("Object files used for analysis are",Obj_Obs,"and",Obj_Model))
print("-----------------------------------------------------------------")
##################################################################
source("CRA_Err_Decomp_Features.R")
CRA_Err_Decomp("${ActINPUT_OBS}","${ActINPUT_MODEL}",Obj_Obs,Obj_Model,thr=${Threshold},objnum=1,"${OUT_PREFIX}",Ngrids=${Ngrids})
}
EOF

Rscript Calc_CRA_${date}_${case}_Th${Threshold}.R
cp Calc_CRA_${date}_${case}_Th${Threshold}.R ${OUT_PATH}
mv Rplots.pdf ${OUT_PATH}/CRA_Features_${date}_${case}.pdf
mv ${date}*.txt ${OUT_PATH}
#done
exit
