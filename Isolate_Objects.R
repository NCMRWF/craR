Isolate_Objects<-function(ActOrf,ActMrf,objOBS,objMODEL,thr,minsize,sepdist){

library(raster)
library("SpatialVx")
###
rotate <- function(x) t(apply(x, 2, rev))	
antirotate <- function(x) apply(t(x),2,rev)
################## INPUTS ##############################
ActOrf <- raster(ActOrf)
ActMrf <- raster(ActMrf)
################## Identify the Matched Objects ########
p<-raster::as.matrix(ActOrf)
q<-raster::as.matrix(ActMrf)

p_rot<-rotate(p)
q_rot<-rotate(q)

hold <- make.SpatialVx(p_rot,q_rot)
print(hold)
look1 <- FeatureFinder(hold, thresh = thr,min.size = c(minsize,minsize))  

look2 <- minboundmatch( look1, type = "multiple", mindist =sepdist) 
look <- MergeForce( look2 )
#plot(look)

SummFeats=summary(look)
ObsSummary<-SummFeats$X
ForecastSummary<-SummFeats$Y
num_featO=dim(ObsSummary)
num_featO<-num_featO[1]
print(paste0('No. of identified features in Obs= ',num_featO[1]))
num_featM=dim(ForecastSummary)
num_featM<-num_featM[1]
print(paste0('No. of identified features in Forecast= ',num_featM[1]))

### get the number of matched features ###
MatchObj<-print(look$matches)
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
        assign(nam,SummFeats$X[ii,3])
        nam2=paste("Area_ModelF",ii,sep="")
        assign(nam2,SummFeats$Y[ii,3])
        ind_Obs<-(look$X.labeled)
        ind_Model<-(look$Y.labeled)
        ind_Obs[ind_Obs!=MatchFeatObs[ii]]<-NA
        ind_Model[ind_Model!=MatchFeatModel[ii]]<-NA
        image.plot(ind_Obs)
        image.plot(ind_Model)

	tempObs=as.matrix(ActOrf)
	tempObs<-rotate(tempObs)
	new_Obs<-(ind_Obs*tempObs)/MatchFeatObs[ii]
	image.plot(new_Obs)
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
Obj_Obs=paste(objOBS,'_Object',ii,'.nc',sep="")
Obj_Model=paste(objMODEL,'_Object',ii,'.nc',sep="")
ObObjnc <- writeRaster(orf, filename=Obj_Obs, format="CDF", overwrite=TRUE)
MoObjnc <- writeRaster(mrf, filename=Obj_Model, format="CDF", overwrite=TRUE)
print(paste("Object files used for analysis are",Obj_Obs,"and",Obj_Model))
}
}
