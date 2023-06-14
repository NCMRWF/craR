#!/bin/bash
module load gnu/R/3.4.3
##################################################
#for i in {1..31..1}
#do
#echo $i

date=20210517
RUNPATH=/home/greeshma/TestArea/CRA_Test/CRA_TEST/RMarkdown/CRA_Analysis_Package_Mod/CRA_Package_V1
OUT_PATH=${RUNPATH}
Threshold=40
case=CTL
ObjNum=1 
Rainlimt=850  ## for colorbar
RainlimtScatter=850  ## for scatter plot
OUT_PREFIX=${date}_Thr${Threshold}_${case}_$ObjNum
shpfile=/home/greeshma/TestArea/CRA_Test/CRA_TEST/RMarkdown/CRA_Analysis_Package_Mod/CRA_Package_V1
figname=CRAplot_${date}_Th${Threshold}_${case}_${ObjNum}.png
##################################################
INPUT_OBS=${RUNPATH}/Obs_${date}_Object${ObjNum}.nc
INPUT_MODEL=${RUNPATH}/Model_${date}_Object${ObjNum}.nc
lat1=10
lat2=28
lon1=60
lon2=80

############ From Stat Files #################################
CRA_threshold="${Threshold}mm/day"
################# Obs Vs Actual ############# 
anpts=$(cat ${OUT_PREFIX}_StatsOrig.txt | awk -F " " '{print $1}')
bnpts=$(cat ${OUT_PREFIX}_StatsOrig.txt | awk -F " " '{print $2}')
aavg=$(cat ${OUT_PREFIX}_StatsOrig.txt | awk -F " " '{print $3}')
bavg=$(cat ${OUT_PREFIX}_StatsOrig.txt | awk -F " " '{print $4}')
amax=$(cat ${OUT_PREFIX}_StatsOrig.txt | awk -F " " '{print $5}')
bmax=$(cat ${OUT_PREFIX}_StatsOrig.txt | awk -F " " '{print $6}')
avol=$(cat ${OUT_PREFIX}_StatsOrig.txt | awk -F " " '{print $7}')
bvol=$(cat ${OUT_PREFIX}_StatsOrig.txt | awk -F " " '{print $8}')
################## Actual Vs Shifted ###########################
Orms=$(cat ${OUT_PREFIX}_AllStats.txt | awk -F " " '{print $11}')
Nrms=$(cat ${OUT_PREFIX}_AllStats.txt | awk -F " " '{print $12}')
Ocorr=$(cat ${OUT_PREFIX}_AllStats.txt | awk -F " " '{print $13}')
Ncorr=$(cat ${OUT_PREFIX}_AllStats.txt | awk -F " " '{print $14}')
################### CRA Error Components ###########################
disErr=$(cat ${OUT_PREFIX}_cra_decomp.txt | awk -F " " '{print $2}')
volErr=$(cat ${OUT_PREFIX}_cra_decomp.txt | awk -F " " '{print $3}')
patternErr=$(cat ${OUT_PREFIX}_cra_decomp.txt | awk -F " " '{print $4}')
######################

cd $RUNPATH
###
cat > ${date}_PlotCRA_${case}_${ObjNum}.R << EOF
# Read raster data and isolate objects by thresholds
# Obtain the area of objects.
library(raster)
library("RColorBrewer")
#################################################
orf <- raster("${INPUT_OBS}")
mrf <- raster("${INPUT_MODEL}")
crat=$Threshold ######### Enter Rainfall threshold here (RF in mm/day units)
orft<-orf
mrft<-mrf
png(file="${figname}",width = 520, height = 520, units = "px", pointsize = 12, bg = "white",  res = NA)

layout(matrix(c(1, 2, 3, 0), ncol = 4))
##################################################################
library(rgdal)
shp <- readOGR(dsn = file.path("$shpfile"), stringsAsFactors = F)
mycols=c("azure","darkturquoise","darkseagreen1","darkgreen","yellow","goldenrod1","violetred1","maroon4")
myrange=colorRampPalette(mycols)
mycrange=myrange(${Rainlimt})
mycols = mycrange[c(0:${Rainlimt})]

par(mfrow=c(2,2), family="serif", mai=c(0.5,0.5,0.5,0.5))
plota<-raster::plot(orf,
	xlim = c(${lon1},${lon2}),
	ylim = c(${lat1},${lat2}),
	zlim=c(0,${Rainlimt}),
	col=mycols,
	xlab = "", 
	ylab = "",
	interpolate=TRUE,
	xaxs="i",
	yaxs="i",
	asp=1.00,
	legend=FALSE,
	main="Observation")
#contour(orf,levels=${Threshold},add=TRUE, lwd = 3,lty = "solid")
title(ylab="Latitude",mgp=c(2,1,0),cex.lab=1)
title(xlab="Longitude",mgp=c(2,1,0),cex.lab=1)
plot(shp,add=TRUE)

par(mai=c(0.5,0.5,0.5,0.5))
pdata=read.delim("${OUT_PREFIX}_ScatterData.txt", header = TRUE, sep = "")
pdata[is.na(pdata)] = 0
plot(x = pdata[,1],y = pdata[,2],
    xlim = c(${Threshold},${RainlimtScatter}),
    ylim = c(${Threshold},${RainlimtScatter}),
    ylab="",
    xlab="")
#   main = "Obs vs Model")

abline(coef = c(0,1),lwd = 2)
title(ylab="Model_Shifted",mgp=c(2,1,0),cex.lab=1)
title(xlab="Observation",mgp=c(2,1,0),cex.lab=1)
#####################################################
par(mai=c(0.5,0.5,0.5,0.5))
plotb <-raster::plot(mrf,
	xlim = c(${lon1},${lon2}),
	ylim = c(${lat1},${lat2}),
	zlim=c(0,${Rainlimt}),
	col=mycols,
	xlab = "", 
	ylab = "",
	interpolate=TRUE,
	xaxs="i",
	yaxs="i",
	asp=1.00,
	legend=FALSE,
	main="Model")
#contour(mrf,levels=${Threshold},add=TRUE, lwd = 3,lty = "solid")
title(ylab="Latitude",mgp=c(2,1,0),cex.lab=1)
title(xlab="Longitude",mgp=c(2,1,0),cex.lab=1)
plot(shp,add=TRUE)
											
########
detach("package:raster")
library(gridExtra)
library(grid)
library(formattable)

bg_params=list(fill=c(NA, rep("white",5)), col="white")
tt=ttheme_minimal(core=list(fg_params=list(hjust=0, x=0.10)),
                      rowhead=list(fg_params=list(hjust=0, x=0)))

rows<-c("CRA Threshold=${CRA_threshold}         ")
pushViewport(viewport(x=0.64,y=0.480,height=0.8,width=0.9))
grid.table(rows,theme=tt)

testdf1<-data.frame(Observation=c(${anpts},${aavg},${amax},${avol}),Forecast=c(${bnpts},${bavg},${bmax},${bvol}))
pushViewport(viewport(x=0.58,y=0.37,height=0.8,width=0.5))
rownames(testdf1)<-c("gridpoints >= ${Threshold} mm/day","Average Rainfall (mm/day)","Maximum Rain (mm/day)","Rain Volume (km^3)")
grid.table(testdf1,theme=tt)


testdf2<-data.frame(Original=c(${Orms},${Ocorr}),Shifted=c(${Nrms},${Ncorr}))
rownames(testdf2)<-c("RMS Error (mm/day)","Correleation Coefficient")
pushViewport(viewport(x=0.365,y=0.28,height=1.0,width=1.0))
grid.table(testdf2,theme=tt)


rows2<-c(substitute(paste(bold("Error Decomposition:"))))
pushViewport(viewport(x=0.195,y=0.38,height=0.5,width=0.9))
grid.table(rows2,theme=tt)

testdf3<-data.frame(Original=c(10,7),Shifted=c(8,6))
rows3<-c("Displacement Error","${disErr}%")
rows4<-c("Volume Error", "${volErr}%")
rows5<-c("Pattern Error","${patternErr}%")
ab=rbind(rows3,rows4,rows5)
pushViewport(viewport(x=0.605,y=0.25,height=0.9,width=1.0))
grid.table(ab,theme=tt,rows=NULL)

EOF

Rscript ${date}_PlotCRA_${case}_${ObjNum}.R
cp ${date}_PlotCRA_${case}_${ObjNum}.R ${OUT_PATH}
mv ${figname} ${OUT_PATH}

#done
exit
