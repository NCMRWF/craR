### Read Object files identified by feature finder ###
### and does CRA Analysis                          ###
CRA_Err_Decomp<-function(ActOrf,ActMrf,orf,mrf,thr,objnum,out_prefix,Ngrids){

library(raster)
#library(fields)
##
############# Functions for calculation of RMSE, COG and Statistics ########
RMSE <- function(x, y) { sqrt(mean((x - y)^2)) }
#######
cog <- function(r){
  ptop<-rasterToPoints(r,na.rm=TRUE)
  xpt<-ptop[,1] ; ypt<-ptop[,2] ; rf<-ptop[,3]
  xp<-sum(xpt*rf)/sum(rf)
  yp<-sum(ypt*rf)/sum(rf)
  return(c(xp,yp))
}
#######
objstats<-function(rf){
  xgrd=xres(rf)
  ltmin<-rf@extent[3]
  ltmax<-rf@extent[4]
ptop<-rasterToPoints(rf,na.rm=TRUE)
orfpts<-ptop[,3] ;nptso<-length(ptop[,3])
rfmean<-mean(orfpts) ; rfmax<-max(orfpts) # Getting mean & max from stats
pi=22/7
cosfac=cos(0.5*(ltmin+ltmax)*pi/180.)
factor=xgrd*xgrd*111.*(cosfac*111.)*1.e-6    #to convert to rain volume [km^3]
rfvol=nptso*rfmean*factor
rfarea=(nptso*factor)*1.e6
return(c(nptso, rfmean, rfmax, rfvol, rfarea))
}
###

#### Reading input file ###
ActOrf<-raster(ActOrf)  ## Actual Input Obs
ActMrf<-raster(ActMrf)        ## Actual input model
orf<-raster(orf)        ## Feature from Obs
mrf<-raster(mrf)              ## Feature from Model

##################################################################
xgrd=xres(orf); ygrd=yres(orf); grd=xgrd

CRAmask<-(sum(orf,mrf,na.rm=TRUE))           # get the overlap of model and observation
CRAmask[is.na(orf) & is.na(mrf)]<-NA
CRAmask_orig<-CRAmask      

######### Calulate all stats over the identified feature  ##########
orfstatsThr<-objstats(orf)
mrfstatsThr<-objstats(mrf)
OrfThrnpts<-orfstatsThr[1] ; MrfThrnpts<-mrfstatsThr[1]
OrfThravg<-orfstatsThr[2] ; MrfThravg<-mrfstatsThr[2]
OrfThrmax<-orfstatsThr[3] ; MrfThrmax<-mrfstatsThr[3]
OrfThrvol<-orfstatsThr[4] ; MrfThrvol<-mrfstatsThr[4]
OrfThrarea<-orfstatsThr[5] ; MrfThrarea<-mrfstatsThr[5]  

############## Apply CRA area mask on Obs and Model ##################
a<-ActOrf ; b<-ActMrf
a[is.na(CRAmask)]<-NA ;b[is.na(CRAmask)]<-NA
####   Calculate COG for Obs and Model over CRA Area ###
COGorf<-cog(a) ; COGmrf<-cog(b)
## Calulate all stats over the CRA masked area ################
astats<-objstats(a) ;bstats<-objstats(b)
anpts<-astats[1] ; bnpts<-bstats[1]
aavg<-astats[2] ; bavg<-bstats[2]
amax<-astats[3] ; bmax<-bstats[3]
avol<-astats[4] ; bvol<-bstats[4]
aarea<-astats[5] ; barea<-bstats[5]  
apts=cbind(xyFromCell(a, 1:ncell(a)), values(a)) 
apts<-apts[,3]
bpts=cbind(xyFromCell(b, 1:ncell(b)), values(b)) 
bpts<-bpts[,3]
numapts<-astats[1]

##### Calculate RMSE & Corr before shifting ######
apts[is.na(bpts)]<-NA;bpts[is.na(apts)]<-NA;   
orig_cor<-cor(na.omit(apts),na.omit(bpts)) ; old_cor<-orig_cor    
orig_rms<-RMSE(na.omit(apts),na.omit(bpts)) ; old_rms<-orig_rms
################################################
print("RMS Error Minimization Starts....")
igrd=Ngrids   # MAX NUMBER OF GRIDS TO SCAN IN W TO E & N TO S ; 169 SCANS
idx=0
prev_rms=99999
for (c in -igrd:igrd){
 dxgrd=c*grd
  for (r in -igrd:igrd){
   idx=idx+1
   dygrd=r*grd
   newb<-shift(ActMrf,dx=dxgrd, dy=dygrd)
   newCRAmask_b<-shift(b,dx=dxgrd, dy=dygrd);  

#### Updation of verification mask based on shifted Model field #########
   verif_mask<-merge(newCRAmask_b, CRAmask_orig, tolerance = 0.5)
#   plot(verif_mask)
##
   av<-ActOrf ; abv<-ActMrf ; bv<-newb
   av[is.na(verif_mask)]<-NA ; abv[is.na(verif_mask)]<-NA 
   bv[is.na(verif_mask)]<-NA
   ptav=cbind(xyFromCell(av, 1:ncell(av)), values(av))
   ptabv=cbind(xyFromCell(abv, 1:ncell(abv)), values(abv))
   ptbv=cbind(xyFromCell(bv, 1:ncell(bv)), values(bv))
   avpts<-ptav[,3] ; abvpts<-ptabv[,3] ; bvpts<-ptbv[,3]
   lnavpts<-length(avpts)
   chkpts<-(lnavpts/numapts)*100
   avpts[is.na(bvpts)]<-NA; bvpts[is.na(avpts)]<-NA;
## Calculate the correlation and RMSE after shifting the model field ##
   shift_cor<-cor(na.omit(avpts),na.omit(bvpts)) ; 
   shift_rms<-RMSE(na.omit(avpts),na.omit(bvpts))
   shift_rms[is.na(shift_rms)]<-99999

##### Calculates COG,Corr and CRA Error Components if shift_rms < old_rms (Opt1) ########
    if((shift_rms < prev_rms)) {
     print(shift_rms)
     print(shift_cor)
	iopt=1 # 1. Sq. Err minimization 2. for patcorr maximization

####------ for iopt=2-----------#########
#    if(shift_cor > old_cor){
#     iopt=2 # 1. Sq. Err minimization 2. for patcorr maximization
#    print('Pattern Correl')  
#    print(shift_cor)   
#    print(old_cor)   
#	print(c) 
#	print(r) 
#########################################
     prev_rms<-shift_rms
     avstats<-objstats(av);bvstats<-objstats(abv)
     avnpts<-avstats[1] ; bvnpts<-bvstats[1]
     avavg<-avstats[2] ; bvavg<-bvstats[2]
     avmax<-avstats[3] ; bvmax<-bvstats[3]
     avvol<-avstats[4] ; bvvol<-bvstats[4]
     avarea<-avstats[5] ; bvarea<-bvstats[5]  
     mrfstats_shift<-summary(bvpts) ; mrfmax_shift<-mrfstats_shift[6]   
     mrfmax=bmax 
     
     COGmrf_shift<-cog(bv)
     pattcorr=shift_cor ; RMSshift=shift_rms
     avpts[is.na(bvpts)]<-NA; bvpts[is.na(avpts)]<-NA;
     ll<-avpts ; mm<-abvpts ; kk<-bvpts
     difr=(na.omit(kk)-na.omit(ll)) ; toterrshift=sum(difr*difr)
     avga=mean(na.omit(ll)) ; avgb=mean(na.omit(mm)) ; avgbshift=mean(na.omit(kk))
     sda=sqrt(mean((na.omit(ll)-avga)*(na.omit(ll)-avga)))
     sdb=sqrt(mean((na.omit(mm)-avgb)*(na.omit(mm)-avgb)))  
     sdbshift=sqrt(mean((na.omit(kk)-avgbshift)*(na.omit(kk)-avgbshift)))
      if(chkpts < 50.){
       print('MORE THAN 50% OF POINTS SHIFTED OUT') ; outcode1=1
      }else{
       print('MORE THAN 50% OF POINTS USED') ; outcode1=0
      }
      if(mrfmax_shift != mrfmax){
       print('MAX RAINFALL LOST') ; outcode2=1
      }else{
       print('MAX RAINFALL RETAINED') ; outcode2=0
      }
    }
  }

}
print("RMS Error Minimization Ends....")
Observed_Rain<-ll
Forecast_Rain_Shifted<-kk
rfscat<-cbind(Observed_Rain, Forecast_Rain_Shifted)

############## CRA error decomposition & other statistics ##############
stats0<-c(OrfThrnpts,MrfThrnpts,OrfThravg,MrfThravg,OrfThrmax,MrfThrmax,OrfThrvol,MrfThrvol,OrfThrarea,MrfThrarea)
stats0<-t(round(stats0,digits=2))
stats1<-c(anpts,bnpts,aavg,bavg,amax,bmax,avol,bvol,aarea,barea,orig_rms, orig_cor)
stats1<-t(round(stats1,digits=2))
stats2<-c(avnpts,bvnpts,avavg,bvavg,avmax,bvmax,avvol,bvvol,avarea,bvarea,
          orig_rms,RMSshift,orig_cor,shift_cor)
stats2<-t(round(stats2,digits=2))
RMSorig=orig_rms
MSEtotal=RMSorig*RMSorig               #total error in CRA
##
 if(iopt==1){
  print('Successfully completed the CRA Analysis using squared error minimization method')
  MSEshift=RMSshift*RMSshift           #error left after two blobs superimposed
  MSEdisplacement=MSEtotal-MSEshift    #difference in RMSE before & after shift
  MSEvolume=(avgbshift-avga)*(avgbshift-avga)        #mean intensity difference
  MSEpattern=MSEshift-MSEvolume        #residual fine scale pattern difference
##
 }

#### Error decomposition for Opt2 (Maximizing Pattern Correlation Method)
 if(iopt==2){
  print('Successfully completed the CRA Analysis using Pattern Correlation Method')
  pattcorr=shift_cor
  corrorig=orig_cor
  MSEdisplacement=2.*sda*sdb*(pattcorr-corrorig)        #uses correlation before and after shift
  MSEvolume=(avgb-avga)*(avgb-avga) # mean intensity difference (note this is not identical to MSEvolume for minERR)
  MSEpattern=(sda-sdb)^2. + 2.*sda*sdb*(1.-pattcorr) # conditional bias + imperfect corr. terms
 }
####### Convert CRA Error components to percent ###########
MSEdisplacement=MSEdisplacement/MSEtotal*100.
MSEvolume      =MSEvolume      /MSEtotal*100.
MSEpattern     =MSEpattern     /MSEtotal*100.
MSEcomps<-c(MSEtotal,MSEdisplacement,MSEvolume,MSEpattern)
MSEcomps<-t(round(MSEcomps,digits=2))
##
cogs<-c(COGorf[1],COGorf[2],COGmrf[1],COGmrf[2],COGmrf_shift[1],COGmrf_shift[2])
cogs<-t(round(cogs,digits=2))
xdif_orig=abs(COGorf[1]-COGmrf[1])*111.1 ; ydif_orig=abs(COGorf[2]-COGmrf[2])*111.1
xdif_shift=abs(COGorf[1]-COGmrf_shift[1])*111.1 ; ydif_shift=abs(COGorf[2]-COGmrf_shift[2])*111.1
dist_orig=sqrt((xdif_orig*xdif_orig)+(ydif_orig*ydif_orig))
dist_shift=sqrt((xdif_shift*xdif_shift)+(ydif_shift*ydif_shift))
xydistcog<-c(xdif_orig,ydif_orig,xdif_shift,ydif_shift,dist_orig,dist_shift)
xydistcog<-t(round(xydistcog,digits=2))
######################### Writing Statistics to files ################################################

write.table(rfscat,file=paste(out_prefix,objnum,'ScatterData.txt',sep="_"),sep="   ",col.names=T,row.names=F,append=F)
write.table(stats0,file=paste(out_prefix,objnum,'StatsOrig.txt',sep="_"),sep="   ",col.names=F,row.names=F,append=F)
#write.table(stats1,file=paste(out_prefix,objnum,'rfstats1.txt',sep="_"),sep="   ",col.names=F,row.names=F,append=F)
write.table(stats2,file=paste(out_prefix,objnum,'AllStats.txt',sep="_"),sep="   ",col.names=F,row.names=F,append=F)
write.table(MSEcomps,file=paste(out_prefix,objnum,'cra_decomp.txt',sep="_"),sep="   ",col.names=F,row.names=F,append=F)
write.table(cogs,file=paste(out_prefix,objnum,'cog.txt',sep="_"),sep="   ",col.names=F,row.names=F,append=F)
write.table(xydistcog,file=paste(out_prefix,objnum,'xydistcog.txt',sep="_"),sep="   ",col.names=F,row.names=F,append=F)
}
