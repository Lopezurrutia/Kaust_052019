## rename 's/Rdara/Rdata/' *
setwd("/home/gelu/saudiprochl/refsgrowth/ciclo24/")
load("miarray_datos24_Fig07.Rdata")
gg<-read.table('datos24_Fig07.csv', header=TRUE, sep=',')


## remove outlier
outlier<-which(dimnames(miarray_filtrado)[[3]]=="845")
miarray_filtrado<-miarray_filtrado[,,-outlier]
## remove replicate without beads
checin<-which(dimnames(miarray_filtrado)[[3]]=="1027")
miarray_filtrado[,,checin]<-NA
miarray_filtrado<-miarray_filtrado[,,-checin]

## remove those two samples from the analysis
gg<-gg[-which(gg$FileName==1027),]
gg<-gg[-which(gg$FileName==845),]


## Get decimal hour from hour in text format
dechour<-function(x)
    {
        x<-as.numeric(strsplit(as.character(x),':')[[1]])
        x<-x[1]+x[2]/60
        return(x)
    }
dechora<- apply(as.matrix(gg$hora),1,dechour)

## Rearrange hours so day starts at 5:30
dechora[which(dechora<5)]<-dechora[which(dechora<5)]+24

## ## Rearrange hours so day starts at 9:30
## dechora[which(dechora<9)]<-dechora[which(dechora<9)]+24


## Change replica numbers so all temperatures have replicates 1,2,3
gg$replica[gg$replica==4]<-2
gg$replica[gg$replica==7]<-3

## Arrange data
RepMat<-tapply(gg$FileName,list(gg$temperatura, dechora,gg$replica),mean,na.rm=TRUE)
dim(RepMat) ## 3 temperatures x 12 times x 3 replicates

populations<-c("MIT9301.G1","MIT9301.S","MIT9301.G2")
volbreaks<-seq(-24,2.6,by=0.2)
## volbreaks<-seq(-0.2,13,by=0.05) ## for lineal scale


## Replay this after running loop to adjust breaks to breaks with data
## ValidBins<-apply(histvols,c(5),sum,na.rm=TRUE)
## volbreaks[which(ValidBins>0)]

midbreaks<-(volbreaks[2:length(volbreaks)]+volbreaks[1:(length(volbreaks)-1)])/2
histvols<-array(NA,c(dim(RepMat),length(populations),length(volbreaks)-1))
dim(histvols)
## 3 temperatures x 12 times x 3 replicates x 3 cycle phases x 133 volume bins







## For each file calculate the SSC histograms with frequencies in cells/mL for each population
for (i in 1:dim(histvols)[1]){
    for (j in 1:dim(histvols)[2]){
        for (h in 1:dim(histvols)[3]){
            fcsnumber<-RepMat[i,j,h]
            if (!is.na(fcsnumber)){ ## Skip the lost samples
                print(fcsnumber)
                load(paste("fcs_gated_ciclo24/RdataGated_",fcsnumber,".Rdata",sep=""))
                dilution_factor<-gg[which(gg$FileName==fcsnumber),]$d3
                SSCbeads<-FCSGate$MatFlow$SSC.H[FCSGate$MatFlow$Gates=="beads"]
                for ( pop in 1:length(populations)){
                    population<-populations[pop]
                    SSCG1<-FCSGate$MatFlow$SSC.H[FCSGate$MatFlow$Gates==population]-mean(SSCbeads)
                    DiaG1<-SSCG1*0.874+1.62
                    VolG1<-(4/3)*pi*((DiaG1/2)^3)
                    ## CARE: mean(VolG1) is different than (4/3)*pi*((mean(DiaG1)/2)^3), Antonio uses the latter, the former is more correct
                    freq<-hist(log(VolG1),plot=FALSE,breaks=volbreaks)
                    ## print(range(VolG1))
                    ## freq<-hist((VolG1),plot=TRUE,breaks=volbreaks)
                    histvols[i,j,h,pop,]<-freq$counts*dilution_factor
                }
            }
        }
    }
}


dimnames(histvols)<-list(c(19,26,30),sort(unique(dechora)),seq(1,3),c("G1","S","G2"),midbreaks)





histvols2<-histvols
histvols2[which(histvols2<=0)]<-NA

getnewg1<-function(timesel="17.5",
                   tempsel=3,replicate=1,nexttime=2,phase=1)
{
    selecttime<-which(dimnames(histvols2)[2][[1]]==timesel)
    time1<-histvols2[tempsel,selecttime,replicate,phase,]
    time2<-histvols2[tempsel,selecttime+nexttime,replicate,phase,]
    ylims<-range(c(time1,time2),na.rm=TRUE)
    plot(time1~midbreaks,type="l",ylim=ylims)
    points(time2~midbreaks,col="red",type="l")
    newg1<-time2-time1
    newg1<-newg1/sum(newg1,na.rm=TRUE)
    newg1[newg1<0]<-NA
    totvol<-sum(newg1*midbreaks,na.rm=TRUE)
    totcells<-sum(newg1,na.rm=TRUE)
    meanvol<-exp(totvol/totcells)

    return(list(freqg1=newg1,meanvol=meanvol))
    }

### G1
newg1<-array(NA,dim=c(3,3))
temps<-array(NA,dim=c(3,3))
temps[1,]<-c(19,19,19)
temps[2,]<-c(26,26,26)
temps[3,]<-c(30,30,30)


par(mfrow=c(3,3))
for (replicate in seq(1,3)){
    newg1[1,replicate]<-getnewg1(timesel="17.5",tempsel=1,replicate=replicate,nexttime=2,phase=1)$meanvol
    newg1[2,replicate]<-getnewg1(timesel="17.5",tempsel=2,replicate=replicate,nexttime=2,phase=1)$meanvol
    newg1[3,replicate]<-getnewg1(timesel="17.5",tempsel=3,replicate=replicate,nexttime=2,phase=1)$meanvol
}
plot(newg1~temps)

xs<-c(19,26,30)
ys<-apply(newg1,1,mean)
sdev<-apply(newg1,1,sd)
micolor<-"black"
points(ys~xs,lwd=2,pch=19,type="b",col=micolor)
arrows(xs, ys-sdev, xs, ys+sdev, length=0.05, angle=90, code=3,col=micolor)## sd bars



### G2

newg1<-array(NA,dim=c(3,3))
temps<-array(NA,dim=c(3,3))
temps[1,]<-c(19,19,19)
temps[2,]<-c(26,26,26)
temps[3,]<-c(30,30,30)

par(mfrow=c(3,3))
for (replicate in seq(1,3)){
    newg1[1,replicate]<-getnewg1(timesel="21.5",tempsel=1,replicate=replicate,nexttime=-2,phase=3)$meanvol
    newg1[2,replicate]<-getnewg1(timesel="21.5",tempsel=2,replicate=replicate,nexttime=-2,phase=3)$meanvol
    newg1[3,replicate]<-getnewg1(timesel="21.5",tempsel=3,replicate=replicate,nexttime=-2,phase=3)$meanvol
}
plot(newg1~temps)

xs<-c(19,26,30)
ys<-apply(newg1,1,mean)
sdev<-apply(newg1,1,sd)
micolor<-"black"
points(ys~xs,lwd=2,pch=19,type="b",col=micolor)
arrows(xs, ys-sdev, xs, ys+sdev, length=0.05, angle=90, code=3,col=micolor)## sd bars










## calculate the total log volume of the particles in each replicate x cycle phase
PhaseLogVol<-  sweep(histvols, 5, midbreaks, FUN = "*")
PhaseLogVol<-  apply(PhaseLogVol,c(1,2,3,4),sum)

## calculate the total number of the particles in each replicate x cycle phase
PhaseTotalCells<-apply(histvols,c(1,2,3,4),sum)

## The mean( geometric mean) volume of each replicate x cycle phase is
PhaseMeanVol<-exp(PhaseLogVol/PhaseTotalCells)

##For the total population of Prochl we have to add volumes and densities before mean
dim(PhaseTotalCells)
ProChlLogVol<-apply(PhaseLogVol,c(1,2,3),sum)
ProChlTotalCells<-apply(PhaseTotalCells,c(1,2,3),sum)
ProChlMeanVol<-exp(ProChlLogVol/ProChlTotalCells)

## Now the mean of the three replicates
sizestrue<-apply(ProChlMeanVol,c(1,2),mean,na.rm=TRUE)

##  Repeat total array for each pahse to make division
RepProchlTotal<-rep(ProChlTotalCells,3)
dim(RepProchlTotal)<-dim(PhaseTotalCells)
PercentPhase<-PhaseTotalCells/RepProchlTotal
range(PercentPhase,na.rm=TRUE)


ThreeTemps<-function(DataMatrix,DataMatrix2=NULL,phases=c(1,3),temperatures=seq(1,3),ylims=c(0,1),ylab="",plotTot=TRUE,plotPhase=TRUE){

    MeanPer<-apply(DataMatrix,c(1,2,4),mean,na.rm=TRUE)
    SDPer<-apply(DataMatrix,c(1,2,4),sd,na.rm=TRUE)
    ## ## Calculate total of all phases
    ## SumPer<-apply(DataMatrix,c(1,2,3),sum,na.rm=TRUE)
    ## SumPer[SumPer==0]<-NA
    if (!is.null(DataMatrix2)){
    SDTot<-apply(DataMatrix2,c(1,2),sd,na.rm=TRUE)
    SumTot<-apply(DataMatrix2,c(1,2),mean,na.rm=TRUE)
    }
    xpos<-as.numeric(colnames(MeanPer))
    labelshour<-paste(floor(xpos),":30",sep="")

    
    plot(0~1,ylim=ylims,xlim=c(5.5,27.5),xlab="Time",ylab=ylab,las=1,,xaxt="n",cex.lab=1.5)
    rect(16.8, ylims[1]-0.1,29,ylims[2]+0.1 , col = gray(0.95),lwd=0)
    text(10,ylims[2],"Light")
    text(23,ylims[2],"Dark")
    axis(side=1,at=xpos,labels=labelshour) ## Axis with time of day and hour labels
    box()
    pchs<-c(15,17,19)
    phaselabel=c("G1","S","G2")
    legendsx<-c(5.5,16.5,27.5)
    for (phasen in phases){
        mipch<-pchs[phasen]
        colordec<-1.33-(phasen/3)
        colorstemp<-c("blue","green","red")
        for (i in temperatures){ ##for each temperature
            templabel<-c("19","26","30")
            curvelabel<-paste(phaselabel[phasen],templabel[i],sep="-")
            if(plotPhase){
                micolor<-rep(0,3)
                micolor[i]<-colordec
                xs<-xpos
                ys<-MeanPer[i,,phasen]
                sdev<-SDPer[i,,phasen]
                micolor<-rgb(micolor[3],micolor[2],micolor[1])
                points(ys~xs,lwd=2,pch=mipch,type="b",col=micolor)
                legend(legendsx[i],ylims[1],curvelabel,col=micolor,pch=mipch,border="")
                arrows(xs, ys-sdev, xs, ys+sdev, length=0.05, angle=90, code=3,col=micolor)## sd bars
            }
            if(!is.null(DataMatrix2)){
                xs<-xpos
                ys<-SumTot[i,]
                sdev<-SDTot[i,]
                micolor<-colorstemp[i]
                points(ys~xs,lwd=2,pch=mipch,type="b",col=micolor)
                arrows(xs, ys-sdev, xs, ys+sdev, length=0.05, angle=90, code=3,col=micolor)## sd bars
            }

        }
    }
    
}




ThreeTemps(DataMatrix=PercentPhase,phases=c(1),temperatures=c(1,2,3),ylab="Percent abundance",ylim=c(0.1,0.9))





png("ExampleCycle.png",width=600,height=900)
par(mfrow=c(2,1))
par(mar=c(4.1,5,1,1))


ThreeTemps(DataMatrix=PercentPhase,phases=c(1,2,3),temperatures=c(1),ylab="Percent abundance",ylim=c(0.1,0.9))
ThreeTemps(DataMatrix=PhaseMeanVol,DataMatrix2=ProChlMeanVol,phases=c(1,2,3),temperatures=c(1),ylims=c(0,0.2),ylab=expression("Cell volume ("*mu*"m"^{3}*")"),plotPhase=FALSE)
dev.off()



png("ExampleCycle2.png",width=600,height=900)
par(mfrow=c(2,1))
par(mar=c(4.1,5,1,1))


ThreeTemps(DataMatrix=PercentPhase,phases=c(1,2,3),temperatures=c(1),ylab="Percent abundance",ylim=c(0.1,0.9))
ThreeTemps(DataMatrix=PhaseMeanVol,DataMatrix2=ProChlMeanVol,phases=c(1),temperatures=c(1),ylims=c(0,0.2),ylab=expression("Cell volume ("*mu*"m"^{3}*")"),plotPhase=TRUE)
dev.off()



png("ExampleCycle3.png",width=600,height=900)
par(mfrow=c(2,1))
par(mar=c(4.1,5,1,1))


ThreeTemps(DataMatrix=PercentPhase,phases=c(1),temperatures=c(1),ylab="Percent abundance",ylim=c(0.1,0.9))
ThreeTemps(DataMatrix=PhaseMeanVol,DataMatrix2=ProChlMeanVol,phases=c(1,2,3),temperatures=c(1),ylims=c(0,0.2),ylab=expression("Cell volume ("*mu*"m"^{3}*")"),plotPhase=FALSE)
dev.off()



png("ExampleCycle4.png",width=600,height=900)
par(mfrow=c(2,1))
par(mar=c(4.1,5,1,1))


ThreeTemps(DataMatrix=PercentPhase,phases=c(1),temperatures=c(1,2,3),ylab="Percent abundance",ylim=c(0.1,0.9))
ThreeTemps(DataMatrix=PhaseMeanVol,DataMatrix2=ProChlMeanVol,phases=c(1,2,3),temperatures=c(1),ylims=c(0,0.2),ylab=expression("Cell volume ("*mu*"m"^{3}*")"),plotPhase=FALSE)
dev.off()


png("ExampleCycle5.png",width=600,height=900)
par(mfrow=c(2,1))
par(mar=c(4.1,5,1,1))


ThreeTemps(DataMatrix=PercentPhase,phases=c(1),temperatures=c(1,2,3),ylab="Percent abundance",ylim=c(0.1,0.9))
ThreeTemps(DataMatrix=PhaseMeanVol,DataMatrix2=ProChlMeanVol,phases=c(1,2,3),temperatures=c(1,2,3),ylims=c(0,0.2),ylab=expression("Cell volume ("*mu*"m"^{3}*")"),plotPhase=FALSE)
dev.off()




png("Ciclos24.png",width=900,height=600)
par(mfcol=c(4,3))
for (i in seq(1,3)){
ThreeTemps(DataMatrix=PercentPhase,phases=c(1,2,3),temperatures=c(i),ylab="Percent abundance")
ThreeTemps(DataMatrix=PhaseMeanVol,DataMatrix2=ProChlMeanVol,phases=c(1,2,3),temperatures=c(i),ylims=c(0,0.5),ylab=expression("Cell volume ("*mu*"m"^{3}*")"))
ylims<-range(PhaseTotalCells[,,,3],na.rm=TRUE)
ThreeTemps(DataMatrix=PhaseTotalCells,phases=c(1,2,3),temperatures=c(i),ylab="Percent abundance",ylims=ylims)
ylims<-range(PhaseTotalCells[,,,1],na.rm=TRUE)
ThreeTemps(DataMatrix=PhaseTotalCells,phases=c(1,2,3),temperatures=c(i),ylab="Percent abundance",ylims=ylims)
}
dev.off()





png("G1s.png",width=600,height=600)

ThreeTemps(DataMatrix=PercentPhase,phases=c(1),temperatures=c(1,2,3),ylab="Percent abundance")


dev.off()



ThreeTemps(DataMatrix=ProChlMeanVol,phases=c(1),temperatures=c(1,2,3),ylab="Percent abundance",ylim=c(0.1,0.9))



ThreeTemps(DataMatrix=PhaseMeanVol,phases=c(3),temperatures=c(1,2,3),ylab="Percent abundance",ylim=c(0.1,0.5))





plot(histvols[2,4,1,1,]~midbreaks,xlim=c(-5,0))



## ## COMPARE TO ANTONIOS MATRIX
## indexpops<-c(7,4,6,5)
## pops<-c("All","G1","S","G2")
## ##abundance
## variable<-35
## varname<-dimnames(miarray_filtrado)[[2]][variable]
## ## Counts match perfectly
## i<-1
## sizes<-tapply(miarray_filtrado[indexpops[i],variable,],list(gg$temperatura, dechora,gg$replica),mean,na.rm=TRUE)

## plot(sizes[1,,1])
## points(ProChlTotalCells[1,,1],pch=19,col="blue")

## plot(sizes~ProChlTotalCells)
## abline(0,1)

## ##volume is overestimated by mean instead of geometric mean
## variable<-38
## varname<-dimnames(miarray_filtrado)[[2]][variable]
## i<-1
## sizes<-tapply(miarray_filtrado[indexpops[i],variable,],list(gg$temperatura, dechora,gg$replica),mean,na.rm=TRUE)
## onetemp<-sizes[2,,]
## plot(sizes[2,,])

## plot(sizes~ProChlMeanVol)
## abline(0,1)



temperature<-3
ylims<-range(log(ProChlTotalCells[temperature,,]),na.rm=TRUE)
plot(log(ProChlTotalCells[temperature,,1]),pch=19,ylim=ylims)
points(log(ProChlTotalCells[temperature,,2]),col="blue",pch=19)
points(log(ProChlTotalCells[temperature,,3]),col="red",pch=19)

