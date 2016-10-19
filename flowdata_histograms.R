########Envt
library("plyr", lib.loc="D:/Programs/R/R-3.3.1/library")
library("flowCore", lib.loc="D:/Programs/R/R-3.3.1/library")
################################################################
# Generate histograms of sorted buckets
################################################################
#Read data file
mydir <- 'D:\\Projects\\FlowCytoData\\data'
#datafiles <-c('10102016_Rainbow6a_003.fcs','10102016_Rainbow6a_004.fcs','10102016_Rainbow6a_005.fcs')
datafiles <-c('16092016_rainbow med_002.fcs')


############################# Functions ########################

CLOCKSTEP <- 0.017624999999999998
DISTANCESTEP <- 0.28199999999999997
DISTANCESPAN <- 1155.0719999999999

#Get time elapsed
elapsed <- function(calctime, starttime){
  elapsed <- as.integer(calctime-starttime)
  return(elapsed)
}

#Get drop phase slice
dropphase <- function(dp){
  hex <- 0xF
  dropphase <- 0
  if (dp > 0){
    dropphase <- bitwAnd(dp, hex)
  }
  return(dropphase)
}

#Get previous distance
prevdist <- function(dropPhaseBits){
  prevdist <- bitwAnd(bitwShiftL(dropPhaseBits, 4), 0xfff)
  return(prevdist)
}

#Classifier sorts to a bucket number
bucketnumber <- function(classif){
  hex <- 0x3f
  bucketnumber <- NaN
  if (classif > 0){
    c <- bitwAnd(classif, hex)
    if (c > 0){
      bucketnumber <- log2(c) + 1
    }
    
  }
  return(bucketnumber)
  
}

# sort status is unsorted, sorted or sort abort
sortstatus <- function(classifier,sortresult, dtime){
  if (is.nan(sortresult) && is.nan(classifier)){
    sortstatus <- 'unsorted'
    
  }else if (classifier > 1 && is.nan(sortresult)){
    if (dtime <= 10){
      sortstatus <- 'elect abort'
    }else{
      sortstatus <- 'sort abort'
    }
    
  }else{
    sortstatus <- 'sorted'
  }
  return(sortstatus)
}

#Calculate bucket swing for whole column - either classifier or sortresult
#- between two sequential and sorted events
diffbucket <- function(blist){
  if (class(blist) =='list'){
    blist <- unlist(blist)
  }
  diffbucket <- list()
  for (i in 2:length(blist)){
    a<-blist[i]
    b<-blist[i-1]
    if (!is.nan(a) && !is.nan(b)){
      bindiff <- a-b
      if (abs(bindiff) > 6){
        print(paste("i: ",i, "=", bindiff))
      }
    }else{
      bindiff <- NaN
    }
    
    print(paste("a=",a," b=",b,"a-b=",a-b))
    diffbucket[i] <- bindiff
    
    
  }
  print(length(diffbucket))
  return(diffbucket)
}

########################## MAIN ####################################
#Produce plots
par(mfrow = c(2, 2))
par(cex = 0.6)
par(mar = c(2,2,0,0), oma=c(1,1,1,1))

##################
#for (da in datafiles){
da = '16092016_rainbow med_002.fcs'
testdata1 <- file.path(mydir,da)
x1 = read.FCS(testdata1,alter.names=TRUE,transformation = FALSE)
summary(x1)
#Determine filtering params
numsortways = keyword(x1,c('NUMSORTWAYS'))
numbins = as.integer(unlist(numsortways, use.names = FALSE))
df <- data.frame(exprs(x1)[,1:10])
#Generate delta time and elapsed time
df$TimeCalc <-apply(df[,c('Time.1','Time.2','Time.3')], 1,function(d){0.017625*(((d[3]*65536*65536)+(d[2]*65536)+d[1]))})
deltaTime <-diff(df$TimeCalc)
df$DeltaTime <-c(0,deltaTime)
starttime <- df$TimeCalc[1]
df$ElapsedTime <-mapply(elapsed,df$TimeCalc,starttime)
#Generate bucket numbers and sort status
df$ClassBucket <-mapply(bucketnumber,df$Classifier.Bits)
s<-unlist(df$ClassBucket)
dfres <-hist(s,freq=TRUE,main=paste("All Bucket counts"), col='yellow', plot=TRUE,xlab="Bucket num",xlim=range(0,6), ylim=range(0,3000))
text(dfres$mids, dfres$density, dfres$counts, adj = c(.5, -.5),col = "blue")
deltaBucketAll <- diffbucket(s)

df$SortBucket <-mapply(bucketnumber,df$Sort.Result.Bits)
s<-unlist(df$SortBucket)
dfres1<-hist(s, freq=TRUE,main=paste("Sorted Bucket counts"), col='yellow', plot=TRUE,xlab="Bucket num",xlim=range(0,6), ylim=range(0,3000))
text(dfres1$mids, dfres1$density, dfres1$counts, adj = c(.5, -.5), col = "blue")
deltaBucketSorted <- diffbucket(s)

dfres2<-hist(!is.nan(deltaBucketAll), freq=TRUE,main=paste("Swing Bucket all"), col='red', plot=TRUE,xlab="Bucket num")
text(dfres2$mids, dfres2$density, dfres2$counts, adj = c(.5, -.5), col = "blue")
dfres3<-hist(!is.nan(deltaBucketSorted), freq=TRUE,main=paste("Swing Bucket sorted"), col='red', plot=TRUE,xlab="Bucket num")
text(dfres3$mids, dfres3$density, dfres3$counts, adj = c(.5, -.5), col = "blue")



#}



