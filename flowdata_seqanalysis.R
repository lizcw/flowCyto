#Envt
library("plyr", lib.loc="D:/Programs/R/R-3.3.1/library")
library("flowCore", lib.loc="D:/Programs/R/R-3.3.1/library")
#Read data file
mydir <- 'D:\\Projects\\FlowCytoData\\data'
#datafiles <-c('10102016_Rainbow6a_003.fcs','10102016_Rainbow6a_004.fcs','10102016_Rainbow6a_005.fcs')
datafiles <-c('16092016_rainbow med_002.fcs')
#Produce plots
par(mfrow = c(length(datafiles), 3))
par(cex = 0.6)
par(mar = c(3,3,0,0), oma=c(1,1,1,1))


#Functions
sortstatus <- function(sortresult,roibits, dtime){ 
  if (roibits == 1 && sortresult == 0){
    sortstatus <- 'unsorted'
  }else if (roibits > 1 && sortresult == 0){
    if (dtime <= 10){
      sortstatus <- 'elect abort'
    }else{
      sortstatus <- 'sort abort'
    }
    
  }else{
    sortstatus <- 'sorted'
  }
}
sortbucket <- function(roibits){ 
  if (roibits > 1){
    sortbucket <- as.integer(log2(roibits-1))
  }else{
    sortbucket <- NA
  }
}
diffbucket <- function(r){
  #diffbucket <- c(0)
  for (i in 2:length(r)){
    
    if (r[i] > 1 && r[i-1] > 1){ #exclude unsorted for both
      bin <- as.integer(log2(r[i]-1))
      binprev <- as.integer(log2(r[i-1]-1))
      bindiff <- binprev-bin
      diffbucket <- c(diffbucket, bindiff)
      if (abs(bindiff) > 6){
        print(paste("i: ",i, "=", bindiff))
      }

    }
    
  }
}
elapsed <- function(calctime, starttime){
  elapsed <- as.integer(calctime-starttime)
}

##MAIN
for (da in datafiles){
  testdata1 <- file.path(mydir,da)
  x1 = read.FCS(testdata1,alter.names=TRUE,transformation = FALSE)
  summary(x1)
  #Determine filtering params
  maxinterval = 15  #calculate from ?dropfreq
  swinglimit = 4    #deltabucket threshold (absolute)
  
  
  #Main
  df <- data.frame(exprs(x1)[,1:10])#or subset for testing[1:10,1:10])
  df$TimeCalc <-apply(df[,c('Time.1','Time.2','Time.3')], 1,function(d){0.017625*(((d[3]*65536*65536)+(d[2]*65536)+d[1]))})
  deltaTime <-diff(df$TimeCalc)
  df$DeltaTime <-c(0,deltaTime)
  starttime <- df$TimeCalc[1]
  df$ElapsedTime <-mapply(elapsed,df$TimeCalc,starttime)
  df$SortBucket <-mapply(sortbucket,df$ROI.Bits.1.16)
  #deltaBucket <- diff(df$SortBucket)
  deltaBucket <- diffbucket(df$ROI.Bits.1.16)
  a <- df$Classifier.Bits
  a1 <- a -32768
  deltaBucket <- diffbucket(a1)
  
  df$DeltaBucket <-c(0,deltaBucket)
  df$SortStatus <-mapply(sortstatus,df$Sort.Result.Bits,df$ROI.Bits.1.16, df$DeltaTime)
  #Show all data
  dfres1 <-hist(df$SortBucket,breaks=6, freq=TRUE,main=paste(da), col='yellow', plot=TRUE,xlab="Bucket diff",xlim=range(0,6))
  text(dfres1$mids, dfres1$density, dfres1$counts, adj = c(.5, -.5), col = "blue")
  dfres <-hist(abs(df$DeltaBucket),breaks=6, freq=TRUE,main=paste("Sequential Bucket swings"), col='yellow', plot=TRUE,xlab="Bucket diff",xlim=range(0,6))
  text(dfres$mids, dfres$density, dfres$counts, adj = c(.5, -.5), col = "blue")
  #Filter on deltatime to detect close sequence events
  shortdf <- subset(df,DeltaTime <= maxinterval)
  res <-hist(abs(shortdf$DeltaBucket),breaks=6, freq=TRUE,col='green', add=TRUE)
  text(res$mids, res$density, res$counts, col = "red")
  res$counts[6]
  #Subset swings of 5 or 6 - view sequence of prev/post 4 swings and/or status aborted
  largeswings <- subset(shortdf[,c('Drop.Phase', 'ROI.Bits.1.16','Sort.Result.Bits','ElapsedTime','DeltaTime','SortBucket','DeltaBucket', 'SortStatus')],abs(DeltaBucket) > swinglimit)
  largeswings
  smallswings <- subset(shortdf[,c('Drop.Phase', 'ROI.Bits.1.16','Sort.Result.Bits','ElapsedTime','DeltaTime','SortBucket','DeltaBucket', 'SortStatus')],abs(DeltaBucket) > 0 & abs(DeltaBucket) <= swinglimit)
  smallswings
  #Compare relative counts of status for each
  xsmall <- count(smallswings$SortStatus)
  xlarge <- count(largeswings$SortStatus)
  resultdata <-matrix(c(xsmall[,2],xlarge[,2]),ncol=2,byrow=FALSE)
  colnames(resultdata) <- c('SmallSwings','LargeSwings')
  rownames(resultdata) <- xsmall[,1]
  resulttable <- as.table(resultdata)
  ptable <-prop.table(resulttable,2)
  barplot(ptable, legend=TRUE,beside=FALSE)
}



