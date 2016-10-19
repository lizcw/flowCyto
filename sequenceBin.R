#R script for parsing bit fields
mydir <- 'D:\\Projects\\FlowCytoData\\data'
datafile <-'16092016_rainbow med_002.fcs' #10102016_Rainbow6a_003.fcs'
testdata1 <- file.path(mydir,datafile)
#inbuilt data
testdata1 <- system.file("extdata","0877408774.B08", package="flowCore")

x1 = read.FCS(testdata1,alter.names=TRUE, transformation = FALSE)
summary(x1)


df <- data.frame(exprs(x1)[,1:10])#or subset for testing[1:10,1:10])
r <- df$ROI.Bits.1.16
r2 <- df$ROI.Bits.17.32

df$SortBucket <-mapply(sortbucket,df$ROI.Bits.1.16)
deltaBucket <- diff(df$SortBucket)
df$DeltaBucket <-c(0,deltaBucket)
#Show all data
dfres1 <-hist(df$SortBucket,breaks=6, freq=TRUE,main=paste(da), col='yellow', plot=TRUE,xlab="Bucket diff",xlim=range(0,14))
text(dfres1$mids, dfres1$density, dfres1$counts, adj = c(.5, -.5), col = "blue")
dfres <-hist(abs(df$DeltaBucket),breaks=6, freq=TRUE,main=paste("Sequential Bucket swings"), col='yellow', plot=TRUE,xlab="Bucket diff",xlim=range(0,14))
text(dfres$mids, dfres$density, dfres$counts, adj = c(.5, -.5), col = "blue")

binall <- list()    # every step compared
binsorted <- list() # exclude unsorted (includes sort aborts)
binactual <- list() # exclude all aborts
# Check of Bin Sequence from ROI.Bits.1.16
# use each step vs only selected bins - produce histogram
ctr <-1
for (i in 2:dim(df)[1]){
  
  if (r[i] > 1 && r[i-1] > 1){ #exclude unsorted for both
    bin = as.integer(log2(r[i]-1))
    binprev = as.integer(log2(r[i-1]-1))
    bindiff = binprev-bin
    binsorted[ctr]<-binprev-bin
    if (abs(bindiff) > 6){
      print(paste("i: ",i, "=", bindiff))
    }
    
    ctr <- ctr + 1
    
  }else{
    #print("NA")
  }
  
}

dfres <-hist(unlist(binsorted),breaks=6, freq=TRUE,main=paste("Sequential Bucket swings - exclude unsorted"), col='yellow', plot=TRUE,xlab="Bucket diff",xlim=range(-14,14))
text(dfres$mids, dfres$density, dfres$counts, adj = c(.5, -.5), col = "blue")