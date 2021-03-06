########Envt
library("plyr")
library("ggplot2")
library("flowCore")
################################################################
# Generate histograms of sorted buckets
################################################################
#Read data file
mydir <- 'D:\\lizcw\\Projects\\RProjects'
mydir <- 'D:\\Projects\\FlowCytoData\\data'
datafilelist <-c('16092016_rainbow med_002.fcs',
                 '10102016_Rainbow6a_003.fcs',
                 '10102016_Rainbow6a_004.fcs',
                 '10102016_Rainbow6a_005.fcs',
                 '25102016_rainbow1.fcs',
                 '25102016_rainbow1_001.fcs')



############################# Functions ########################

CLOCKSTEP <- 0.017624999999999998
DISTANCESTEP <- 0.28199999999999997
DISTANCESPAN <- 1155.0719999999999
ABORTLIMIT <- 12 #abort when in drophase to 16

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
  prevdist <- bitwAnd(bitwShiftL(dropPhaseBits, 4), 0xffff)
  return(log2(prevdist))
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
sortstatus <- function(classifier,sortresult, drop){
  if (is.nan(sortresult) && is.nan(classifier)){
    sortstatus <- 'unsorted'
    
  }else if (classifier > 1 && is.nan(sortresult)){
    # NOT ENOUGH INFO TO DECIDE
    #if (drop < ABORTLIMIT){
    #   sortstatus <- 'elect abort'
    # }else{
    #   sortstatus <- 'sort abort'
    # }
    sortstatus <- 'abort'
  }else{
    sortstatus <- 'sorted'
  }
  return(sortstatus)
}

#Calculate bucket swing for whole column - either classifier or sortresult
#- between two sequential and sorted events
# - all - save with NANs to add back to dataframe, else just list for histogram
# Within time limit only - timeinterval is corresponding deltatime
diffbucket <- function(blist,all, timeinterval, timemax){
  if (class(blist) =='list'){
    blist <- unlist(blist)
  }
  diffbucket <- list()
  ctr <- 0
  for (i in 2:length(blist)){
    a<-blist[i]
    b<-blist[i-1]
    if (!is.nan(a) && !is.nan(b)){
      bindiff <- a-b
      #print(paste("a=",a," b=",b,"a-b=",a-b))
      ctr <- ctr + 1
    }else{
      bindiff <- NaN
    }
    if (all){
      if (timeinterval[i] <= timemax){
        diffbucket[i] <- bindiff
      }else{
        diffbucket[i] <- NaN
      }
      
    }else if(!is.nan(bindiff) && timeinterval[i] <= timemax ){
      diffbucket[ctr] <- bindiff
    }
    
  }
  print(paste("Diffbucket=",length(diffbucket)))
  return(diffbucket)
}

#Plot multiple plots on one page
# from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)
multiplot <- function(p1,p2,p3,p4=NULL, cols=1, layout=NULL, pagetitle="Histograms") {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(p1,p2,p3,p4)
  
  numPlots = length(plots)
 
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(3,2)))
    grid.text(pagetitle, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
    print(p1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
    if (!is.null(p4)){
      print(p3, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
      print(p4, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
    }else{
      print(p3, vp = viewport(layout.pos.row = 3, layout.pos.col = 1:2))
    }
    
   
  }
}

#With 2 IN SEQUENCE
filter2Sequential <- function(d){
  keep <- c() #TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)
  c <- 1
  for(i in 1:(nrow(d)-2)){
    if (i == c){
      row1 <- rownames(d)[c]
      row2 <- rownames(d)[c+1]
      if (as.numeric(row2) - as.numeric(row1) == 1 ){ 
        print(paste("i=",i," row1=", row1," row2=", row2," row3=", row3))
        keep <- c(keep, TRUE,TRUE,TRUE)
        c <- c+2
      }else{
        keep <- c(keep, FALSE)
        c <- c+1
      }
    }
  }
  return(d[keep,])
}


#With 3 IN SEQUENCE
filter3Sequential <- function(d){
  keep <- c() #TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)
  c <- 1
  for(i in 1:(nrow(d)-3)){
    if (i == c){
      row1 <- rownames(d)[c]
      row2 <- rownames(d)[c+1]
      row3 <- rownames(d)[c+2]
      if ((as.numeric(row2) - as.numeric(row1) == 1 ) && (as.numeric(row3) - as.numeric(row2) == 1 )){ 
        print(paste("i=",i," row1=", row1," row2=", row2," row3=", row3))
        keep <- c(keep, TRUE,TRUE,TRUE)
        
        c <- c+3
      }else{
        keep <- c(keep, FALSE)
        c <- c+1
      }
    }
  }
  return(d[keep,])
}

#Extract triplet bins in sequence per column data
extractTriplets <- function(dfcol){
  t <-c()
  i <- 1
  while(i < (length(dfcol)-3)){
    t <-rbind(t,dfcol[c(i,i+1,i+2)])
    i<-i+3
  }
  return(t)
}

aacount <- function(a){
  b<-length(a)
  return(b)
}

########################## MAIN ####################################


##################
#for (da in datafilelist){
da = datafilelist[5]
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
df$DeltaTime <-c(deltaTime,0) #post drop interval
#df$DeltaTime <-c(0,deltaTime) #previous drop interval
starttime <- df$TimeCalc[1]
df$ElapsedTime <-mapply(elapsed,df$TimeCalc,starttime)
timemax <- max(df$DeltaTime)
#Generate drop phase
df$DropphaseNum <- mapply(dropphase, df$Drop.Phase)
#Generate bucket numbers and sort status
df$ClassBucket <-mapply(bucketnumber,df$Classifier.Bits) #Classified to bin
df$SortBucket <-mapply(bucketnumber,df$Sort.Result.Bits) #exclude sort aborts
deltaBucketAll <- unlist(diffbucket(df$ClassBucket,TRUE, df$DeltaTime, timemax))
df$DeltaBucket <-c(0,deltaBucketAll)
deltaBucketSorted <- unlist(diffbucket(df$SortBucket,TRUE, df$DeltaTime, timemax))
df$DeltaBucketSorted <-c(0,deltaBucketSorted)
timemax25 <- quantile(df$DeltaTime)[[2]] #25% in this time interval
deltaBucketShort <- unlist(diffbucket(df$SortBucket,TRUE, df$DeltaTime, timemax25))
df$DeltaBucketShort <-c(0,deltaBucketShort)
df$SortStatus <-mapply(sortstatus,df$ClassBucket,df$SortBucket, df$DropphaseNum)
df$Prev <-mapply(prevdist, df$Drop.Phase)
shortdf <- subset(df,DeltaTime <= timemax25)


#############################################################
#Extract sequence data
df1 <- df[df$DeltaTime<=43,] #Within one drop limit
df1 <- df1[df1$SortStatus=="sorted",]
df2 <-filter3Sequential(df1)
deltatriplets<- extractTriplets(df2$DeltaBucket)
deltatriplets
#Write to file
fname <-paste0(da,"_sequences.csv")
write.table(deltatriplets, file.path(mydir,fname), sep=",")
#Pull out longer sequences
df3<-which(diff(as.numeric(rownames(df1)))==1)
a<-split(df3, cumsum(c(0, diff(df3) != 1))) ##INCORRECT - SKIPPING INFO
a1 <- mapply(aacount,a)
title <- paste("Sequence lengths for ", da)
m <- count(a1)

#Select sequences with 6 in a row
#t1 <- matrix(c(a[which(a1==1)],a[which(a1==2)],a[which(a1==3)],a[which(a1==4)],a[which(a1==5)],a[which(a1==6)]),ncol=6,byrow=FALSE)




#####pAUSE IN LOOP
#cat ("Press [enter] to continue")
#line <- readline()
#}

