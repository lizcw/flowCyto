#R script for parsing bit fields
mydir <- 'D:\\Projects\\FlowCytoData\\data'
datafile <-'16092016_rainbow med_002.fcs'
testdata1 <- file.path(mydir,datafile)
x1 = read.FCS(testdata1,alter.names=TRUE, transformation = "linearize-with-PnG-scaling")
summary(x1)
df <- data.frame(exprs(x1)[,1:10])#or subset for testing[1:10,1:10])
a <- df$Classifier.Bits #ie first bit for PHASE OK
r <- df$ROI.Bits.1.16
a[1]
n <-4
# for (i in 1:10){
#   dropphase <- intToBits(a[i]) 
#   print(class(dropphase))
#   dp <-packBits(dropphase, "integer")
#   
#   print(dropphase)
#   print(dp)
# }

# Comparison of Bin from Classfier.Bits vs ROI.Bits.1.16 - LATTER IS BEST
a1 <- a -32768
for (i in 1:dim(df)[1]){
  pwr = a[i]-32768
  bin = log2(a1[i])+1
  if (r[i] > 1){
    bindiff = as.integer(log2(r[i]-1))
    if (bin-bindiff != 0){
      print(paste("a1: ",a[i], "=", bin, " roi=", bindiff))
    }
    
  }else{
    bindiff = 'NA'
  }
  #print(paste("a: ",a[i], "=", log2(pwr)+1))
  
  
}

