library("MALDIquant")
library("MALDIquantForeign")
#example <- system.file("exampledata", package="MALDIquantForeign")
#spec <- import(file.path(example,"brukerflex"), verbose = FALSE)
#spec <- transformIntensity(spec, method="sqrt")
#spec <- smoothIntensity(spec, method="SavitzkyGolay", halfWindowSize=10)
#base <- estimateBaseline(spec [[1]], method="SNIP", iterations=100)
#spec <- removeBaseline(spec, method="SNIP", iterations=100)
#spec <- calibrateIntensity(spec, method="TIC")
#spec <- alignSpectra(spec,SNR=2, tolerance=0.002, warpingMethod="lowess")

#appp <- sapply(spec, function(y)metaData(y)$sampleName)
#sample <- factor(unlist(appp))

exampleDirectory <- system.file("/Users/audreycui01/Documents/Data Preprocessing/Ovarian Cancer Dataset/Ovarian Cancer")

#filePath <- file.path("/Users/audreycui01/Documents/Data Preprocessing/Ovarian Cancer Dataset/Ovarian Cancer", "Ovarian Cancer daf-0614.csv")
cancerDir <- file.path("/Users/audreycui01/Documents/Data Preprocessing/Ovarian Cancer Dataset/OvarianCD_PostQAQC/Cancer")
controlDir <- file.path("/Users/audreycui01/Documents/Data Preprocessing/Ovarian Cancer Dataset/OvarianCD_PostQAQC/Normal")
temp = list.files(path = controlDir)
temp2 = list.files(path = cancerDir)
List = list()
counter = 1
max = 0
massList = list()
massList[[1]]<-0
for (j in 1:(length(temp)+length(temp2)))
{
  print("HERE")
  if (j<=length(temp))
    file <- file.path(controlDir, temp[j])
  else 
    file <- file.path(cancerDir, temp2[j-length(temp)])
  spectra <- importTxt(file, verbose=FALSE)
  spectra <- transformIntensity(spectra, method="sqrt")
  spectra <- smoothIntensity(spectra, method="SavitzkyGolay", halfWindowSize=10)
  baseline <- estimateBaseline(spectra [[1]], method="SNIP", iterations=100)
  spectra <- removeBaseline(spectra, method="SNIP", iterations=100)
  spectra <- calibrateIntensity(spectra, method="TIC")
  spectra <- alignSpectra(spectra, halfWindowSize=20,SNR=5, tolerance=0.05, warpingMethod="lowess")
  app <- sapply(spectra, function(x)metaData(x))
  samples <- factor(unlist(app))
  
  avgSpectra <- averageMassSpectra(spectra, labels=samples, method="mean")
  
  noise <- estimateNoise(avgSpectra[[1]])
  peaks <- detectPeaks(avgSpectra, method="MAD", halfWindowSize=20, SNR=5)
  peaks <- binPeaks(peaks, tolerance=0.005)
  peaks <- filterPeaks(peaks, minFrequency=0.2)
  
  featureMatrix <- intensityMatrix(peaks, avgSpectra)
  
  
  for (i in 1:ncol(featureMatrix))
  {
    theMass = as.double(dimnames(featureMatrix)[[2]][i])
    
    done <- FALSE
    if (lengths(massList)==1 && massList[[1]] == 0)
    {
      massList[[1]] <- theMass
      done <- TRUE 
    }
    
    for (k in 1:length(massList))
    {
      if (abs(theMass-massList[[k]]) < 0.001 || done == TRUE)
      {
        break
      }
      else if (theMass < massList[[k]])
      {
        massList <- append(massList, theMass, k-1)
        break
      }
      
      else if (k==length(massList))
      {
        massList[[k+1]]<- theMass
        break
      }
      
    }
    
  }
  
}

miMatrix <- matrix(ncol = 1+length(massList), nrow = length(temp)+length(temp2))
colnames(miMatrix) = c(massList, "Outcome") 
for (j in 1:(length(temp)+length(temp2)))
{
  if (j<=length(temp))
  {
    file <- file.path(controlDir, temp[j])
    miMatrix[j, ncol(miMatrix)]= as.integer(0)
    
  }
  
  else 
  {
    file <- file.path(cancerDir, temp2[j-length(temp)])
    miMatrix[j, ncol(miMatrix)]= as.integer(1)
  }
  
  spectra <- importTxt(file, verbose=FALSE)
  spectra <- transformIntensity(spectra, method="sqrt")
  spectra <- smoothIntensity(spectra, method="SavitzkyGolay", halfWindowSize=10)
  baseline <- estimateBaseline(spectra [[1]], method="SNIP", iterations=100)
  spectra <- removeBaseline(spectra, method="SNIP", iterations=100)
  spectra <- calibrateIntensity(spectra, method="TIC")
  spectra <- alignSpectra(spectra, halfWindowSize=20,SNR=5, tolerance=0.05, warpingMethod="lowess")
  app <- sapply(spectra, function(x)metaData(x))
  samples <- factor(unlist(app))
  
  avgSpectra <- averageMassSpectra(spectra, labels=samples, method="mean")
  
  noise <- estimateNoise(avgSpectra[[1]])
  peaks <- detectPeaks(avgSpectra, method="MAD", halfWindowSize=20, SNR=5)
  peaks <- binPeaks(peaks, tolerance=0.005)
  peaks <- filterPeaks(peaks, minFrequency=0.2)
  
  featureMatrix <- intensityMatrix(peaks, avgSpectra)
  index <- 1
  for (i in 1:length(featureMatrix))
  {
    theMass = as.double(dimnames(featureMatrix)[[2]][i])
    for (k in index:length(massList))
    {
      if(abs(theMass - massList[[k]]) < 0.001)
      {
        miMatrix[j, k] <- featureMatrix[1,i]
        index <- k+1
        break
      }
      else 
        miMatrix[j, k] <- 0
    }
    
  }
  if (index<length(massList))
  {
    for(i in index:length(massList))
      miMatrix[j, i] <- 0
  }
  
}
print(miMatrix)
write.csv(miMatrix, file = "allMI3.csv")