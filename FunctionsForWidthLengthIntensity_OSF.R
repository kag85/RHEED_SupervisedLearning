#Function for choosing peaks
library(tidyverse)
library(zoo)
library(imager)
library(EBImage)
library(pracma)
library(splines)

findPeaks <- function(path, i){
  EBimg <- readImage(paste0(path,i))
  nmaskt = thresh(EBimg, w = 15, h = 15, offset = 0.1)
  nmaskf = fillHull(opening(nmaskt, makeBrush(5, shape='disc')))
  img003 <- as.array(combine(getFrame(nmaskt, 3), getFrame(nmaskf, 3)), all = TRUE)
  #display(img003)
  imgBWbefore <- as.data.frame(bwlabel(img003))
  #display(normalize(imgBWbefore))
  IDb <- rownames(imgBWbefore)
  peaksDF1b <- cbind(IDb, imgBWbefore)
  
  #make bwlabels table long
  temp1b <- peaksDF1b %>%
    pivot_longer(-IDb, names_to = "VN", values_to = "value")
  
  #what <- temp1b %>%
  #  filter(value == 5)
  temp1b$VN <- gsub("V", "", temp1b$VN)
  temp1b <- temp1b %>%
    filter(VN < 600)
  
  tempNum <- temp1b %>%
    group_by(value) %>%
    tally() %>%
    filter(n > 50)
  
  tempSig <- inner_join(tempNum, temp1b)
  tempSig$IDb <- as.character(tempSig$IDb)
  tempSig$IDb <- as.numeric(tempSig$IDb)
  tempSig$VN <- as.numeric(tempSig$VN)
  
  tempCmaxb <- tempSig %>%
    filter(value > 0) %>%
    group_by(value) %>%
    summarize(Cmax = max(as.numeric(VN))) %>%
    arrange(Cmax) %>%
    head(1)
  # Pull out Cmax at the line drawn to crop
  #remove value corresponding to Cmax from the table
  
  #arrange the remaining spots and take the boundaries of their widths
  temp2b <- tempSig %>%
    filter(value > 0) %>%
    group_by(value) %>%
    summarize(Rmax = max(IDb), Rmin = min(IDb)) %>%
    anti_join(tempCmaxb)%>%
    mutate(avg = (Rmax +Rmin)/2)
  
  #table(temp1b$value)
  
  tempDFb <- data.frame(c(temp2b$Rmax, temp2b$Rmin)) %>%
    arrange(c.temp2b.Rmax..temp2b.Rmin.)
  
  spotBoundsb <- as.data.frame(rollmean(tempDFb$c.temp2b.Rmax..temp2b.Rmin., 2))
  if (dim(spotBoundsb)[1] > 1){
    #finds even rows (for the boundaries between spots):
    spotBounds1b <- seq(2, dim(spotBoundsb)[1], by = 2)
    #reports the boundaries:
    spotBounds2b <- as.data.frame(spotBoundsb[c(spotBounds1b), 1])
    
    spotBounds2b <- spotBounds2b %>%
      mutate(xVals = `spotBoundsb[c(spotBounds1b), 1]`) %>%
      select(xVals)
  }else {
    spotBounds2b <- data.frame(matrix(NA, nrow = 1, ncol = 1)) %>%
      rename(xVals = 'matrix.NA..nrow...1..ncol...1.')
  }
  tempCmaxb <- tempCmaxb %>%
    mutate(yvals = Cmax*(1))
  both <- list(spotBounds2b, tempCmaxb, temp2b)
  return(both)
}
##############################################################################################################

path <- "s:/ecg-members/Gliebe/Data/SRO/Examples_Frames/"

#A function that reads in all frames, selects a range of x values and then
#smooths the intensity curve across x to get values such as the locatoin
#of the maximum and the width of the peak
getwidthx <- function(temp) {
  temp <- temp %>%
    filter(x >= 200, x <= 400)
  model <- smooth.spline(x = temp$x, y = temp$value, df = 50)
  result <- pracma::findpeaks(model$y, minpeakheight = .2) 
  dfresult <- as.data.frame(result)
  
  if (length(dfresult != 0)) {
    widthsx <- dfresult %>%
      mutate(width = V4 - V3) %>%
      mutate(peakdist = (V2 - lag(V2, k = 1))) #%>%
    #select(widthsx, V2)
  } else {
    widthsx = NA
  }
  return(widthsx)
}

#The for loop below runs the function for an entire folder of frames (that comprise
#one video) and give each frame a number so that I can represent time across the
#videos.  The frames are then binded together.

files <- list.files(path)
df <- data.frame()

frame <- 1
for (i in files){
#i <- "thumb_0001.jpg"
  print(i)
  both <- findPeaks(path, i)
  spotBounds2b <- both[[1]]
  tempCmaxb <- both[[2]]
  avg <- both[[3]]
  temp1 <- load.image(paste0(path,i))
  temp1 <- as.data.frame(temp1)
  
  xtest <- temp1 %>%
    filter(y >= tempCmaxb$yvals) %>%
    group_by(y) %>%
    do(result = getwidthx(.))
  
  xtest <- unnest(data = xtest, cols = c(result)) 
  
  xtest$frame <- frame
  xtest$file <- i
  if (length(names(xtest)) != 10)
  {xtest$result <- NA}
  
  xtest <- xtest[complete.cases(xtest[,2]),]
  shift <- min(avg$avg) - min(xtest$V2)
  dftry <- xtest %>%
    mutate(try = V2 + shift)
  dftry3 <- data.frame()
  if (is.na((spotBounds2b)[1,1]) != TRUE){ 
    if (dim(spotBounds2b)[1] > 1) {
      tempDFb <- mean(spotBounds2b$xVals)
      
      orderingPOS <- spotBounds2b %>%
        mutate(diff = tempDFb - xVals) %>%
        filter(diff > 0) 
      
      odds <- seq(1, dim(orderingPOS)[1]*2, by = 2)
      orderingPOS <- orderingPOS %>%
        arrange(diff) %>%
        mutate(ID = odds)
      orderingNEG <- spotBounds2b %>%
        mutate(diff = tempDFb - xVals) %>%
        filter(diff < 0)
      evens <- seq(2, dim(orderingNEG)[1]*2, by = 2)
      orderingNEG <- orderingNEG %>%
        arrange(desc(diff)) %>%
        mutate(ID = evens)
      full <- rbind(orderingNEG, orderingPOS) %>%
        arrange(ID)
      
      
      dfnew <- dftry
      for (i in 1:dim(full)[1]){
        
        dfnew[, ncol(dfnew) + 1] <- full$xVals[i] - dfnew$try
        names(dfnew)[ncol(dfnew)] <- paste0("dist_", i)
      }
      
      test <- dfnew[, 12:(11 + dim(full)[1])]
      test$peakNum <- 0
      dim(test)[1]
      
      for (i in 1:dim(test)[1]) {
        if ((sum(test[i,] < 0)) > 0){
          x <- t(test[i,])
          test$peakNum[i] <- which(grepl(max(x[x < 0]), x))
        } else {
          test$peakNum[i] <- (dim(full)[1] + 1)
        }
      }
      dftry2 <- cbind(dftry, test) 
      dftry2 <- dftry2 %>%
        select(y, V1, V2, V3, V4, width, peakdist, frame, file, try, peakNum)
      frame <- frame + 1
      dftry2$frame <- frame
      dftry2$file <- i
      dftry3 <- rbind(dftry3, dftry2) 
    }else{
      frame = frame + 1
      dftry3 <- dftry %>%
        mutate(peakNum = 1)
      dftry3 <- dftry3 %>%
        select(y, V1, V2, V3, V4, width, peakdist, frame, file, try, peakNum)
      #for (i in 1:dim(dftry)[1]){
      #  if (dftry$try[i] > spotBounds2b$xVals
    }
  }else{
    dftry3 <- dftry %>%
      mutate(peakNum = 1)
    frame <- frame + 1
    dftry3 <- dftry3 %>%
      select(y, V1, V2, V3, V4, width, peakdist, frame, file, try, peakNum)
  }
  
  
  df <- rbind(df, dftry3)
  
}

idf <- df #%>%
  #dplyr::select(-result)

dfnum <- idf
dfnum <- dfnum[complete.cases(dfnum[,3]),]


dfstat <- df %>%
  #select(y, width, peaknum, file) %>%
  group_by(file, peakNum) %>%
  summarise(mean(width))
#The above finds the mean of the widths along a peak.  I also could have used 
#The maximum width across a peak

dfstat2 <- df %>%
  group_by(file, peakNum) %>%
  summarise(max(width))

joinddf <- left_join(df, dfstat)
joinddf <- left_join(joinddf, dfstat2)

joinddf <- joinddf[complete.cases(joinddf[,10:12]),]
write.csv(joinddf, "s:/ecg-members/Gliebe/Data/SRO/Example1OSF.csv")
#this saves a file of all the intensity distribution data across x holding y constant.

#The below function and loop is similar to the one above, except now it
#is working along the y axis to obtain lengths of peaks
getwidthy <- function(temp) {
  modely <- smooth.spline(x=temp$y, y=temp$value, df= 60)
  result <- try((pracma::findpeaks(modely$y, minpeakheight = .1)), silent = TRUE) 
  dfresulty <- as.data.frame(result)
  
  if (length(dfresulty != 0)) {
    widthsy <- dfresulty %>%
      mutate(widthsy = V4 - V3)
  } else {
    widthsy = NA
  }
}

frame <- 0
files <- list.files(path)
dfy <- data.frame()
for (i in files){
#i <- "thumb0001.jpg"
  print(i)
  both <- findPeaks(path, i)
  spotBounds2b <- both[[1]]
  tempCmaxb <- both[[2]]
  avg <- both[[3]]
  temp1 <- load.image(paste0(path,i))
  temp1 <- as.data.frame(temp1)
  
  ytest <- temp1 %>%
    filter(y >= tempCmaxb$Cmax) %>%
    group_by(x) %>%
    do(result = getwidthy(.)) 
  ytest <-unnest(data = ytest, cols = c(result))
  ytest <- ytest[complete.cases(ytest[,2]),]
  ytest$frame <- frame
  ytest$file <- i
  
  dfytry3 <- data.frame()
  
  
  if (is.na((spotBounds2b)[1,1]) != TRUE){ 
    if (dim(spotBounds2b)[1] > 1) {
      tempDFb <- mean(spotBounds2b$xVals)
      
      orderingPOS <- spotBounds2b %>%
        mutate(diff = tempDFb - xVals) %>%
        filter(diff > 0) 
      
      odds <- seq(1, dim(orderingPOS)[1]*2, by = 2)
      orderingPOS <- orderingPOS %>%
        arrange(diff) %>%
        mutate(ID = odds)
      orderingNEG <- spotBounds2b %>%
        mutate(diff = tempDFb - xVals) %>%
        filter(diff < 0)
      evens <- seq(2, dim(orderingNEG)[1]*2, by = 2)
      orderingNEG <- orderingNEG %>%
        arrange(desc(diff)) %>%
        mutate(ID = evens)
      full <- rbind(orderingNEG, orderingPOS) %>%
        arrange(ID)
      
      
      dfnew <- ytest
      for (i in 1:dim(full)[1]){
        
        dfnew[, ncol(dfnew) + 1] <- full$xVals[i] - dfnew$x
        names(dfnew)[ncol(dfnew)] <- paste0("dist_", i)
      }
      
      test <- dfnew[, 9:(8 + dim(full)[1])]
      test$peakNum <- 0
      
      for (i in 1:dim(test)[1]) {
        if ((sum(test[i,] < 0)) > 0){
          x <- t(test[i,])
          test$peakNum[i] <- which(grepl(max(x[x < 0]), x))
        } else {
          test$peakNum[i] <- (dim(full[1]) + 1)
        }
      }
      dftry2 <- cbind(ytest, test) 
      dftry2 <- dftry2 %>%
        select(x, V1, V2, V3, V4, widthsy, frame, file, peakNum)
      frame <- frame + 1
      dftry2$frame <- frame
      dftry2$file <- i
      dfytry3 <- rbind(dfytry3, dftry2)
      dfy <- rbind(dfy, dfytry3)
    }else{
      frame = frame + 1
      dfytry3 <- ytest %>%
        mutate(peakNum = 1)
      dfy <- rbind(dfy, dfytry3)
    }
  }else{
    dfytry3 <- ytest %>%
      mutate(peakNum = 1)
    dfy <- rbind(dfy, dfytry3)
    frame <- frame + 1
  }
  
  
  dfy <- rbind(dfy, dfytry3)
  
}

dfynum <- dfy %>%
  group_by(frame, peakNum) %>%
  mutate(maxL = max(widthsy))
#dfynum = dfynum[!duplicated(dfynum[,6:7]),]

write.csv(dfynum, "s:/ecg-members/Gliebe/Data/SRO/Example_Y_OSF.csv")
#this just saves a file of all of the information about the intensity distribution across
#y holding x constant.




#The below condenses the previous data frames into a single one containing the
#maximum width, length, and intensities for each diffraction spot of each
#frame.  The final file is combine2 and you can save it to your file path.
#There are other data frames that are saved along the way to ensure no progress
# is lost; however, these are not necessary and can be commented out.

reduce <- read.csv("s:/ecg-members/Gliebe/Data/SRO/Example1OSF.csv")
#library(tidyverse)
take2 <- reduce %>%
  dplyr::group_by(frame, peakNum) %>%
  dplyr::summarize(max(width)) %>%
  dplyr::rename(width = 'max(width)')
take3 <- merge(take2, reduce, by = c("frame", "peakNum", "width"))
take3 <- take3 %>% arrange(frame, peakNum)  
take3 <- take3 %>%
  dplyr::select(width, frame, peakNum)
take4 <- unique(take3)
reduce <- reduce %>%
  dplyr::select(-X) %>%
  rename( xIntensity = V1, xMaxLoc = V2, xPeakStart = V3, xPeakEnd = V4, Xwidth = width)

reduced <- reduce %>%
  arrange(frame, peakNum) %>%
  mutate(time = (4/5)*frame)
maxXInt <- reduced %>%
  group_by(frame, peakNum) %>%
  summarise(max(xIntensity)) %>%
  rename(xIntensity = 'max(xIntensity)')
reducedI <- inner_join(reduced, maxXInt)
trial <- merge(reducedI, take3, by = c("frame", "peakNum"))
reducedI <- reducedI %>% arrange(frame, peakNum)
finalI <- inner_join(reducedI, take4)

reducey <- read.csv("s:/ecg-members/Gliebe/Data/SRO/Example_Y_OSF.csv")

reducey <- reducey %>%
  select(-X) %>%
  rename(yIntensity = V1, yMaxLoc = V2, yPeakStart = V3, yPeakEnd = V4)
reducedy <- reducey[reducey$widthsy==reducey$maxL,]
reducedy <- reducedy %>%
  arrange(frame, peakNum) %>%
  mutate(time = (4/5)*frame)
less <- reducedy %>%
  select(-yMaxLoc, -yPeakStart, -yPeakEnd, -maxL, -file, -x, -yIntensity)
less <- unique(less)
maxYInt <- reducedy %>%
  group_by(frame, peakNum) %>%
  summarise(max(yIntensity)) %>%
  rename(yIntensity = 'max(yIntensity)')
reducedyI <- inner_join(less, maxYInt)
reducedyI <- reducedyI %>% arrange(frame, peakNum)

combine2 <- inner_join(reducedI, reducedyI)
combine2 <- combine2 %>%
  arrange(frame, peakNum) %>%
  select(-y) 

combine2 <- combine2 %>%
  mutate(sample = "Example")# %>%
  #mutate(mat = "LNTO")


write.csv(combine2, "s:/ecg-members/Gliebe/Data/SRO/Example_Result.csv")


