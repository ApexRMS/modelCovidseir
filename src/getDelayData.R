library(rsyncrosim)
library(tidyr)
library(data.table)
library(dplyr)
library(rightTruncation)

myScenario <- scenario()
env <- ssimEnvironment()

inputSheet <- datasheet(myScenario, "modelCovidseir_ReportingDelays")
downloadURL <- inputSheet$delayURL

rdaFileName <- paste(env$TransferDirectory, "delayData.rda", sep='/')
download.file(
    gsub(' ', '', downloadURL, fixed = TRUE),
    destfile=rdaFileName,
    quiet = TRUE
)

delayData <- get(load(rdaFileName)) %>%
    data.table() %>%
    mutate_at(vars(time_to_report), as.integer) %>%
    subset(
        "2020-08-01" < symptom_onset_date &
            reported_date < "2020-12-01"
    )

delayDataFilename <- paste(env$TransferDirectory, "DelayData.csv", sep='/')
write.csv(delayData, delayDataFilename)

outputSheet <- datasheet(myScenario, "modelCovidseir_DelayOutput", empty = T)
outputSheet <- transform(outputSheet, DownloadDateTime = as.character(DownloadDateTime))

outputSheet = addRow(outputSheet, c(delayDataFilename, as.character(Sys.time())))
saveDatasheet(myScenario, outputSheet, "modelCovidseir_DelayOutput")

delayDatasheet <- datasheet(myScenario, "modelCovidseir_DelayData")
delayDatasheet[nrow(delayData), ] <- NA
delayDatasheet$dateReported <- delayData$reported_date
delayDatasheet$dateSymptom <- delayData$symptom_onset_date
delayDatasheet$reportingGap <- delayData$time_to_report
saveDatasheet(myScenario, delayDatasheet, "modelCovidseir_DelayData")

# set the first day so that we can calculate the time that passed in days
hDay0 <- min(delayData$symptom_onset_date)
# set the number of rows and columns of the matrix
numRows <- (max(delayData$symptom_onset_date)-hDay0)[[1]]
numCols <- (max(delayData$reported_date)-hDay0)[[1]]

# fill the Hnr matrix with the case data as described in the previous long comment
Hnr <- matrix(0, nrow=numRows, ncol=numCols)
for(i in 1:nrow(delayData))
{
    onset <- (delayData[i]$symptom_onset_date - hDay0)[[1]]
    repor <- (delayData[i]$reported_date - hDay0)[[1]];
    Hnr[onset, repor] <- Hnr[onset, repor]+1
}

mleRes = suppressWarnings(nlm(
    f = negLL_Weibull_counts_matrix,
    p = c(2, 7),
    h_nr = Hnr
))

# the parameters for the distribution
mleShape <- mleRes$estimate[1]
mleScale <- mleRes$estimate[2]

wParams <- datasheet(myScenario, "modelCovidseir_WeibullParameters")
wParams[1,] <-NA
wParams$mleShape <- mleShape
wParams$mleScale <- mleScale
saveDatasheet(myScenario, wParams, "modelCovidseir_WeibullParameters")

write.csv(wParams, paste(env$TransferDirectory, "modelCovidseir_WeibullParameters.csv", sep="/"), row.names=FALSE)

# smooth Weibull curve to be superimposed
weibullX <- seq(min(delayData$time_to_report), max(delayData$time_to_report), by=0.05)
theDist <- data.table(x = weibullX, y = nrow(delayData)*dweibull(weibullX, shape=mleShape, scale=mleScale))

weibullDist <- datasheet(myScenario, "modelCovidseir_DelayDist")
weibullDist[length(weibullX),] <-NA
weibullDist$x <- theDist$x
weibullDist$y <- theDist$y
saveDatasheet(myScenario, weibullDist, "modelCovidseir_DelayDist")
