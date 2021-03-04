library(rsyncrosim)
library(tidyr)
library(data.table)
library(dplyr)
library(rightTruncation)

myScenario <- scenario()
env <- ssimEnvironment()

inputSheet <- datasheet(myScenario, "modelCovidseir_DelayFileInfo")

outputSheet <- datasheet(myScenario, "modelCovidseir_DelayFileInfo", empty = T)
outputSheet[1,] <- NA
outputSheet <- transform(outputSheet, DelayDownloadDateTime = as.character(DelayDownloadDateTime), WeibullDownloadDateTime = as.character(WeibullDownloadDateTime))

rdaFileName <- paste(env$TempDirectory, "delayData.rda", sep='/')
download.file(
    gsub(' ', '', inputSheet$delayURL, fixed = TRUE),
    destfile=rdaFileName,
    quiet = TRUE
)

delayData <- get(load(rdaFileName)) %>%
    data.table() %>%
    mutate_at(vars(time_to_report), as.integer)

delayDataFilename <- paste(env$TransferDirectory, "CaseReportingDelayData.csv", sep='/')
write.csv(delayData, delayDataFilename)

outputSheet$DelayDataFile <-delayDataFilename
outputSheet$delayURL <- inputSheet$delayURL
outputSheet$DelayDownloadDateTime <- as.character(Sys.time())

delayDatasheet <- datasheet(myScenario, "modelCovidseir_RawDelayData", empty=T)
delayDatasheet[nrow(delayData), ] <- NA
delayDatasheet$dateReported <- delayData$reported_date
delayDatasheet$dateSymptom <- delayData$symptom_onset_date
delayDatasheet$reportingGap <- delayData$time_to_report
saveDatasheet(myScenario, delayDatasheet, "modelCovidseir_RawDelayData")

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

wParams <- datasheet(myScenario, "modelCovidseir_WeibullParameters", empty=T)
wParams[1,] <- NA
wParams$mleShape <- mleShape
wParams$mleScale <- mleScale
saveDatasheet(myScenario, wParams, "modelCovidseir_WeibullParameters")

weibullFilename <- paste(env$TempDirectory, "modelCovidseir_WeibullParameters.csv", sep="/")

write.csv(wParams, weibullFilename, row.names=FALSE)

outputSheet$WeibullDataFile <- weibullFilename
outputSheet$WeibullDownloadDateTime <- as.character(Sys.time())

saveDatasheet(myScenario, outputSheet, "modelCovidseir_DelayFileInfo")
