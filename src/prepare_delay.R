# prepare_delay.R

# TODO - add header information

library(rsyncrosim)
library(tidyr)
library(data.table)
library(dplyr)
library(rightTruncation)

# Set the epi package project-scoped strings for jurisdiction and variables
jurisdictionBC <- "Canada - British Columbia"

env <- ssimEnvironment()
myScenario <- scenario()

# Add the required variables and jurisdictions to the SyncroSim project
saveDatasheet(myScenario, data.frame(Name = jurisdictionBC), "epi_Jurisdiction")   # Will ignore if exists already

# Get the inputs
inputSheet <- datasheet(myScenario, "modelCovidseir_DelayInputs")

# will write the delay data to this sheet
outputSheet <- datasheet(myScenario, "modelCovidseir_DelayOutputs", empty = T)
outputSheet[1,] <- NA
outputSheet <- transform(outputSheet, DelayDownloadDateTime = as.character(DelayDownloadDateTime))

# read the rda file from the Andersen Github
rdaFileName <- paste(env$TempDirectory, "delay_data.rda", sep='/')
download.file(
    gsub(' ', '', inputSheet$delayURL, fixed = TRUE),
    destfile=rdaFileName,
    quiet = TRUE
)
# mutate the data
delayData <- get(load(rdaFileName)) %>%
    data.table() %>%
    mutate_at(vars(time_to_report), as.integer)
# write  the delay data to file and log the name/time
delayDataFilename <- paste(env$TransferDirectory, "individual_reporting_delay_data.csv", sep='/')
write.csv(delayData, delayDataFilename)

outputSheet$DelayDataFile <-delayDataFilename
# outputSheet$delayURL <- inputSheet$delayURL
outputSheet$DelayDownloadDateTime <- as.character(Sys.time())

# save to output sheet
delayDatasheet <- datasheet(myScenario, "modelCovidseir_DelayData", empty=T)
delayDatasheet[nrow(delayData), ] <- NA
delayDatasheet$dateReported <- delayData$reported_date
delayDatasheet$dateSymptom <- delayData$symptom_onset_date
delayDatasheet$reportingGap <- delayData$time_to_report
saveDatasheet(myScenario, delayDatasheet, "modelCovidseir_DelayData")

# Generate Weibull distribution

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

# data sheets and file information
wParams <- datasheet(myScenario, "modelCovidseir_DelayWeibull", empty=T)
wParams[1,] <- NA
wParams$mleShape <- mleShape
wParams$mleScale <- mleScale
saveDatasheet(myScenario, wParams, "modelCovidseir_DelayWeibull")

weibullFilename <- paste(env$TempDirectory, "reporting_delay_weibull_parameters.csv", sep="/")

write.csv(wParams, weibullFilename, row.names=FALSE)

outputSheet$WeibullDataFile <- weibullFilename
# outputSheet$WeibullDownloadDateTime <- as.character(Sys.time())

saveDatasheet(myScenario, outputSheet, "modelCovidseir_DelayOutputs")

