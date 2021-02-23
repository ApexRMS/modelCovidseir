# getdelayData.R
# Download the symptom-onset-to-case-reporting delay data from the "rightTruncation" github
# and gather a list of user-provided break points for the model fitting

# Inputs:
#	Download URL
#	Table of changes to social distancing f Parameter
#		<date>	<beta mean>	<beta sd>
# Outputs:
#	delayData.rda - the (large) raw rda file retrieved from the rightTruncation github
#	delayDatasheet - a table of case reporting data, each row representing a unique case
#		<symptom onset date>	<case reporting date>	<gap (measured in days)>
#	DelayData.csv - CSV file of the case delay data, filtered by date for fitting the Model
#	SocialDistancngSegments.csv - CSV of the social distancing breakpoints  for fitting

library(rsyncrosim)
library(tidyr)
library(data.table)
library(dplyr)

currentScenario <- scenario()
env <- ssimEnvironment()

grabInput <- datasheet(currentScenario, "covidSEIR_ReportingDelays")

rdaFileName <- paste(env$TransferDirectory, "delayData.rda", sep='/')
download.file(
    gsub(" ", "", grabInput$delayURL, fixed = TRUE),
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

delayDataFilename <- paste(env$TransferDirectory, "DelayData.csv", sep="/")
write.csv(delayData, delayDataFilename)

delayDatasheet <- datasheet(currentScenario, "covidSEIR_DelayData")
delayDatasheet[nrow(delayData), ] <- NA
delayDatasheet$dateReported <- delayData$reported_date
delayDatasheet$dateSymptom <- delayData$symptom_onset_date
delayDatasheet$reportingGap <- delayData$time_to_report
saveDatasheet(currentScenario, delayDatasheet, "covidSEIR_DelayData")

runControl <- datasheet(currentScenario, "epi_RunControl")
runControl[1,] <- NA
runControl$MinimumTimestep = min(delayData$symptom_onset_date)
runControl$MaximumTimestep = max(delayData$reported_date)
runControl$MinimumIteration = 0
runControl$	MaximumIteration = 0
runControl$ModelHistoricalDeaths = TRUE
saveDatasheet(currentScenario, runControl, "epi_RunControl")

fSegments <- datasheet(currentScenario, "covidSEIR_FSegments")
segmentsFilename <- paste(env$TransferDirectory, "SocialDistancingSegments.csv", sep="/")
write.csv(fSegments, segmentsFilename)
