library(tidyr)
library(data.table)
library(dplyr)
library(rsyncrosim)

library(covidseir)

future::plan(future::multisession)

myScenario <- scenario()
env <- ssimEnvironment()

genParams <- datasheet(myScenario, "modelCovidseir_General")
projParams <- datasheet(myScenario, "modelCovidseir_ProjParams")
runControl <- datasheet(myScenario, "epi_RunControl")

caseData <- datasheet(myScenario, "epi_DataSummary") %>% transform(Timestep = as.Date(Timestep))

totalDuration <- (as.Date(runControl$MaximumTimestep) - min(caseData$Timestep))[[1]]
daysToProject <-(as.Date(runControl$MaximumTimestep) - as.Date(max(caseData$Timestep)))[[1]]

# change to
theFit <-readRDS(sprintf("%s\\%s.rds", env$TempDirectory, genParams$RDS))

if((daysToProject<=0) | (is.na(daysToProject))) stop()

projParams <- datasheet(myScenario, "modelCovidseir_ProjParams")

simData <- data.table(matrix(ncol=3, nrow=(totalDuration+1)*runControl$MaximumIteration))
names(simData) <-c("Iteration", "Date", "simCases")
simData[, Iteration:=as.integer(Iteration)]
simData[, Date:=as.Date(Date)]
simData[, simCases:=as.numeric(simCases)]

# in the case of a huge number of iterations,we'll have a table already set up
rowCounter <- 1
firstDay <- min(caseData$Timestep, na.rm=TRUE)
lastDay <- totalDuration
lut <- dplyr::tibble(
    day = seq_len(lastDay),
    date = seq(firstDay, firstDay + length(day) - 1, by = "1 day")
)

simData <- datasheet(myScenario, name = "epi_DataSummary", empty = T, optional = T, lookupsAsFactors = F)
simData <- transform(simData, Timestep = as.Date(Timestep))

for(index in 1:runControl$MaximumIteration)
{
    theProj <- covidseir::project_seir(
        theFit,
        iter = 1:1,
        forecast_days = daysToProject, # number of days into the future to predict
        f_fixed_start = max(theFit$days) + projParams$StartChange, # the day at which to start changing f
        f_multi = rep(1.1, daysToProject - projParams$StartChange + 1),
        f_multi_seg = projParams$FSegment, # which f segment to use
        parallel = (projParams$Parallel=="Yes")
    )

    tidyProj <- covidseir::tidy_seir(theProj, resample_y_rep = projParams$Resampling)
    tidyProj <- dplyr::left_join(tidyProj, lut, by = "day")

    for(projRow in 1:nrow(theProj))
    {
        simData <- addRow(simData, value=c(
            Variable = "Cases",
            Jurisdiction = "Canada - British Columbia",
            Timestep = tidyProj[projRow,]$date,
            Value = tidyProj[projRow,]$y_rep_mean,
            Iteration = index
        ))
    }

    # startRow <- min(which(is.na(simData$Iteration)))
    # simData[startRow:(startRow+totalDuration), Iteration:=index]
    # simData[startRow:(startRow+totalDuration), Date:=tidyProj$date]
    # simData[startRow:(startRow+totalDuration), simCases:=tidyProj$y_rep_mean]

}

simData <-transform(simData, TImestep = as.character(Timestep))

saveDatasheet(myScenario, simData, name = "epi_DataSummary", append = F)
saveDatasheet(myScenario, simData, name = "modelCovidseir_SimData", append = F)


# simData <-na.omit(simData)
# theFinal <- datasheet(myScenario, "epi_DataSummary", empty=T)
# theFinal[nrow(simData),] <- NA
# theFinal$Variable <- "Cases"
# theFinal$Jurisdiction <- "Canada - British Columbia"
# theFinal$Value <-simData$simCases
# theFinal$Timestep <-simData$Date
# theFinal$Iteration <-simData$Iteration
# saveDatasheet(myScenario, theFinal, "epi_DataSummary")
#
# simDataFilename <-sprintf("%s\\simData.csv", env$TransferDirectory)
# write.csv(simData, simDataFilename, row.names=F)
