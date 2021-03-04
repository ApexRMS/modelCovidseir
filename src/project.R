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

maxIteration <- runControl$MaximumIteration
if(length(maxIteration)==0){ maxIteration <- 10 }
if(is.na(maxIteration)){ maxIteration <- 10 }

caseData <- datasheet(myScenario, "epi_DataSummary") %>% transform(Timestep = as.Date(Timestep))

maxTimeStep <- runControl$MaximumTimestep
if(length(maxTimeStep)==0){ maxTimeStep <- max(caseData$Timestep) + 45 }
if(is.na(maxTimeStep)){ maxTimeStep <- max(caseData$Timestep) + 45 }

totalDuration <- (as.Date(maxTimeStep) - min(caseData$Timestep))[[1]]
daysToProject <-(as.Date(maxTimeStep) - as.Date(max(caseData$Timestep)))[[1]]

if((daysToProject<=0) | (is.na(daysToProject))) stop()

theFit <- readRDS(sprintf("%s\\%s.rds", env$TempDirectory, genParams$RDS))

projParams <- datasheet(myScenario, "modelCovidseir_ProjParams")

# in the case of a huge number of iterations,we'll have a table already set up
firstDay <- min(caseData$Timestep, na.rm=TRUE)
lastDay <- totalDuration + 10
lut <- dplyr::tibble(
    day = seq_len(lastDay),
    date = seq(firstDay, firstDay + length(day) - 1, by = "1 day")
)

simData <- datasheet(myScenario, "epi_DataSummary", empty = T, optional = T, lookupsAsFactors = F)
# simData$Iteration <- numeric(0)
# simData$Timestep <- numeric(0)
simData <- simData %>% transform(Timestep=as.Date(Timestep))
simData$AgeMin <- NULL
simData$AgeMax <- NULL
simData$Sex <- NULL

envBeginSimulation(maxIteration)

for(index in 1:maxIteration)
{
    envReportProgress(1, index)
    theProj <- covidseir::project_seir(
        theFit,
        iter = 1:1,
        forecast_days = daysToProject,
        f_fixed_start = max(theFit$days) + projParams$StartChange,
        f_multi = rep(1.1, daysToProject - projParams$StartChange + 1),
        f_multi_seg = projParams$FSegment,
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
    envStepSimulation()
}

envEndSimulation()

simData <- transform(simData, Timestep=as.IDate(Timestep))

saveDatasheet(myScenario, simData, name = "epi_DataSummary", append = F)
