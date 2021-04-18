library(tidyr)
library(data.table)
library(dplyr)
library(rsyncrosim)

library(covidseir)

# allowing parallel processing
future::plan(future::multisession)

myScenario <- scenario()
env <- ssimEnvironment()

genParams <- datasheet(myScenario, "modelCovidseir_General")
projParams <- datasheet(myScenario, "modelCovidseir_ProjParams")
# runControl <- datasheet(myScenario, "epi_RunControl")

# if max iteration not given in the table by the user, we set it here
# maxIteration <- runControl$MaximumIteration
# if(length(maxIteration)==0){ maxIteration <- 10 }
# if(is.na(maxIteration)){ maxIteration <- 10 }
maxIteration <- 20

# BC case data
caseData <- data.table(datasheet(myScenario, "epi_DataSummary"))
caseData <- caseData[Variable == "Cases - Daily", .SD, .SDcols=c("Timestep", "Value")][order(Timestep)]
caseData[, Timestep:=as.Date(Timestep)]

# if the end day of the projection is not given by the user, default to 45 days
# maxTimeStep <- runControl$MaximumTimestep
# minTimeStep <- runControl$MinimumTimestep

# # TODO - sort out minimum timestep - should check for valid value from user
# minTimeStep <- max(caseData$Timestep)
# 
# if(length(maxTimeStep)==0){ maxTimeStep <- max(caseData$Timestep) + 28 }
# if(is.na(maxTimeStep)){ maxTimeStep <- max(caseData$Timestep) + 28 }

maxTimeStep <- max(caseData$Timestep) + 28

# calculate the duration and extension of the simulation
totalDuration <- (as.Date(maxTimeStep) - min(caseData$Timestep))[[1]]
daysToProject <-(as.Date(maxTimeStep) - as.Date(max(caseData$Timestep)))[[1]]

# if no projection is requested, fail
if((daysToProject<=0) | (is.na(daysToProject))) stop()

# importing the fit object from the fitting step
theFit <- readRDS(sprintf("%s\\%s.rds", env$TempDirectory, genParams$RDS))

# parameters used for the projection function
projParams <- datasheet(myScenario, "modelCovidseir_ProjParams")

# in the case of a huge number of iterations,we'll have a table already set up
firstDay <- min(caseData$Timestep, na.rm=TRUE)
lastDay <- totalDuration + 10
lut <- dplyr::tibble(
    day = seq_len(lastDay),
    date = seq(firstDay, firstDay + length(day) - 1, by = "1 day")
)

# read the standard output table
simData <- data.table(Iteration=numeric(), Timestep=character(), Variable=character(), Jurisdiction=character(), Value=numeric(), TransformerID=character())
simData[, Timestep:=as.Date(Timestep)]

theProj <- covidseir::project_seir(
    theFit,
    iter = 1:maxIteration,
    forecast_days = daysToProject,
    f_fixed_start = max(theFit$days) + projParams$StartChange,
    f_multi = rep(projParams$Multiplic, daysToProject - projParams$StartChange + 1),
    f_multi_seg = projParams$FSegment,
    imported_cases = projParams$Imports,
    imported_window = projParams$ImportWindow,
    parallel = (projParams$Parallel=="Yes")
)

epiDatasummary <- theProj %>% 
    mutate(Timestep=day+min(caseData$Timestep)) %>% 
    select(Timestep, y_rep, .iteration) %>% 
    rename(Iteration=.iteration, Value=y_rep) %>% 
    mutate(Jurisdiction="Canada - British Columbia") %>% 
    mutate(Variable = "Cases - Daily") %>% 
    mutate(TransformerID = "covidseir model: Project") %>%
    data.table()

saveDatasheet(myScenario, epiDatasummary, name = "epi_DataSummary")
