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

runMin <- projParams$MinIteration
runMax <- projParams$MaxIteration

fitFilename <- sprintf("%s\\%s.rds", env$TransferDirectory, genParams$RDS)
theFit <- readRDS(fitFilename)

projFilename <- paste(env$TransferDirectory, "proj.rda", sep='\\')

day_start_change <- 1;
proj2 <- covidseir::project_seir(
    theFit, # the fit object loaded
    iter = runMin:runMax, # MCMC iterations
    forecast_days = projParams$DaysProject, # number of days into the future to predict
    f_fixed_start = max(theFit$days) + projParams$StartChange, # the day at which to start changing f
    f_multi = rep(1.1, projParams$DaysProject - projParams$StartChange + 1),
    f_multi_seg = projParams$FSegment, # which f segment to use
    parallel = (projParams$Parallel=="Yes")
);

theProj <- covidseir::project_seir(theFit, iter = runMin:runMax)

projectionFilename <- sprintf("%s\\%s_projection.rds", env$TransferDirectory, genParams$RDS)
saveRDS(theProj,projectionFilename)

tidyProj <- covidseir::tidy_seir(theProj, resample_y_rep = 20)
first_day <- min(caseData$Timestep, na.rm=TRUE)
last_day <- nrow(caseData) # + days_to_project
lut <- dplyr::tibble(
    day = seq_len(last_day),
    date = seq(first_day, first_day + length(day) - 1, by = "1 day")
)
tidyProj <- dplyr::left_join(tidyProj, lut, by = "day")

tidyProjFilename <- sprintf("%s\\%s_tidyproj.rds", env$TransferDirectory, genParams$RDS)
saveRDS(tidyProj, file=tidyProjFilename)

projTable <- datasheet(myScenario, "modelCovidseir_TidyProj")
projTable[nrow(tidyProj),] <- NA

projTable$Date <- tidyProj$date
projTable$Day <- tidyProj$day

projTable$YRepMean <- tidyProj$y_rep_mean
projTable$YRep05 <- tidyProj$y_rep_0.05
projTable$YRep25 <- tidyProj$y_rep_0.25
projTable$YRep50 <- tidyProj$y_rep_0.50
projTable$YRep75 <- tidyProj$y_rep_0.75
projTable$YRep95 <- tidyProj$y_rep_0.95

projTable$MuMean <- tidyProj$mu_mean
projTable$Mu05 <- tidyProj$mu_0.05
projTable$Mu25 <- tidyProj$mu_0.05
projTable$Mu50 <- tidyProj$mu_0.05
projTable$Mu75 <- tidyProj$mu_0.05
projTable$Mu95 <- tidyProj$mu_0.05

saveDatasheet(myScenario, projTable, "modelCovidseir_TidyProj")
