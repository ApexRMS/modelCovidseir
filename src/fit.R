# fit.R

# TODO - add header information

library(tidyr)
library(data.table)
library(dplyr)
library(rsyncrosim)

library(covidseir)

myScenario <- scenario()
env <- ssimEnvironment()

# get the case data

# TODO: run on multiple user-selected jurisdictions
# datasheet(myScenario, "modelCovidseir_FitRunJurisdictions")$Jurisdictions

dataSummary <- data.table(datasheet(myScenario, "epi_DataSummary"))
caseData <- dataSummary[Variable == 'Cases - Daily' & Jurisdiction == 'Canada - British Columbia', .SD, .SDcols=c("Timestep", "Value")][order(Timestep)]
caseData[, Timestep:=as.Date(Timestep)]

if(nrow(caseData) == 0){ stop() }

# estimated fraction of cases sampled
proportionTested <- c(rep(0.14, 13), rep(0.21, 38))
proportionTested <- c(proportionTested, rep(0.37, nrow(caseData) - length(proportionTested)))

# the starting segment is naturally #1
currentSegment <- 1
# the values of the f parameter in each segment, being estimations, are given as beta mean and sd pairs
beta_means_f_param <- c()
beta_sds_f_param <- c()

fSeg <- c(0, rep(currentSegment, nrow(caseData) - 1))
{
    # at the start of the simulation, say the distanced individuals contribute f=0.4 to the force of infection
    currentSegment <- currentSegment + 1
    beta_means_f_param <- c(beta_means_f_param, 0.4)
    beta_sds_f_param <- c(beta_sds_f_param, 0.2)
}
day_second_rise <- which(caseData$Timestep == lubridate::ymd("2020-10-01"))
if(length(day_second_rise))
{
    # if restrictions were loosened around 1st Oct, f goes up in this segment
    fSeg[seq(day_second_rise, length(fSeg))] <- currentSegment
    currentSegment <- currentSegment + 1
    beta_means_f_param <- c(beta_means_f_param, 0.95)
    beta_sds_f_param <- c(beta_sds_f_param, 0.2)
}
day_they_got_smarter <- which(caseData$Timestep == lubridate::ymd("2020-11-13"))
if(length(day_they_got_smarter))
{
    # case numbers rose meteorically, so restrictions were tightened on 15 Nov, bringing f down to 0.45
    fSeg[seq(day_they_got_smarter, length(fSeg))] <- currentSegment
    currentSegment <- currentSegment + 1
    beta_means_f_param <- c(beta_means_f_param, 0.6)
    beta_sds_f_param <- c(beta_sds_f_param, 0.2)
}
day_going_up_again <- which(caseData$Timestep == lubridate::ymd("2021-02-15"))
if(length(day_going_up_again))
{	# loosened restrictions on 15 Dec, but not like before - f rises to 0.7
    fSeg[seq(day_going_up_again, length(fSeg))] <- currentSegment
    currentSegment <- currentSegment + 1
    beta_means_f_param <- c(beta_means_f_param, 0.7)
    beta_sds_f_param <- c(beta_sds_f_param, 0.2)
}

genParams <- datasheet(myScenario, "modelCovidseir_General")
runControl <- datasheet(myScenario, "modelCovidseir_FitRunSettings")
wParams <- datasheet(myScenario, "modelCovidseir_DelayWeibull")

numIter <- runControl$MaximumIteration
if(length(numIter) == 0){ numIter <- 2000 }
if(is.na(numIter)){ numIter <- 2000 }

theFit <- covidseir::fit_seir(
    daily_cases = caseData$Value,
    samp_frac_fixed = proportionTested,
    f_seg = fSeg,
    f_prior = cbind(fSegments$PriorMean, fSegments$PriorSD), # MUST BE A MATRIX, NOT A TABLE
    R0_prior =  c(log(genParams$R0PriorMean), genParams$R0PriorSD),
    e_prior = c(genParams$EPriorMean, genParams$EPriorSD),
    chains = 4,
    iter = numIter,
    N_pop = genParams$NPop,
    i0_prior = c(log(genParams$I0PriorMean), genParams$I0PriorSD),
    delay_shape = wParams$MleShape,
    delay_scale = wParams$MleScale,
    time_increment = 0.1,
    save_state_predictions = TRUE,
    fit_type = if(genParams$FitType==1) "NUTS" else if(genParams$FitType==2) "VB" else "optimizing"
)

# save the fit object to a file
fitFilename <- sprintf("%s\\%s.rds", env$TempDirectory, genParams$RDS)
saveRDS(theFit, file=fitFilename)

# log the time and location info for the output file
fitFileInfo <- datasheet(myScenario, "modelCovidseir_FitFileInfo")
fitFileInfo[1,] <- NA
fitFileInfo$FitDataFile <-fitFilename
fitFileInfo$MadeDateTime <- as.character(Sys.time())
saveDatasheet(myScenario, fitFileInfo, "modelCovidseir_FitFileInfo")

# save the posterior parameter values
genPosteriors <- datasheet(myScenario, "modelCovidseir_PostsGeneral")
genPosteriors[numIter,] <- NA
genPosteriors$Iteration <- 1:numIter
genPosteriors$I0Post <- theFit$post$i0
genPosteriors$R0Post <- theFit$post$R0
genPosteriors$EPost <- theFit$post$e
genPosteriors$StartDeclinePost <- theFit$post$start_decline
genPosteriors$EndDeclinePost <- theFit$post$end_decline
genPosteriors$PhiPost <- theFit$post$phi[,1]
saveDatasheet(myScenario, genPosteriors, "modelCovidseir_PostsGeneral")

# committing the contact rate posteriors to their own data sheet
contactPosts <- datasheet(myScenario, "modelCovidseir_PostsContactRates", empty=T)
contactPosts[numIter*nrow(fSegments)+1,] <- 0

rowCounter <- 1
for(segNumber in 1:ncol(theFit$post$f_s))
{
    for(i in 1:numIter)
    {
        contactPosts[rowCounter,] = list(
            Iteration = i,
            Segment = segNumber,
            Day = fSegments[segNumber,]$BreakDay,
            Posterior = theFit$post$f_s[,segNumber][i]
        )
        rowCounter <- rowCounter + 1
    }
}

contactPosts <- data.table(contactPosts)
contactPosts[, Day:=as.Date(Day, format = "%Y-%m-%d")]
contactPosts <- na.omit(contactPosts)
saveDatasheet(myScenario, contactPosts,  "modelCovidseir_PostsContactRates")

theProj <- covidseir::project_seir(theFit, iter = 1:50)
tidyProj <- covidseir::tidy_seir(theProj, resample_y_rep = 20)

