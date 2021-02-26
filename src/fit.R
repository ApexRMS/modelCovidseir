library(tidyr)
library(data.table)
library(dplyr)
library(rsyncrosim)

library(covidseir)

# 	delayDatasheet - a table of case reporting data, each row representing a unique case
# 		<symptom onset date>	<case reporting date>	<gap (measured in days)>
#       Example:    2019-01-01  0.4     0.2
#                     2020-10-01  0.95    0.2
#                     2020-11-15  0.45    0.2
#                     2020-12-15  0.7     0.2

myScenario <- scenario()
env <- ssimEnvironment()

ymd <- lubridate::ymd

##################### DUMMY DATA #########################

fitParams <- data.table(NPop=5071000, TimeIncrement=100, NChains=2, NIterations=100, FitType="optimizing", MCMCSeed=42)
fitPriors <-data.table(R0Prior="2.6,0.2", EPrior="0.8,0.5", StartDecline="15,0.05", EndDecline="22, 0.05", I0Prior="8,1")

sampFrac <- c(rep(0.14, 13), rep(0.21, 38))
sampFrac <- c(sampFrac, rep(0.37, nrow(caseData) - length(sampFrac)))

wParams <- data.table(mleShape=1.6768, mleScale=5.9)

caseData <- datasheet(myScenario, "epi_DataSummary")
caseData$Timestep <- as.Date(caseData$Timestep)

fSegments <- data.table(Date=c("2019-01-01", "2020-10-01", "2020-11-15", "2020-12-15"), Mean=c(0.4, 0.95, 0.45, 0.7), SD=c(0.2, 0.2, 0.2, 0.2))

#######################################################


# fSegments <- data.table(datasheet(myScenario, "modelCovidseir_FSegments"))
fSegments$Mean <- sapply(fSegments$Mean, function(x) min(x,1))
fSegments$Date <- as.Date(fSegments$Date)
segmentsFilename <- paste(env$TransferDirectory, "SocialDistancingSegments.csv", sep='/')
write.csv(fSegments, segmentsFilename)

# is there a way of giving the user dropdown menus? I seen this in another package - it would be clutch for the Karlen model selsction too

currentSegment <- 1
fSeg <- c()
for(date in caseData$Timestep)
{
    if(currentSegment < nrow(fSegments))
    if(date == fSegments[currentSegment+1, Date])
    {
        currentSegment <- currentSegment + 1
    }
    fSeg <-c(fSeg, currentSegment)
}
fSeg[1] <- 0

# fitParams <-datasheet(myScenario, "modelCovidseir_FitParameters")
# fitPriors <- datasheet(myScenario, "modelCovidseir_FitPriors")
# # wParams <- datasheet(myScenario)

getMean <- function(x){ return(as.numeric(strsplit(x, ',')[[1]])[1]) }
getLogMean <- function(x){ return(log(getMean(x))) }
getSD <- function(x){ return(as.numeric(strsplit(x, ',')[[1]])[2]) }

fit <- covidseir::fit_seir(
    # obs_model="NB2",
    daily_cases = caseData$Value,
    samp_frac_fixed = sampFrac,
    f_seg = fSeg,
    R0_prior =  c(getLogMean(fitPriors$R0Prior), getSD(fitPriors$R0Prior)),
    f_prior = fSegments[, .SD, .SDcols=c("Mean", "SD")],
    e_prior = c(getMean(fitPriors$EPrior), getSD(fitPriors$EPrior)),
    start_decline_prior = c(getLogMean(fitPriors$StartDecline), getSD(fitPriors$StartDecline)),
    end_decline_prior = c(getLogMean(fitPriors$EndDecline), getSD(fitPriors$EndDecline)),
    chains = fitParams$NChains,
    iter = fitParams$NIterations,
    N_pop = fitParams$NPop,
    i0_prior = c(getLogMean(fitPriors$I0Prior), getSD(fitPriors$I0Prior)),
    delay_shape = wParams$mleShape, # shape_MLE,
    delay_scale = wParams$mleScale, # scale_MLE,
    time_increment = fitParams$TimeIncrement,
    #######################
    save_state_predictions = TRUE,
    fit_type = fitParams$FitType
)
