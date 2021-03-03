library(tidyr)
library(data.table)
library(dplyr)
library(rsyncrosim)

library(covidseir)

myScenario <- scenario()
env <- ssimEnvironment()

wParams <- datasheet(myScenario, "modelCovidseir_WeibullParameters")

caseData <- data.table(datasheet(myScenario, "epi_DataSummary"))
caseData$Timestep <- as.Date(caseData$Timestep)

fSegments <- datasheet(myScenario, "modelCovidseir_ContactRateFractions") %>% arrange(BreakDay) %>% data.table

fSeg <- c()
currentSegment <- 1
for(daysSince0 in 1:nrow(caseData))
{
    if(currentSegment < nrow(fSegments)){
    if(daysSince0 == fSegments[currentSegment+1, BreakDay])
    {  currentSegment <- currentSegment + 1 }}

    fSeg <- c(fSeg, currentSegment)
}
fSeg[1] <- 0

sampFrac <- datasheet(myScenario, "modelCovidseir_SamplingFractions") %>% arrange(Day) %>% data.table()

proportionTested <- c()
currentRow <- 1
for(daysSince0 in 1:nrow(caseData))
{
    if(currentRow < nrow(sampFrac))
    if(daysSince0 == sampFrac[currentRow+1, Day])
    { currentRow <- currentRow+1 }

    proportionTested <- c(proportionTested, sampFrac[currentRow, Proportion])
}

genParams <- datasheet(myScenario, "modelCovidseir_General")

runControl <- datasheet(myScenario, "epi_RunControl")

theFit <- covidseir::fit_seir(
    obs_model="NB2",
    daily_cases = caseData$Value,
    samp_frac_fixed = proportionTested,
    # f_seg = fSeg,
    # f_prior = cbind(fSegments$PriorMean, fSegments$PriorSD), # MUST BE A MATRIX, NOT A TABLE
    R0_prior =  c(log(genParams$R0PriorMean), genParams$R0PriorSD),
    e_prior = c(genParams$EPriorMean, genParams$EPriorSD),
    start_decline_prior = c(log(genParams$StartDeclinePriorMean), genParams$StartDeclinePriorSD),
    end_decline_prior = c(log(genParams$EndDeclinePriorMean), genParams$EndDeclinePriorSD),
    chains = 4,
    iter = runControl$MaximumIteration,
    N_pop = genParams$NPop,
    i0_prior = c(log(genParams$I0PriorMean), genParams$I0PriorSD),
    delay_shape = wParams$mleShape,
    delay_scale = wParams$mleScale,
    time_increment = 0.1,
    save_state_predictions = TRUE,
    fit_type = if(genParams$FitType==1) "NUTS" else if(genParams$FitType==2) "VB" else "optimizing"
)

fitFilename <- sprintf("%s\\%s.rds", env$TempDirectory, genParams$RDS)
saveRDS(theFit, file=fitFilename)

fitPosteriors <- datasheet(myScenario, "modelCovidseir_FitPosteriors")
fitPosteriors[runControl$MaximumIteration,] <- NA
fitPosteriors$Iteration <- 1:runControl$MaximumIteration
fitPosteriors$I0Post <- theFit$post$i0
fitPosteriors$R0Post <- theFit$post$R0
fitPosteriors$EPost <- theFit$post$e
fitPosteriors$StartDeclinePost <- theFit$post$start_decline
fitPosteriors$EndDeclinePost <- theFit$post$end_decline
fitPosteriors$PhiPost <- theFit$post$phi[,1]
fitPosteriors$FSeg1Post <- theFit$post$f_s[,1]
# fitPosteriors$FSeg2Post <-theFit$post$f_s[,2]
# fitPosteriors$FSeg3Post <-theFit$post$f_s[,3]
# fitPosteriors$FSeg4Post <-theFit$post$f_s[,4]
fitPosteriors$FSeg2Post <- NULL
fitPosteriors$FSeg3Post <- NULL
fitPosteriors$FSeg4Post <- NULL
saveDatasheet(myScenario, fitPosteriors, "modelCovidseir_FitPosteriors")
