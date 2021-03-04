library(tidyr)
library(data.table)
library(dplyr)
library(rsyncrosim)

library(covidseir)

myScenario <- scenario()
env <- ssimEnvironment()

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

wParams <- datasheet(myScenario, "modelCovidseir_WeibullParameters")

numIter <- runControl$MaximumIteration
if(length(numIter) == 0){ numIter <- 2000 }
if(is.na(numIter)){ numIter <- 2000 }

theFit <- covidseir::fit_seir(
    obs_model="NB2",
    daily_cases = caseData$Value,
    samp_frac_fixed = proportionTested,
    f_seg = fSeg,
    # f_prior = cbind(fSegments$PriorMean, fSegments$PriorSD), # MUST BE A MATRIX, NOT A TABLE
    R0_prior =  c(log(genParams$R0PriorMean), genParams$R0PriorSD),
    e_prior = c(genParams$EPriorMean, genParams$EPriorSD),
    start_decline_prior = c(log(genParams$StartDeclinePriorMean), genParams$StartDeclinePriorSD),
    end_decline_prior = c(log(genParams$EndDeclinePriorMean), genParams$EndDeclinePriorSD),
    chains = 4,
    iter = numIter,
    N_pop = genParams$NPop,
    i0_prior = c(log(genParams$I0PriorMean), genParams$I0PriorSD),
    delay_shape = wParams$mleShape,
    delay_scale = wParams$mleScale,
    time_increment = 0.25,
    save_state_predictions = TRUE,
    fit_type = if(genParams$FitType==1) "NUTS" else if(genParams$FitType==2) "VB" else "optimizing"
)

fitFilename <- sprintf("%s\\%s.rds", env$TempDirectory, genParams$RDS)
saveRDS(theFit, file=fitFilename)

fitFileInfo <- datasheet(myScenario, "modelCovidseir_FitFileInfo")
fitFileInfo[1,] <- NA
fitFileInfo$FitDataFile <-fitFilename
fitFileInfo$MadeDateTime <- as.character(Sys.time())
saveDatasheet(myScenario, fitFileInfo, "modelCovidseir_FitFileInfo")

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

contactPosts <- datasheet(myScenario, "modelCovidseir_PostsContactRates", empty=T)
contactPosts[numIter*nrow(fSegments)+1,] <- NA

rowCounter <- 1
for(segNumber in 1:ncol(theFit$post$f_s))
{
    for(i in 1:numIter)
    {
        contactPosts[rowCounter,] <- list(
            Iteration = i,
            Segment = segNumber,
            Day = fSegments[segNumber,]$BreakDay,
            Posterior = theFit$post$f_s[,segNumber][i]
        )
        rowCounter <- rowCounter + 1
    }
}

contactPosts <- na.omit(contactPosts)
saveDatasheet(myScenario, contactPosts,  "modelCovidseir_PostsContactRates")
