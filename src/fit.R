library(tidyr)
library(data.table)
library(dplyr)
library(rsyncrosim)

library(covidseir)

myScenario <- scenario()
env <- ssimEnvironment()

wParams <- fread(paste(env$TransferDirectory, "modelCovidseir_WeibullParameters.csv", sep="/"))

# caseData <- data.table(datasheet(myScenario, "epi_DataSummary"))

##################### TEMP CODE UNTIL EPI IS RECOMPILED AND RERELEASED ############################

caseData <- fread("http://www.bccdc.ca/Health-Info-Site/Documents/BCCDC_COVID19_Regional_Summary_Data.csv") %>%
    separate(col=Date, into=c("date"), sep=' ') %>%
    mutate_at(vars(date), as.IDate) %>%
    subset(HA=="All" & Province=="BC") %>%
    select(date, Cases_Reported) %>%
    rename(value=Cases_Reported) %>%
    rename(Timestep=date, Value=value) %>%
    dplyr::as_tibble()

###################################################################################################

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

theFit <- covidseir::fit_seir(
    obs_model="NB2",
    daily_cases = caseData$Value,
    samp_frac_fixed = proportionTested,
    f_seg = fSeg,
    R0_prior =  c(log(genParams$R0PriorMean), genParams$R0PriorSD),
    # f_prior = cbind(fSegments$PriorMean, fSegments$PriorSD), # MUST BE A MATRIX, NOT A TABLE
    e_prior = c(genParams$EPriorMean, genParams$EPriorSD),
    start_decline_prior = c(log(genParams$StartDeclinePriorMean), genParams$StartDeclinePriorSD),
    end_decline_prior = c(log(genParams$EndDeclinePriorMean), genParams$EndDeclinePriorSD),
    chains = 4,
    iter = 100,
    N_pop = genParams$NPop,
    i0_prior = c(log(genParams$I0PriorMean), genParams$I0PriorSD),
    delay_shape = wParams$mleShape,
    delay_scale = wParams$mleScale,
    time_increment = 0.1,
    save_state_predictions = TRUE,
    fit_type = if(genParams$FitType==1) "NUTS" else if(genParams$FitType==2) "VB" else "optimizing"
)

fitFilename <- sprintf("%s\\%s.rds", env$TransferDirectory, genParams$RDS)
saveRDS(theFit, file=fitFilename)
