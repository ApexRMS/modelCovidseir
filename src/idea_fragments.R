############################################### FRAGMENTS ##########################################
# fragments after changes in direction or unaccepted ideas
####################################################################################################

# Reading segnemt data indexed by date instead (currently Integer days since start)
#
# fSegments <- data.table(datasheet(myScenario, "modelCovidseir_FSegments"))
fSegments$Mean <- sapply(fSegments$Mean, function(x) min(x,1))
fSegments$Date <- as.Date(fSegments$Date)
segmentsFilename <- paste(env$TransferDirectory, "SocialDistancingSegments.csv", sep='/')
write.csv(fSegments, segmentsFilename)

# sample data
fSegments <- data.table(Date=c("2019-01-01", "2020-10-01", "2020-11-15", "2020-12-15"), Mean=c(0.4, 0.95, 0.45, 0.7), SD=c(0.2, 0.2, 0.2, 0.2))
# creating the segments by date
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

# reading sampling fractions indexed by date instead (currently Integer days since start)
currentSegment <- 1
sFrac <- c()
for(date in caseData$Timestep)
{
    if(currentSegment < nrow(sampFrac))
    if(date == sampFrac[currentSegment+1, Date])
    {
        currentSegment <- currentSegment + 1
    }
    sFrac <-c(sFrac, sampFrac[currentSegment, Proportion])
}
fSeg[1] <- 0

# functions for getting the means and sd of the prior distributions given as comma separated strings
## the log means aren't looged before entry
getMean <- function(x){ return(as.numeric(strsplit(x, ',')[[1]])[1]) }
getLogMean <- function(x){ return(log(getMean(x))) }
getSD <- function(x){ return(as.numeric(strsplit(x, ',')[[1]])[2]) }

# automating sample fraction caluclation using BC CDC lab data and an underascertainment parameter
underascertainmentRatio <- 8
labData <- fread("http://www.bccdc.ca/Health-Info-Site/Documents/BCCDC_COVID19_Dashboard_Lab_Information.csv")

# check to make sure that both pieces of data have the same attached dates
sampFrac <- min(
    labData[Date %in% caseData$Timestep & Region=="BC", New_Tests*Positivity/100]/(underascertainmentRatio*caseData$Value),
    caseData$sampfrac
)

# broken code - in case the dates don't match up,check manually in a for loop
for(date in sapply(caseData$Timestep, function(x) toString(x)))
{
    print(date)
    perint()
    numer <- labData[(Region=="BC") & (Date==date), (New_Tests*Positivity/100)]
    denom <- caseData[Timestep==date, Value]*underascertainmentRatio

    caseData[Timestep==date, sampfrac] <- min(1, numer/denom)
}
