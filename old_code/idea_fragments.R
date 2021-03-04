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

# smooth Weibull curve to be superimposed
weibullX <- seq(min(delayData$time_to_report), max(delayData$time_to_report), by=0.05)
theDist <- data.table(x = weibullX, y = nrow(delayData)*dweibull(weibullX, shape=mleShape, scale=mleScale))

weibullDist <- datasheet(myScenario, "modelCovidseir_DelayDist")
weibullDist[length(weibullX),] <-NA
weibullDist$x <- theDist$x
weibullDist$y <- theDist$y
saveDatasheet(myScenario, weibullDist, "modelCovidseir_DelayDist")

subset(
       "2020-08-01" < symptom_onset_date &
           reported_date < "2020-12-01"
   )

   <!-- <datafeed name="FitRawOutput" displayName="Raw Output" dataScope="Scenario" isOutput="True">
     <datasheets>
       <datasheet name="TidyRawOutput" displayName="Raw Output">
         <columns>
           <column name="ScenarioID" dataType="Integer"/>
           <column name="Date" dataType="Date" displayName="Date"/>
           <column name="Day" dataType="Integer" displayName="Days elapsed"/>
           <column name="YRepMean" dataType="Double" displayName=""/>
           <column name="YRep05" dataType="Double" displayName=""/>
           <column name="YRep25" dataType="Double" displayName=""/>
           <column name="YRep50" dataType="Double" displayName=""/>
           <column name="YRep75" dataType="Double" displayName=""/>
           <column name="YRep95" dataType="Double" displayName=""/>
           <column name="MuMean" dataType="Double" displayName=""/>
           <column name="Mu05" dataType="Double" displayName=""/>
           <column name="Mu25" dataType="Double" displayName=""/>
           <column name="Mu50" dataType="Double" displayName=""/>
           <column name="Mu75" dataType="Double" displayName=""/>
           <column name="Mu95" dataType="Double" displayName=""/>
         </columns>
       </datasheet>
     </datasheets>
   </datafeed> -->

   <column name="FSeg1Post" dataType="Double" displayName="Contact rate fraction 1"/>
   <column name="FSeg2Post" dataType="Double" displayName="Contact rate fraction 2"/>
   <column name="FSeg3Post" dataType="Double" displayName="Contact rate fraction 3"/>
   <column name="FSeg4Post" dataType="Double" displayName="Contact rate fraction 4"/>


   datafeed name="FitFileInfo" displayName="Fit Data RDS" dataScope="Scenario">
   	<datasheets>
   		<datasheet name="FitFileInfo" displayName="Fit File Data" isSingleRow="True">
   			<columns>
   				<column name="InputDatasheetID" dataType="Integer" isPrimary="True"/>
   				<column name="ScenarioID" dataType="Integer"/>
   				<column name="FitDataFile" displayName="Fit Information Object File" dataType="String" isExternalFile="True"/>
   				<column name="MadeDateTime" displayName="Download Date/Time" dataType="Date"/>
   			</columns>
   		</datasheet>
   	</datasheets>
   </datafeed>

   <!-- <item name="Variable" datasheet="SimData" column="Value" variableSourceColumn="Variable" prefixFolderName="False" appendTo="epi_Variables"/> -->
   <!-- <layout name="modelCovidseir_Charts" configurationSheet="RunControl" xAxisLabelFormat="yyyy-MMM-dd">
   	<group name="Variables" dataSheet="DataSummary" filter="Jurisdiction">
   		<item name="Variable" dataSheet="DataSummary" column="Value" variableSourceColumn="Variable" prefixFolderName="False"/>
   	</group>
   </layout> -->

   # simData <- datasheet(myScenario, name = "epi_DataSummary", empty = T, optional = T, lookupsAsFactors = F)
# simData <- transform(simData, Timestep = as.Date(Timestep))
#
# simData <- data.table(matrix(ncol=8, nrow=maxIteration*totalDuration))
# names(simData) <-
# simData[, Iteration:=as.integer(Iteration)]
