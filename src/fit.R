# if the rightTruncation and covidseir R packages aren't installed, use the following code
# library(devtools)
# library(pkgbuild)
# library(remotes)
#
# if(!require("rightTruncation")){ devtools::install_github("andrew-edwards/rightTruncation") }
# if(!require("covidseir")){ remotes::install_github("seananderson/covidseir", build_vignettes = TRUE) }
#-------------------------------------------------------------------------------------------------------

library(tidyr)
library(data.table)
library(dplyr)
library(rsyncrosim)

library(rightTruncation)
library(covidseir)

# is there a way of giving the user dropdown menus? I seen this in another package - it would be clutch for the Karlen model selsction too

TYPE_OF_FIT <- "optimizing"

# set the first day so that we can calculate the time that passed in days
hDay0 <- min(delayData$symptom_onset_date)
# set the number of rows and columns of the matrix
numRows <- (max(delayData$symptom_onset_date)-hDay0)[[1]]
numCols <- (max(delayData$reported_date)-hDay0)[[1]]

# fill the Hnr matrix with the case data as described in the previous long comment
Hnr <- matrix(0, nrow=numRows, ncol=numCols)
for(i in 1:nrow(delayData))
{
	onset <- (delayData[i]$symptom_onset_date - hDay0)[[1]]
	repor <- (delayData[i]$reported_date - hDay0)[[1]];
	Hnr[onset, repor] <- Hnr[onset, repor]+1
}

mleRes = suppressWarnings(nlm(
	f = negLL_Weibull_counts_matrix,
	p = c(2, 7),
	h_nr = Hnr
))

# the parameters for the distribution
mleShape <- mleRes$estimate[1]
mleScale <- mleRes$estimate[2]

wParams <- datasheet(currentScenario, "modelCovidseir_WeibullParameters")
wParams[1,] <-NA
wParams$mleShape <- mleShape
wParams$mleScale <- mleScale
saveDatasheet(currentScenario, wParams, "modelCovidseir_WeibullParameters")

# smooth Weibull curve to be superimposed
weibullX <- seq(min(delayData$time_to_report), max(delayData$time_to_report), by=0.05)
theDist <- data.table(x = weibullX, y = nrow(delayData)*dweib ull(weibullX, shape=mleShape, scale=mleScale))

weibullDist <- datasheet(currentScenario, "modelCovidseir_DelayDist")
weibullDist[length(weibullX),] <-NA
weibullDist$x <- theDist$x
weibullDist$y <- theDist$y
saveDatasheet(currentScenario, weibullDist, "modelCovidseir_DelayDist")
