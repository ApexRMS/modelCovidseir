rm(list=ls());

setwd("C:/Users/User/Dropbox/samplepackage");

# library(data.table)
# library(dplyr)
#
# # table Differences showing that all case counts are not equal
# # I reckon we can check the regional count
# # not sure where the delay data is coming from, but it doesn['t account for all cases
#
# # download.file(
# #   "https://github.com/andrew-edwards/rightTruncation/raw/master/data/delay_data_2021_01_05.rda",
# #   destfile="delay_data.rda", quiet=TRUE
# # );
#
# data.table(fread("case_data.csv")) %>%
#     subset(Province=="BC" & HA=="All" & HSDA=="All") %>%
#     select(Date, Cases_Reported) %>%
#     rename(orig_New_Cases=Cases_Reported) -> BC_cases;
#
# data.table(fread("BCCDC_COVID19_Dashboard_Lab_Information.csv")) %>%
#     subset(Region=="BC" & Date>=min(BC_cases$Date)) %>%
#     mutate(New_Cases=ceiling(Positivity/100*New_Tests)) %>%
#     select("Date", "New_Tests", "New_Cases") %>%
#     rename(lab_New_Tests=New_Tests, lab_New_Cases=New_Cases) %>%
#     data.table() -> BC_lab;
#
# data.table(get(load("delay_data.rda"))) %>%
#     group_by(reported_date) %>%
#     summarise(delay_New_Cases=n(), .groups='drop') %>%
#     rename(Date=reported_date) %>%
#     data.table() -> delay_data;
#
# earliest_common_date <- max(min(BC_cases$Date), min(BC_lab$Date), min(delay_data$Date));
# latest_common_date <- min(max(BC_cases$Date), max(BC_lab$Date), max(delay_data$Date));
#
# BC_cases <- BC_cases[Date>=earliest_common_date & Date<=latest_common_date];
# BC_lab <- BC_lab[Date>=earliest_common_date & Date<=latest_common_date];
# delay_data <- delay_data[Date>=earliest_common_date & Date<=latest_common_date];
#
# data.table(Date=sort(unique(c(BC_cases$Date, BC_lab$Date, delay_data$Date)))) %>%
#     BC_cases[., on=.(Date), nomatch=0L] %>%
#     delay_data[., on=.(Date), nomatch=0L] %>%
#     BC_lab[., on=.(Date), nomatch=0L] -> Differences
#
# # searching the smooth time series of the case data for change points from various tests and seeing how they sqaure with Jens' MountainMath data

library(trend)
library(ggplot2)

data.table(fread("case_data.csv")) %>%
    subset(Province=="BC" & HA=="All" & HSDA=="All") %>%
    select(Date, Cases_Reported_Smoothed) %>%
    subset(Date>="2020-08-01") %>%
    rename(Count=Cases_Reported_Smoothed) -> BC_cases;

change_objs <- c();
complete <- FALSE;

skip <- 20;
TS <- rollapply(BC_cases$Count, width=skip, mean, align="right");

i <- 40;
j <- i + 15;

while(j <= length(TS))
{
    # writeLines(sprintf("\ti=%i, j=%i", i, j));

    obj <- br.test(TS[i:j]);

    if(obj$p.value<0.05 & obj$p.value>0)
    {
        # writeLines("\t\tgot object")
        change_objs <- c(change_objs, obj$estimate[[1]]+i+skip);
        i <- j;
    }

    if(j>=length(TS)){ break;}

    j <- j + 1;
}


pl_time_series <- ggplot(BC_cases, aes(Date, Count)) +
    geom_line(colour="black", size=1) +
    theme_bw() +
    scale_y_continuous(expand=c(0,0)) +
    theme(
        axis.title = element_text(size=15),
        axis.text = element_text(size=10),
    );

# calculating the gradient
# green for positive reopening events, red for restrictions
point_colour <- function(point)
{
    grad <- (BC_cases[poi+1]$Count - BC_cases[poi-1]$Count)/(BC_cases[poi+1]$Date - BC_cases[poi-1]$Date)
    print(poi+skip)
    if(grad > 0) return("red");
    if(grad < 0) return("green");
    return("black");
}


for(poi in change_objs)
{
    pl_time_series <- pl_time_series +  geom_vline(xintercept=BC_cases[poi]$Date, linetype="dashed", colour=point_colour(poi+skip), size=1)
}

# pl_time_series <- ggplot(BC_cases, aes(Date, Count)) +

# pl_time_series <- ggplot(data.table(1:length(TS), TS), aes(x=V1, y=TS)) +
#     geom_line(colour="black", size=1) +
#     theme_bw() +
#     scale_y_continuous(expand=c(0,0)) +
#     theme(
#         axis.title = element_text(size=15),
#         axis.text = element_text(size=10),
#     );
#
# # calculating the gradient
# # green for positive reopening events, red for restrictions
# point_colour <- function(point)
# {
#     grad <- (BC_cases[poi+1]$Count - BC_cases[poi-1]$Count)/(BC_cases[poi+1]$Date - BC_cases[poi-1]$Date)
#     print(poi+skip)
#     if(grad > 0) return("red");
#     if(grad < 0) return("green");
#     return("black");
# }
#
#
# # for(poi in change_objs)
# # {
# #     pl_time_series <- pl_time_series +  geom_vline(xintercept=BC_cases[poi+skip]$Date, linetype="dashed", colour=point_colour(poi+skip), size=1)
# # }
#
# plot(pl_time_series)


plot(pl_time_series)


