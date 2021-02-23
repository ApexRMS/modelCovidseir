rm(list=ls());

library(rsyncrosim)
library(covidseir)
library(rightTruncation)
library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)

ymd <- lubridate::ymd

options(width=Sys.getenv("COLUMNS"));

# write.csv(fread("http://www.bccdc.ca/Health-Info-Site/Documents/BCCDC_COVID19_Regional_Summary_Data.csv"), file="Temp_Cases.csv", row.names=FALSE);
# write.csv(fread("http://www.bccdc.ca/Health-Info-Site/Documents/BCCDC_COVID19_Dashboard_Lab_Information.csv"), file="Testing_Temp.csv", row.names=FALSE);
#
# download.file("https://github.com/andrew-edwards/rightTruncation/raw/master/data/delay_data_2021_01_05.rda", destfile="delay_data.rda")

# delay_data <- data.table(get(load("delay_data.rda"))) %>%
# 				select(-time_to_report) %>%
# 				mutate_if(is.character, as.Date) %>%
# 				subset(symptom_onset_date>="2020-03-15");
#
# day0 <- min(delay_data$symptom_onset_date);
#
# num_rows <- (max(delay_data$symptom_onset_date)-day0)[[1]];
# num_cols <- (max(delay_data$reported_date)-day0)[[1]];
#
# start_time <- Sys.time()
#
# H_nr <- matrix(0, nrow=num_rows, ncol=num_cols);
# for(i in 1:nrow(delay_data))
# {
# 	onset <- (delay_data[i]$symptom_onset_date - day0)[[1]];
# 	repor <- (delay_data[i]$reported_date - day0)[[1]];
# 	H_nr[onset, repor] <- H_nr[onset, repor]+1;
# }
#
# MLE_res = nlm(f = negLL_Weibull_counts_matrix, p = c(3, 15), h_nr = H_nr);
#
# k_MLE <- MLE_res$estimate[1]
# lambda_MLE <- MLE_res$estimate[2]

Temp_Cases <- data.table(fread("Temp_Cases.csv"));
Testing_Temp <- data.table(fread("Testing_Temp.csv"))

# get all the cases for BC
BC_cases <- Temp_Cases %>%
	subset(HA=="All" & Province=="BC") %>%
	separate(col=Date, into=c("Day"), sep=' ')

BC_cases[, Day:=as.Date(Day)];

BC_testing <- Testing_Temp # %>%
	# subset(Date>="2020-03-01");

BC_testing <- BC_testing[Region=="BC"];
names(BC_testing)[which(names(BC_testing)=="Date")] <- "Day";

BC_testing[, Population:=5071336];
BC_testing[, Untested:=Population-Total_Tests]; # assuming tha teveryone is tested once, so one test per person
BC_testing[, Turn_Around:=NULL];
BC_testing[, P_tested:=New_Tests/Untested];
BC_testing[, P_infected:=cumsum(Positivity*New_Tests/100)*8/Population]; # adjusting for underreporting
# BC_testing[, P_infected_given_tested:=cumsum(Positivity*New_Tests/100)/Total_Tests];
BC_testing[, P_infected_given_tested:=Positivity/100];
BC_testing[, P_tested_given_infected:=(P_infected_given_tested*P_tested)/P_infected];

stop()


BC_dat <- merge(BC_cases, BC_testing, by="Day");
BC_dat[, c("HA", "HSDA", "Region", "Turn_Around", "Total_Tests", "Cases_Reported_Smoothed"):=NULL];
names(BC_dat)[which(names(BC_dat)=="New_Tests")] <- "Tests_Done";

stop()



stop()

BC_dat <- tibble(BC_dat[Day>"2020-03-01"] );

# by-day estimate of the fraction of positive cases sampled/detected
# Based on estimation with hospital data in other model:
samp_frac <- c(rep(0.14, 13), rep(0.21, 38))
samp_frac <- c(samp_frac, rep(0.37, nrow(BC_dat) - length(samp_frac)));

# modelling the imposition and commission of contact rates
# basically functions as a list of indices of some vector of values
# instead of having a vector of rates c(1.2, 2.1, 3.4, 4.7, 4.7, 3.4, 3.4)
# let r1=1.2, r2=2.1, r3=3.4, r4=4.7 and then say
# c(r1, r2, r3, r4, r4, r3, r3)
f_seg <- c(0, rep(1, nrow(BC_dat) - 1))
day_new_f <- which(BC_dat$Day == ymd("2020-05-01"))
f_seg[seq(day_new_f, length(f_seg))] <- 2
day_ch <- which(BC_dat$Day == ymd("2020-06-01"))
f_seg[seq(day_ch, length(f_seg))] <- 3
# so this vector represents the presnece of three diferent contact rate fraction

start_time <- Sys.time();

# chains=4

fit <- covidseir::fit_seir(
	daily_cases = BC_cases$Cases_Reported, # the daily cases from the input matrix without the attached date - I reckon the dates just go into the graph at the end
	# obs_model=c('NB2), "Poisson"), # type of obersation model
	# forecast_days = 0, # number of days into the future to forecast - suggested to use zero here and project_seir() for predictions
	# time_increment = 0.25, # time increment for the ODEs and Weibull delay-model integration units. larger = faster asnd less accurate
	samp_frac_fixed = samp_frac, # a vector of sampled fractions
	# samp_frac_type = c("fixed", "estimated", "rw", "segmented"), # fixed, estimated, contrained random walk or segmented
	# samp_frac_seg = NULLl # vector of sample fraction segment indexes
	f_seg = f_seg, # a vector of segment indexes. should start at 0 to represent the "no social distancing" value f0, isn't used until after the end_decline day
	# days_back = 45, # number of days to go back for the Weibull case-delay integration. large enough that the results don't change
	R0_prior = c(log(2.6), 0.2), # lognormal log men and sd for R0
	# phi_prior = 1,  SD of ‘1/sqrt(phi) ~ Normal(0, SD)’ prior, where NB2(mu, phi) and ‘Var(Y) = mu + mu^2 / phi’. See the Stan Prior Choice Recommendations. If of length 2 will be treated as lognormal prior on phi.
	f_prior = cbind(c(0.4, 0.5, 0.6), c(0.2, 0.2, 0.2)), # beta mean and SD for the f parameters
	e_prior = c(0.8, 0.05), # beta mean and sd for the e parameter - the fraction of people who are social distancing
	# samp_freq_prior = c(0.4, 0.2),
	start_decline_prior = c(log(15), 0.1), # lognormal log mean and sd for the parameter reporesenting the day that osocial distancing starts ramping up
	end_decline_prior = c(log(22), 0.1), # lognormal log mean and sd for the day that the phyical distancing ramp finishes
	# use_ramp = TRUE, # logical. use the ramp?
	# rw_sigma = 0.1, # fised sd of the optional samp_frac random walks
	# seed = 42, # MCMC seed for rstan::stan()
	chains = 4, # number of MCMC chains for rstan::stan()
	iter = 500, # number of posterior samples
	N_pop = 5.1e6, # BC population
	# pars = c(D = 5, k1 = 1/5, k2 = 1, q = 0.05, ud = 0.1, ur = 0.02, f0 = 1), # named numeric vector of fixed parameter values
	i0_prior = c(log(8), 1), # lognormal log mean and sd of the number people infected at the initial point
	# state_0 = c(E1_frac = 0.4, E2_frac = 0.1, I_frac = 0.5, Q_num = 0, R_num = 0, E1d_frac = 0.4, E2d_frac = 0.1, Id_frac = 0.5, Qd_num = 0, Rd_num = 0), # initial state of the ODE
	# save_state_predictions = FALSE, # whether to include the state predictions in the results
	# delay_scale = 9.85, # Weibull scale parameter for the delay in reporting
	# delay_shape = 1.73, # Weibull shape parameter for the delay in reporting
	# ode_control = c(1e-07, 1e-06, 1e+06), # control options for the Stan ODE solver: relative difference, absolute difference, max iterations
	# fit_type = "optimizing" # Stan sampling/fitting algorithm to use
	# init = c("prior_random", "optimising"), # initialisation type (draw random;y from the prior or try to use the MAP estimmate)
	# init_list = NULL,
	# X= NULL
);

print(Sys.time() - start_time);

# choose to use parallel prcessing for the projections
future::plan(future::multisession)

# generate a lut for adding the dates back into the fit before we plot it
first_day <- day0;
last_day <- 300 # how many days to create dates for
lut <- dplyr::tibble(
	day = seq_len(last_day),
	date = seq(first_day, first_day + length(day) - 1, by = "1 day")
);

# ##############################################################################
# VISUALIZING THE FIT OF THE MODEL
# ##############################################################################

# take the fitted mode and calculate projections
proj_for_viz <- covidseir::project_seir(fit, iter = 1:50);

# take the posterior samples, resample from the negative bninomial observation model for smoother predictions and transform the output into a tidy data frame for plotting
tidy_proj_for_viz <- covidseir::tidy_seir(proj_for_viz, resample_y_rep = 20);

# # join a date column back on by creating a lut - the current table only has columns for the numberic days since the simulation started
first_day <- min(dat$date);
last_day <- 300 # how many days to create dates for
lut <- dplyr::tibble(
	day = seq_len(last_day),
	date = seq(first_day, first_day + length(day) - 1, by = "1 day")
);
# join the two to get a table with the date column on it for plotting
tidy_proj_for_viz <- dplyr::left_join(tidy_proj_for_viz, lut, by = "day");

# plot the model fit
model_fit_viz <- covidseir::plot_projection(tidy_proj, obs_dat=dat, value_column="value", date_column="date");
# plot the model resuduals
residuals_viz <- covidseir::plot_residuals(tidy_proj, obs_dat=dat, obj=fit);
# plot the randomized quantile residuals
set.seed(1);
quantile_residuals_viz <- covidseir::plot_residuals(tidy_proj, obs_dat = dat, obj = fit, type = "quantile");
# plot a histogram of the randimised residuals
resid <- covidseir::plot_residuals(tidy_proj, obs_dat = dat, obj = fit, type = "quantile", return_residuals = TRUE);
# Q-Q plot
qqnorm(resid)
qqline(resid)

# ######################################################################
# PROJECTION
# ######################################################################

# the number of days to project into the future
days_project <- 45;
# the day on which social distancing is relaxed
day_start_reduction <- 5;
# generate a new prediction into the future
proj2 <- covidseir::project_seir(
	obj = fit, # the previous model fit based on the data we have
	forecast_days = days_project, # number of days into the future we're willing to go
	f_fixed_start = max(fit$days) + day_start_reduction, # the day at which social distancing restrictions are loosened
	# f_fixed = NULL,
	f_multi = rep(0.67, days_project - day_start_reduction + 1), # the multiplicative factor by which the social distancing is reduced
	f_multi_seg = 3, # which f segment to use
	iter = 1:50, # using the first 50 posterior samples of the fit (defaults to using all teh values)
	# return_states = FALSE # returning the states form the differential modelling
	# imported_cases = 0, # number of imported cases
	# imported_window = 1 # number of days over which to distribute the imported cases
	# parallel = TRUE # run the computations in parallel rather than serial
	# X = obj$stan_data$X # optional model matrix acting additively on log expected cases
);
# tidy the results
tidy_proj2 <- covidseir::tidy_seir(proj2, resample_y_rep = 30);
# add the date column back on
tidy_proj2 <- dplyr::left_join(tidy_proj2, lut, by = "day");

##############################################################################################################

# Calculating the threshold for increase
# this can be uses to plot the histogram
threshold <- covidseir::get_threshold(
	obj = fit, # the fit that you worked before
	iter = 1:30, # vector of MCMC iterations
	# forecast_days = 25, # number of dauys to forecast into the future
	fs = seq(0.3, 0.8, length.out=4) # contact rate fractions to test to find the threshold
	# show_plot = TRUE # show the scatter plot
	# window_check = 25 # the number of days to use fromt he lst day forecasted
);

##############################################################################################################

#> Finding the MAP estimate.
print(Sys.time() - start_time)

##################################################################################

# get the scejario that's currently running
# myScenario <- scenario();

# # Load scenario's input datasheet from SyncroSim library into R data table
# DTin <- data.table(datasheet(myScenario, name="helloworld_InputDatasheet"));
#
# # Extract model inputs from this R dataframe and then do calculations
# x = DTin[, x];
# a = DTin[, a];
# y = x * a;
#
# # Setup an empty R dataframe ready to accept output in SyncroSim datasheet format
# DTout = datasheet(myScenario, name="helloworld_OutputDatasheet")
#
# # Copy output into this R dataframe
# DTout[1:length(y),"y"] <- y;
#
# # Save this R dataframe back to the SyncroSim library's output datasheet
# saveDatasheet(myScenario, data=DTout, name="helloworld_OutputDatasheet");











###########################################################################################

# Sampling with the NUTS HMC sampler.                                                                                                                                                                                                                                                                                       SAMPLING FOR MODEL 'seir' NOW (CHAIN 1).                                                                                                                     Chain 1:                                                                                                                                                     Chain 1: Gradient evaluation took 0.078125 seconds                                                                                                           Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 781.25 seconds.                                                                  Chain 1: Adjust your expectations accordingly!                                                                                                               Chain 1:                                                                                                                                                     Chain 1:                                                                                                                                                     Chain 1: Iteration:   1 / 500 [  0%]  (Warmup)                                                                                                               Chain 1: Iteration:  50 / 500 [ 10%]  (Warmup)                                                                                                               Chain 1: Iteration: 100 / 500 [ 20%]  (Warmup)                                                                                                               Chain 1: Iteration: 150 / 500 [ 30%]  (Warmup)                                                                                                               Chain 1: Iteration: 200 / 500 [ 40%]  (Warmup)                                                                                                               Chain 1: Iteration: 250 / 500 [ 50%]  (Warmup)                                                                                                               Chain 1: Iteration: 251 / 500 [ 50%]  (Sampling)                                                                                                             Chain 1: Iteration: 300 / 500 [ 60%]  (Sampling)                                                                                                             Chain 1: Iteration: 350 / 500 [ 70%]  (Sampling)                                                                                                             Chain 1: Iteration: 400 / 500 [ 80%]  (Sampling)                                                                                                             Chain 1: Iteration: 450 / 500 [ 90%]  (Sampling)                                                                                                             Chain 1: Iteration: 500 / 500 [100%]  (Sampling)                                                                                                             Chain 1:                                                                                                                                                     Chain 1:  Elapsed Time: 5175.27 seconds (Warm-up)                                                                                                            Chain 1:                5594.44 seconds (Sampling)                                                                                                           Chain 1:                10769.7 seconds (Total)                                                                                                              Chain 1:                                                                                                                                                                                                                                                                                                                  SAMPLING FOR MODEL 'seir' NOW (CHAIN 2).                                                                                                                     Chain 2:                                                                                                                                                     Chain 2: Gradient evaluation took 0.046875 seconds                                                                                                           Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 468.75 seconds.                                                                  Chain 2: Adjust your expectations accordingly!                                                                                                               Chain 2:                                                                                                                                                     Chain 2:                                                                                                                                                     Chain 2: Iteration:   1 / 500 [  0%]  (Warmup)                                                                                                               Chain 2: Iteration:  50 / 500 [ 10%]  (Warmup)                                                                                                               Chain 2: Iteration: 100 / 500 [ 20%]  (Warmup)                                                                                                               Chain 2: Iteration: 150 / 500 [ 30%]  (Warmup)                                                                                                               Chain 2: Iteration: 200 / 500 [ 40%]  (Warmup)                                                                                                               Chain 2: Iteration: 250 / 500 [ 50%]  (Warmup)                                                                                                               Chain 2: Iteration: 251 / 500 [ 50%]  (Sampling)                                                                                                             Chain 2: Iteration: 300 / 500 [ 60%]  (Sampling)                                                                                                             Chain 2: Iteration: 350 / 500 [ 70%]  (Sampling)                                                                                                             Chain 2: Iteration: 400 / 500 [ 80%]  (Sampling)                                                                                                             Chain 2: Iteration: 450 / 500 [ 90%]  (Sampling)                                                                                                             Chain 2: Iteration: 500 / 500 [100%]  (Sampling)                                                                                                             Chain 2:                                                                                                                                                     Chain 2:  Elapsed Time: 4620.41 seconds (Warm-up)                                                                                                            Chain 2:                6178.02 seconds (Sampling)                                                                                                           Chain 2:                10798.4 seconds (Total)                                                                                                              Chain 2:                                                                                                                                                                                                                                                                                                                  SAMPLING FOR MODEL 'seir' NOW (CHAIN 3).                                                                                                                     Chain 3:                                                                                                                                                     Chain 3: Gradient evaluation took 0.046875 seconds                                                                                                           Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 468.75 seconds.                                                                  Chain 3: Adjust your expectations accordingly!                                                                                                               Chain 3:                                                                                                                                                     Chain 3:                                                                                                                                                     Chain 3: Iteration:   1 / 500 [  0%]  (Warmup)                                                                                                               Chain 3: Iteration:  50 / 500 [ 10%]  (Warmup)                                                                                                               Chain 3: Iteration: 100 / 500 [ 20%]  (Warmup)                                                                                                               Chain 3: Iteration: 150 / 500 [ 30%]  (Warmup)                                                                                                               Chain 3: Iteration: 200 / 500 [ 40%]  (Warmup)                                                                                                               Chain 3: Iteration: 250 / 500 [ 50%]  (Warmup)                                                                                                               Chain 3: Iteration: 251 / 500 [ 50%]  (Sampling)                                                                                                             Chain 3: Iteration: 300 / 500 [ 60%]  (Sampling)                                                                                                             Chain 3: Iteration: 350 / 500 [ 70%]  (Sampling)                                                                                                             Chain 3: Iteration: 400 / 500 [ 80%]  (Sampling)                                                                                                             Chain 3: Iteration: 450 / 500 [ 90%]  (Sampling)                                                                                                             Chain 3: Iteration: 500 / 500 [100%]  (Sampling)                                                                                                             Chain 3:                                                                                                                                                     Chain 3:  Elapsed Time: 4911.16 seconds (Warm-up)                                                                                                            Chain 3:                6522.81 seconds (Sampling)                                                                                                           Chain 3:                11434 seconds (Total)                                                                                                                Chain 3:                                                                                                                                                                                                                                                                                                                  SAMPLING FOR MODEL 'seir' NOW (CHAIN 4).                                                                                                                     Chain 4:                                                                                                                                                     Chain 4: Gradient evaluation took 0.046875 seconds                                                                                                           Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 468.75 seconds.                                                                  Chain 4: Adjust your expectations accordingly!                                                                                                               Chain 4:                                                                                                                                                     Chain 4:                                                                                                                                                     Chain 4: Iteration:   1 / 500 [  0%]  (Warmup)                                                                                                               Chain 4: Iteration:  50 / 500 [ 10%]  (Warmup)                                                                                                               Chain 4: Iteration: 100 / 500 [ 20%]  (Warmup)                                                                                                               Chain 4: Iteration: 150 / 500 [ 30%]  (Warmup)                                                                                                               Chain 4: Iteration: 200 / 500 [ 40%]  (Warmup)                                                                                                               Chain 4: Iteration: 250 / 500 [ 50%]  (Warmup)                                                                                                               Chain 4: Iteration: 251 / 500 [ 50%]  (Sampling)                                                                                                             Chain 4: Iteration: 300 / 500 [ 60%]  (Sampling)                                                                                                             Chain 4: Iteration: 350 / 500 [ 70%]  (Sampling)                                                                                                             Chain 4: Iteration: 400 / 500 [ 80%]  (Sampling)                                                                                                             Chain 4: Iteration: 450 / 500 [ 90%]  (Sampling)                                                                                                             Chain 4: Iteration: 500 / 500 [100%]  (Sampling)                                                                                                             Chain 4:                                                                                                                                                     Chain 4:  Elapsed Time: 4989.12 seconds (Warm-up)                                                                                                            Chain 4:                5431 seconds (Sampling)                                                                                                              Chain 4:                10420.1 seconds (Total)                                                                                                              Chain 4:                                                                                                                                                     Time difference of 12.06554 hours                                                                                                                            Error in eval(ei, envir) : object 'dat' not found                                                                                                            In addition: There were 50 or more warnings (use warnings() to see the first 50)                                                                             > delay_data                                                                                                                                                        reported_date symptom_onset_date                                                                                                                          1:    2020-04-28\         2020-04-25                                                                                                                         2:    2020-04-28         2020-04-24                                                                                                                          3:    2020-04-26         2020-04-05                                                                                                                          4:    2020-05-06         2020-05-02                                                                                                                          5:    2020-05-10         2020-05-09                                                                                                                         ---                                                                                                                                                       45022:    2021-01-05         2021-01-01                                                                                                                      45023:    2021-01-05         2021-01-03                                                                                                                      45024:    2021-01-05         2021-01-01                                                                                                                      45025:    2021-01-05         2021-01-02                                                                                                                      45026:    2021-01-05         2021-01-03                                                                                                                      > delay_data <- data.table(get(load("delay_data.rda")))                                                                                                      > write.csv(delay_data, "Delay_Data.csv")                                                                                                                    > write.csv(delay_data, "Delay_Data.csv", row.names=FALSE)                                                                                                   > delay_data                                                                                                                                                        reported_date symptom_onset_date time_to_report                                                                                                           1:    2020-03-15         2020-03-07         8 days                                                                                                           2:    2020-04-28         2020-04-25         3 days                                                                                                           3:    2020-04-28         2020-04-24         4 days                                                                                                           4:    2020-04-26         2020-04-05        21 days                                                                                                           5:    2020-03-24         2020-03-13        11 days                                                                                                          ---                                                                                                                                                       45508:    2021-01-05         2021-01-01         4 days                                                                                                       45509:    2021-01-05         2021-01-03         2 days                                                                                                       45510:    2021-01-05         2021-01-01         4 days                                                                                                       45511:    2021-01-05         2021-01-02         3 days                                                                                                       45512:    2021-01-05         2021-01-03         2 days                                                                                                       > write.csv(delay_data, "Delay_Data.csv", row.names=FALSE)                                                                                                   > delay_data <- data.rable(fread("Delay_Data.csv"))                                                                                                          Error in data.rable(fread("Delay_Data.csv")) :                                                                                                                 could not find function "data.rable"                                                                                                                       > delay_data <- data.table(fread("Delay_Data.csv"))                                                                                                          > q()                                                                                                                                                        Save workspace image? [y/n/c]: c                                                                                                                             >                                                                                                                                                            >                                                                                                                                                            > delay_dat                                                                                                                                                  Error: object 'delay_dat' not found                                                                                                                          > delay_data                                                                                                                                                        reported_date symptom_onset_date time_to_report                                                                                                           1:    2020-03-15         2020-03-07              8                                                                                                           2:    2020-04-28         2020-04-25              3                                                                                                           3:    2020-04-28         2020-04-24              4                                                                                                           4:    2020-04-26         2020-04-05             21                                                                                                           5:    2020-03-24         2020-03-13             11                                                                                                          ---                                                                                                                                                       45508:    2021-01-05         2021-01-01              4                                                                                                       45509:    2021-01-05         2021-01-03              2                                                                                                       45510:    2021-01-05         2021-01-01              4                                                                                                       45511:    2021-01-05         2021-01-02              3                                                                                                       45512:    2021-01-05         2021-01-03              2                                                                                                       > str(delay_data)                                                                                                                                            Classes ‘data.table’ and 'data.frame':  45512 obs. of  3 variables:                                                                                           $ reported_date     : IDate, format: "2020-03-15" "2020-04-28" "2020-04-28" "2020-04-26" ...                                                                 $ symptom_onset_date: IDate, format: "2020-03-07" "2020-04-25" "2020-04-24" "2020-04-05" ...                                                                 $ time_to_report    : int  8 3 4 21 11 15 4 1 3 1 ...                                                                                                        - attr(*, ".internal.selfref")=<externalptr>                                                                                                                > MLE_res                                                                                                                                                    $minimum                                                                                                                                                     [1] 111722                                                                                                                                                                                                                                                                                                                $estimate                                                                                                                                                    [1] 1.587609 5.775352                                                                                                                                                                                                                                                                                                     $gradient                                                                                                                                                    [1]  0.005270410 -0.001398411                                                                                                                                                                                                                                                                                             $code                                                                                                                                                        [1] 1                                                                                                                                                                                                                                                                                                                     $iterations                                                                                                                                                  [1] 16                                                                                                                                                                                                                                                                                                                    > matrix(1,2)                                                                                                                                                     [,1]                                                                                                                                                    [1,]    1                                                                                                                                                    [2,]    1                                                                                                                                                    > matrix({1},{2})                                                                                                                                                 [,1]                                                                                                                                                    [1,]    1                                                                                                                                                    [2,]    1                                                                                                                                                    > data.frame(1,2)_                                                                                                                                           Error: unexpected input in "data.frame(1,2)_"                                                                                                                > data.frame(1,2)                                                                                                                                              X1 X2                                                                                                                                                      1  1  2                                                                                                                                                      > data.table(fread("C:/Users/User/Documents/samplepackage/Delay_Data.csv"))                                                                                  Error in fread("C:/Users/User/Documents/samplepackage/Delay_Data.csv") :                                                                                       File 'C:/Users/User/Documents/samplepackage/Delay_Data.csv' does not exist or is non-readable. getwd()=='/mnt/c/Users/User/Documents/samplepackage'        > fit                                                                                                                                                        Inference for Stan model: seir.                                                                                                                              4 chains, each with iter=500; warmup=250; thin=1;                                                                                                            post-warmup draws per chain=250, total post-warmup draws=1000.                                                                                                                                                                                                                                                                           mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff Rhat                                                                                    R0             3.63    0.01 0.22  3.20  3.49  3.64  3.78  4.06   375 1.01                                                                                    i0             0.04    0.00 0.03  0.01  0.02  0.03  0.05  0.12   364 1.01                                                                                    e              0.81    0.00 0.05  0.70  0.77  0.81  0.84  0.90   351 1.00                                                                                    f_s[1]         0.38    0.00 0.05  0.26  0.34  0.38  0.41  0.47   299 1.00                                                                                    f_s[2]         0.39    0.00 0.05  0.26  0.36  0.39  0.42  0.48   315 1.00                                                                                    f_s[3]         0.53    0.00 0.04  0.45  0.51  0.54  0.56  0.61   286 1.00                                                                                    start_decline 25.67    0.25 4.86 17.90 22.25 25.09 28.43 37.58   380 1.00                                                                                    end_decline   51.37    0.16 3.12 44.33 49.73 51.56 53.49 57.02   381 1.00                                                                                    phi[1]         5.57    0.02 0.56  4.52  5.19  5.55  5.90  6.81   744 1.00                                                                                                                                                                                                                                                 Samples were drawn using NUTS(diag_e) at Thu Jan 14 14:22:44 2021.                                                                                           For each parameter, n_eff is a crude measure of effective sample size,                                                                                       and Rhat is the potential scale reduction factor on split chains (at                                                                                         convergence, Rhat=1).
