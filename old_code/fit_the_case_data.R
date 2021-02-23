rm(list=ls());

# the only solution that worked for me:
# "tinker"'s answer
# https://stackoverflow.com/questions/49196697/how-to-get-the-directory-of-the-executing-script-in-r

library(base)
library(rstudioapi)

get_directory <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file <- "--file="
  rstudio <- "RStudio"

  match <- grep(rstudio, args)
  if (length(match) > 0) {
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  } else {
    match <- grep(file, args)
    if (length(match) > 0) {
      return(dirname(normalizePath(sub(file, "", args[match]))))
    } else {
      return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
    }
  }
}
path <- paste(get_directory(), "/output/", "make_projections.R", sep = "");
file_directory <- paste(as.character(head(strsplit(path, '/')[[1]], -2)), collapse='/');
setwd(file_directory);


if(!require("pacman")){ install.packages("pacman", repos="https://utstat.toronto.edu/cran/"); }

###################################################################################################
# LOAD ALL THE REQUIRED LIBRARIES
###################################################################################################

pacman::p_load(
	tidyr,
	data.table,
	ggplot2,
	dplyr,
	devtools,
	pkgbuild,
	remotes
	# tools/
);

ymd <- lubridate::ymd;

if(!require("rightTruncation")){ devtools::install_github("andrew-edwards/rightTruncation"); }
library(rightTruncation);

if(!require("covidseir")){ remotes::install_github("seananderson/covidseir", build_vignettes = TRUE); }
library(covidseir);

##############################################################################################
# choices are "NUTS", "VB" or "optimizing"
# use optimizing for fast runs (~1-2 minutes)
# NUTS takes ~6-12 hours for the plain fit, same for fitting with covariates
##############################################################################################
# TYPE_OF_FIT <- "NUTS";
TYPE_OF_FIT <- "optimizing";

##############################################################################################
# whether to generate new predictions, or just use old files stored
# FALSE - use old files
# TRUE - run the fits again
# will give a warning if the file doesn't exist
##############################################################################################
USE_EXISTING_FILES <- FALSE;

##############################################################################################
# the beginning of the time series: the fitting starts at this date
# I chose this date because it corresponds to a low case number after the original April case bump
# we also don't get much out of modelling case data before this point
##############################################################################################
Day0 <- "2020-07-01";

text_size <- 15;

if(!(file.exists("delay_data.rda") & USE_EXISTING_FILES))
{
	writeLines("\n\tDownloading the case data and writing to file...\n");

	# ##################################################################################################
	# One of the parameters of the model is $\mu_r$, representing the expected number of observed cases on
	# 	day $r$. Since there is a delay between symptom onset and reporting of the case, we calculate this
	# 	expected value by an integral, where $w(s)$ is a density function for a delay term $s$. This delay
	# 	density function takes the form of a Weibull distribution with shape and scale parameters we can
	# 	get from reporting delay data. We'll get the data and perform MLE to get these two numbers
	#
	# The data we get form the github page is an .rda giving the time of symptom onset, the time the case
	# 	was reported, and the number of days between those two. It looks like this:
	#
	# delay_data
	# A tibble: 2,355 x 3
	# <date>        <date>             <drtn>
	# 1 2020-06-02    2020-05-31          2 days
	# 2 2020-03-15    2020-03-07          8 days
	# 3 2020-04-28    2020-04-25          3 days
	#
	# Each row represents as single case, so that row 3 represents someone who lost their taste on 28 April
	# 	and had their case reported 3 days later.

    # https://github.com/andrew-edwards/rightTruncation/raw/master/data/delay_data_2021_01_05.rda

    # http://www.bccdc.ca/Health-Info-Site/Documents/BCCDC_COVID19_Regional_Summary_Data.csv

	download.file(
		"https://github.com/andrew-edwards/rightTruncation/raw/master/data/delay_data_2021_01_05.rda",
		destfile="delay_data.rda", quiet=TRUE
	);

	writeLines("\tImporting and filtering case data, getting delay distribution parameters...\n");

	# load the file and get the delay data
	# I change the "time_to_report" column to an integer, getting rid of the " days" part
	# I also take the cases between the 1 October and 1 December
	# 1 October is where case numbers start to rise
	# 1 December is after the mid-November case peak
	# this is also mirrored in their paper (page 4)
	# though it was early days for them, I reckon this is still reasonable
	delay_data <- data.table(get(load("delay_data.rda"))) %>%
	mutate_at(vars(time_to_report), as.integer) %>%
	mutate_at(vars(symptom_onset_date, reported_date), as.Date) %>%
	subset(
		"2020-08-01"< symptom_onset_date &
		reported_date < "2020-12-01"
	);
	save(delay_data, file="delay_data.rda");

	# Below, the actual MLE is done by the function nlm from the R "stats" package. It takes three arguments:
	# 	1) f - the function of the data that is to be minimized
	# 	2) p - initial parameter values (shape, scale) for which to calculate NLL
	# 	3) h_nr - a square matrix argument to be passed to the function f

	# We'll call earliest date present in the table (first date on symptom onset) "h_day_0"; remember that we've
	# 	chosen 1 October 2020 above. We transform the delay_data table to make the h_nr matrix so that the row
	# 	number represents which day since 1 October at which the symptoms started, the column number represents
	# 	which day (since 1 October) that the case was reported, and the value at that position represents the
	# 	number of cases present in the data set with that combination of dates. For example, h_nr[8,10] gives
	# 	the number of people that started experiencing symptoms on 8 October (1 October + 7 days, since R is
	# 	1-indexed) and reported the case on 10 October. The entries of the matrix indicate time as well as gap.
	#	The package provides a function to turn the delay_data tibble into a matrix and another function to work
		# directly with a tibble but the below manual approach works.

	# set the first day so that we can calculate the time that passed in days
	h_day_0 <- min(delay_data$symptom_onset_date);
	# set the number of rows and columns of the matrix
	num_rows <- (max(delay_data$symptom_onset_date)-h_day_0)[[1]];
	num_cols <- (max(delay_data$reported_date)-h_day_0)[[1]];

	# fill the H_nr matrix with the case data as described in the previous long comment
	H_nr <- matrix(0, nrow=num_rows, ncol=num_cols);
	for(i in 1:nrow(delay_data))
	{
		onset <- (delay_data[i]$symptom_onset_date - h_day_0)[[1]];
		repor <- (delay_data[i]$reported_date - h_day_0)[[1]];
		H_nr[onset, repor] <- H_nr[onset, repor]+1;
	}

	# Given that we only have data to the current date, we can think of the data as being truncated on the
	# 	right (treating time as moving left-to-right). This is where the "rightTruncation" R package comes in.
	# 	Below, the actual MLE is done by the function nlm from the R "stats" package. The "rightTruncation"
	# 	package provides the function "negLL_Weibull_counts_matrix", which calculates the log-likelihood of
	# 	the shape and scale parameters of the count data fed to it. p is the vector of values for which the
	# 	NLL is calculated. It usually gives warnings for the NAs produced during estimation, so I've
	# 	suppressed those. I reckon we could just use fitdistr() from the MASS package if the data could be
	# 	considered complete.

	MLE_res = suppressWarnings(nlm(
		f = negLL_Weibull_counts_matrix,
		p = c(2, 7),
		h_nr = H_nr
	));

	# save the output to be read later
	save(MLE_res, file="RUN_MLE_res.rda");
}

##########################################################
# Plotting the distribution of the delay case data and the Weibull distribution with the parameters we get above

delay_data <- data.table(get(load("delay_data.rda")));
MLE_res2 <- get(load("RUN_MLE_res.rda"));

# the parameters for the distribution
shape_MLE <- MLE_res$estimate[1];
scale_MLE <- MLE_res$estimate[2];

# smooth Weibull curve to be superimposed
weibullPoints <- seq(min(delay_data$time_to_report), max(delay_data$time_to_report), by=0.05);
the_dist <- data.table(cbind(weibullPoints, nrow(delay_data)*dweibull(weibullPoints, shape=shape_MLE, scale=scale_MLE)));

delay_data_fit <- ggplot(delay_data) +
	geom_histogram(aes(time_to_report), binwidth=1, fill="grey", colour="black") +
	geom_line(data=the_dist, aes(weibullPoints, V2), colour="black", size=1) +
	labs(x="Days from symptom onset to reporting", y='Frequency') +
	theme_bw() +
	theme(
		axis.title=element_text(size=text_size),
		axis.text=element_text(size=text_size-5),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank()
	) +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0));
	ggsave(plot=delay_data_fit, file="covidseir_reporting_delay.png", height=6, width=6, units="in", dpi=600);

writeLines("\tRetrieving case data...\n");

if(!(file.exists("case_specific_data.rda") & USE_EXISTING_FILES))
{
	# getting the BC regional summary case data csv  from the BC CDC website
	download.file(
		"http://www.bccdc.ca/Health-Info-Site/Documents/BCCDC_COVID19_Regional_Summary_Data.csv",
		destfile = "case_data.csv", quiet=TRUE
	);

	# Filtering the downloaded case data. We'll take the combined data from all health authorities, with all
	# 	cases occurring after Day0, and renaming the Date and Case_Reported columns as "date" and "value"
	# 	respectively.
	BC_cases <- data.table(read.csv("case_data.csv")) %>%
		separate(col=Date, into=c("date"), sep=' ') %>%
		mutate_at(vars(date), as.IDate) %>%
		subset(HA=="All" & Province=="BC") %>%
		subset(date>=Day0) %>%
		select(date, Cases_Reported) %>%
		rename(value=Cases_Reported) %>%
		dplyr::as_tibble();

	write.csv(BC_cases, file="specific_case_data.csv", row.names=FALSE);
}

BC_cases <- data.table(fread("specific_case_data.csv"));

# Creating a lut (tibble) so that we can join a date column on to future output object for plotting
first_day <- min(BC_cases$date, na.rm=TRUE);
last_day <- 300;
lut <- dplyr::tibble(
	day = seq_len(last_day),
	date = seq(first_day, first_day + length(day) - 1, by = "1 day")
);

# Plotting the case data from Day0 until the most recent report in the BC CDC data
pl_time_series <- ggplot(BC_cases, aes(date, value)) +
	geom_line(colour="blue") +
	geom_point() +
	theme_bw() +
	labs(x="Date", y="Reported Cases") +
	ylim(c(0, 1.1*max(BC_cases$value))) +
	theme(
		axis.title=element_text(size=text_size),
		axis.text=element_text(size=text_size-5)
	) +
	scale_x_date(date_breaks = "1 month");
	ggsave(plot=pl_time_series, file="covidseir_case_data.png", height=6, width=22, units="in", dpi=600);

writeLines("\tPreamble to fitting\n");

# estimated fraction of cases sampled
samp_frac <- c(rep(0.14, 13), rep(0.21, 38));
samp_frac <- c(samp_frac, rep(0.37, nrow(BC_cases) - length(samp_frac)));

# $f$ represents the amount contributed by socially distanced individuals to the force of infection. A very bad
# 	analogy/illustrative misinterpretation: $f$ is similar to the proportion of an individual's normal contacts
# 	that is maintained during the simulation. $f=1$ represents no social distancing (you communicate with
# 	100% of your social contacts, and $f=0.2$ represents social distancing to the extent where only 20% of
# 	social contacts are maintained). f=1, no distancing, f=0 complete quarantine (for distanced individuals in
# 	the model). The way it's implemented here is on a time-varying basis; there are segments modulating the
# 	value of the f parameter over time. Say there was one range where protocols were relaxed and a following
# 	period where distancing measures were enforced by some policy change, these two values of f would be
# 	different. Here, as a rough initial fit, I've made a few segments matching large trends in the data; these
# 	segments are numbered consecutively. Naturally, these can be refined to match various types of data, for
# 	instance Google transit data, tracked policy changes, etc.
#
# The segments form a lut; the number $n$ of the segment (segment 1, segment 4, etc) is substituted for the beta
# 	mean and sd of the nth f-segment.

# the starting segment is naturally #1
current_segment <- 1;
# the values of the f parameter in each segment, being estimations, are given as beta mean and sd pairs
beta_means_f_param <- c();
beta_sds_f_param <- c();

f_seg <- c(0, rep(current_segment, nrow(BC_cases) - 1))
{
	# at the start of the simulation, say the distanced individuals contribute f=0.4 to the force of infection
	current_segment <- current_segment + 1;
	beta_means_f_param <- c(beta_means_f_param, 0.4);
	beta_sds_f_param <- c(beta_sds_f_param, 0.2);
}
day_second_rise <- which(BC_cases$date == ymd("2020-10-01"));
if(length(day_second_rise))
{
	# if restrictions were loosened around 1st Oct, f goes up in this segment
	f_seg[seq(day_second_rise, length(f_seg))] <- current_segment;
	current_segment <- current_segment + 1;
	beta_means_f_param <- c(beta_means_f_param, 0.95);
	beta_sds_f_param <- c(beta_sds_f_param, 0.2);
}
day_they_got_smarter <- which(BC_cases$date == ymd("2020-11-15"));
if(length(day_they_got_smarter))
{
	# case numbers rose meteorically, so restrictions were tightened on 15 Nov, bringing f down to 0.45
	f_seg[seq(day_they_got_smarter, length(f_seg))] <- current_segment;
	current_segment <- current_segment + 1;
	beta_means_f_param <- c(beta_means_f_param, 0.45);
	beta_sds_f_param <- c(beta_sds_f_param, 0.2);
}
day_going_up_again <- which(BC_cases$date == ymd("2020-12-15"));
if(length(day_going_up_again))
{	# loosened restrictions on 15 Dec, but not like before - f rises to 0.7
	f_seg[seq(day_going_up_again, length(f_seg))] <- current_segment;
	current_segment <- current_segment + 1;
	beta_means_f_param <- c(beta_means_f_param, 0.7);
	beta_sds_f_param <- c(beta_sds_f_param, 0.2);
}

# names of data sets pulled into the sim
fit_filename <- sprintf("RUN_%s_fit.rda", TYPE_OF_FIT);
proj_filename <- sprintf("RUN_%s_proj.rda", TYPE_OF_FIT);
tidy_proj_filename <- sprintf("RUN_%s_tidy_proj.rda", TYPE_OF_FIT);

if(!(file.exists(fit_filename) & USE_EXISTING_FILES))
{
	writeLines("\tFitting the data...\n");
	start_time <- Sys.time();

	# covidseir::fit_seir fits an SEIR compartmental model to the case data supplied, using the
	# 	distancing parameters, delay distribution parameters, ODE solver parameters and other rstan
	# 	parameters. Naturally, not all parameters are used here; quite a few were left at the
	# 	function's default values.
	#
	# We can also add data from hospitalisations Here

	fit <- covidseir::fit_seir(
		daily_cases = BC_cases$value, # the number of daily cases
		obs_model="NUTS", # observation model
		samp_frac_fixed = samp_frac, # proportion of cases sampled
		f_seg = f_seg, # segments of f values
		R0_prior = c(log(2.6), 0.2), # prior for the basic reproductive rate
		f_prior = cbind(beta_means_f_param, beta_sds_f_param), # beta means and sds of the f segment values
		e_prior = c(0.8, 0.05),	# beta mean and sd for the e parameter representing the proportion of people social distancing
		chains = 4, # number of MCMC chains for rstan
		iter = 500, # number of MCMC iterations for rstan
		N_pop = 5.1e6, # model population size
		i0_prior = c(log(10), 1), # number of infected cases on the first day
		delay_shape = shape_MLE, # shape parameter of the reporting delay distribution
		delay_scale = scale_MLE, # scale parameter for the reporting delay distribution
		time_increment = 0.1, # time increments for the ODE solver and delay distribution model integration
		########################
		# type of fit - argument to rstan - "optimising" runs in a few seconds, and overlaps the result of longer, more
		#	resource-intensive runs, so these short ones are good for now
		fit_type = TYPE_OF_FIT
	);

	save(fit, file=fit_filename);
	print(Sys.time() - start_time);
}

future::plan(future::multisession);

writeLines("\n\tPreparing data for graphs...\n");

# load the fit object from a file
fit <- get(load(fit_filename));

# Anderson's vignette shows him plotting the fit by doing a 0-day projection and saving that object
if(!(file.exists(proj_filename) & USE_EXISTING_FILES))
{
	proj <- covidseir::project_seir(fit, iter = 1:50);
	save(proj, file=proj_filename);
}

proj <- get(load(proj_filename));

# we can smooth the results by resampling
if(!(file.exists(tidy_proj_filename) & USE_EXISTING_FILES))
{
	tidy_proj <- covidseir::tidy_seir(proj, resample_y_rep = 20);
	tidy_proj <- dplyr::left_join(tidy_proj, lut, by = "day");
	save(tidy_proj, file=tidy_proj_filename);
}

tidy_proj <- get(load(tidy_proj_filename));

# plot the fit (NOT a projection here, despite the name)
projection_viz <- covidseir::plot_projection(tidy_proj, obs_dat=BC_cases, value_column="value", date_column="date") +
	labs(x="Date", y="Reported Cases") +
	theme_bw() +
	theme(
		axis.title=element_text(size=text_size),
		axis.text=element_text(size=text_size-5),
		panel.grid.major = element_line(size=0),
		panel.grid.minor = element_line(size=0)
	) +
	scale_y_continuous(expand=c(0,0)) +
	scale_x_date(date_breaks = "1 month");
	ggsave(projection_viz, file=sprintf("covidseir_%s_data_fit.png", TYPE_OF_FIT), height=6, width=22, dpi=600);

# # plot the model residuals
# residuals_viz <- covidseir::plot_residuals(tidy_proj, obs_dat=BC_cases, obj=fit, date_column="date") +
# 	theme_bw() +
# 	labs(x="Date", y="Residual") +
# 	theme(
# 		axis.title=element_text(size=text_size),
# 		axis.text=element_text(size=text_size-5)
# 	) +
# 	scale_y_continuous(expand=c(0,0)) +
# 	scale_x_date(date_breaks = "1 month");
# 	ggsave(residuals_viz, file=sprintf("covidseir_%s_residual.png", TYPE_OF_FIT), height=6, width=9, dpi=600);
#
# # plot quantile residuals
# set.seed(1);
# quantile_residuals_viz <- covidseir::plot_residuals(tidy_proj, obs_dat=BC_cases, obj=fit, type="quantile", date_column="date") +
# 	theme_bw() +
# 	labs(x="Date", y="Quantile Residual") +
# 	theme(
# 		axis.title=element_text(size=text_size),
# 		axis.text=element_text(size=text_size-5),
# 		panel.grid.major = element_line(size=0),
# 		panel.grid.minor = element_line(size=0)
# 	);
# 	ggsave(quantile_residuals_viz, file=sprintf("covidseir_%s_quantile_residual.png", TYPE_OF_FIT), height=6, width=9, dpi=600);
#
# # get the residuals themselves for a histogram
# set.seed(1)
# raw_residuals <- plot_residuals(tidy_proj, obs_dat=BC_cases, obj=fit, type="quantile", return_residuals=TRUE, date_column="date", value_column="value");
#
# # histogram of the residuals
# pl_residuals <- ggplot(data.table(val=raw_residuals), aes(val)) +
# 	geom_histogram(binwidth=0.5, fill="grey", colour="black") +
# 	theme_bw() +
# 	labs(x="Residual", y="Frequency") +
# 	theme(
# 		axis.title=element_text(size=text_size),
# 		axis.text=element_text(size=text_size-5),
# 		panel.grid.major = element_line(size=0),
# 		panel.grid.minor = element_line(size=0)
#
# 	);
# 	ggsave(pl_residuals, file=sprintf("covidseir_%s_quantile_residual.png", TYPE_OF_FIT), height=6, width=9, dpi=600);
#
# # Q-Q plot
# pl_normal_qq <- ggplot(data=as.data.table(qqnorm(raw_residuals)), aes(sample=y)) +
# 	stat_qq() +
# 	stat_qq_line() +
# 	theme_bw() +
# 	labs(x="Theoretical Quantiles", y="Sample Quantiles") +
# 	theme(
# 		axis.title=element_text(size=text_size),
# 		axis.text=element_text(size=text_size-5),
# 		panel.grid.major = element_line(size=0),
# 		panel.grid.minor = element_line(size=0)
# 	);
# 	ggsave(pl_normal_qq, file=sprintf("covidseir_%s_normal_qq_plot.png", TYPE_OF_FIT), height=4, width=4, dpi=600);
#
# ###################################################################################################################
#
# adding covariates, in this case, days of the week to the fit
writeLines("\nFitting with covariates...\n")
start_time <- Sys.time();

fit_cov_filename <- sprintf("RUN_%s_fit_cov.rda", TYPE_OF_FIT);
proj_cov_filename <- sprintf("RUN_%s_proj_cov.rda", TYPE_OF_FIT);
tidy_proj_cov_filename <- sprintf("RUN_%s_tidy_proj_cov.rda", TYPE_OF_FIT);

dow <- data.frame(day_of_week = rep(gl(7, 1), 999)[-c(1:6)]);
dow <- dow[seq_len(nrow(BC_cases)), , drop = FALSE];

if(!(file.exists(fit_cov_filename) & USE_EXISTING_FILES))
{
	X <- model.matrix(~ 0 + day_of_week, dow);

	fit_cov <- fit_seir(
		daily_cases = BC_cases$value,
		samp_frac_fixed = samp_frac,
		f_seg = f_seg,
		R0_prior = c(log(2.6), 0.2),
		f_prior =cbind(beta_means_f_param, beta_sds_f_param),
		e_prior = c(0.8, 0.05),
		iter = 500,
		N_pop = 5.1e6,
		i0_prior = c(log(10), 1),
		delay_shape = shape_MLE,
		delay_scale = scale_MLE,
		time_increment = 0.1,
		X = X, # <- the model matrix for covariate fit
		###########################
		fit_type = TYPE_OF_FIT,
	);

	save(fit_cov, file=fit_cov_filename);
	print(Sys.time() - start_time);
}

# # Rinse and repeat the previous steps, but with the new covidseir objects
#
# fit_cov <- get(load(fit_cov_filename));
#
# if(!(file.exists(proj_cov_filename) & USE_EXISTING_FILES))
# {
# 	proj_cov <- covidseir::project_seir(fit_cov, iter = 1:40);
# 	save(proj_cov, file=proj_cov_filename);
# }
#
# proj_cov <- get(load(proj_cov_filename));
#
# if(!(file.exists(tidy_proj_cov_filename) & USE_EXISTING_FILES))
# {
# 	tidy_proj_cov <- covidseir::tidy_seir(proj_cov);
# 	tidy_proj_cov <- dplyr::left_join(tidy_proj_cov, lut, by = "day");
# 	save(tidy_proj_cov, file=tidy_proj_cov_filename);
# }
#
# tidy_proj_cov <- get(load(tidy_proj_cov_filename));
#
# writeLines("\n\tPlotting graphs with covariates...\n");
#
# projection_viz_cov <- covidseir::plot_projection(tidy_proj_cov, obs_dat=BC_cases, value_column="value", date_column="date") +
# 	labs(x="Date", y="Reported Cases (added covariates)") +
# 	theme_bw() +
# 	theme(
# 		axis.title=element_text(size=text_size),
# 		axis.text=element_text(size=text_size-5),
# 		panel.grid.major = element_line(size=0),
# 		panel.grid.minor = element_line(size=0)
# 	) +
# 	scale_y_continuous(expand=c(0,0)) +
# 	scale_x_date(date_breaks = "1 month");
# 	ggsave(projection_viz_cov, file=sprintf("covidseir_%s_data_fit_covariates.png", TYPE_OF_FIT), height=6, width=22, dpi=600);
#
# MEOECC <- purrr::map_dfr(1:7, ~tibble(dow=.x, b=fit_cov$post[[paste0("beta[", .x, "]")]]));
# pl_MEOECC <- ggplot(MEOECC, aes(dow, exp(b), group = dow)) +
# 	geom_violin(colour="black", fill="blue") +
# 	labs(x="Day of the week", y="Multiplicative effect on expected case counts") +
# 	theme_bw() +
# 	scale_x_continuous(breaks = 1:7, labels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) +
# 	theme(
# 		axis.title=element_text(size=text_size),
# 		axis.text=element_text(size=text_size-5),
# 		panel.grid.major.x = element_line(size=0),
# 		panel.grid.minor.x = element_line(size=0)
# 	);
# 	ggsave(plot=pl_MEOECC, file=sprintf("covidseir_%s_violins.png", TYPE_OF_FIT), height=6, width=12, units="in", dpi=600);
#
# writeLines("\tFinished\n");
