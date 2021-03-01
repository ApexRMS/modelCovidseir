rm(list=ls());

# options(width=Sys.getenv("COLUMNS"));

if(!require("pacman")){ install.packages("pacman", repos="https://utstat.toronto.edu/cran/"); }

###################################################################################################
# LOAD ALL THE REQUIRED LIBRARIES
###################################################################################################

pacman::p_load(tidyr, data.table, ggplot2, dplyr, devtools, pkgbuild, remotes);


if(!require("rightTruncation")){ devtools::install_github("andrew-edwards/rightTruncation"); }
library(rightTruncation);

if(!require("covidseir")){ remotes::install_github("seananderson/covidseir", build_vignettes = TRUE); }
library(covidseir);

# # setwd("C:/Users/User/Documents/samplepackage/");
# setwd("/mnt/c/Users/User/Dropbox/samplepackage/");

ymd <- lubridate::ymd;
# # download the file from github
# # each line of the data.table below represents a single case: when the person's symptoms becan, when their case was reported, and the gap between those two
download.file(
	"https://github.com/andrew-edwards/rightTruncation/raw/master/data/delay_data_2021_01_05.rda",
	destfile="delay_data.rda"
);

writeLines("\tImporting and filtering case data, getting delay distribution parameters...\n");

delay_data <- data.table(get(load("delay_data.rda"))) %>%
	mutate_at(vars(time_to_report), as.integer) %>%
	subset(
		"2020-08-01"< symptom_onset_date &
		reported_date < "2020-12-01"
	);

# here, I'm creating a matrix out of the delay data that will them go into the nlm function
# h_day_0 represents the earliest symptom onset in the data (as we filtered it)
# the row number represents which day since 1st October someone got ill
# the column number represents the days since
h_day_0 <- min(delay_data$symptom_onset_date);
num_rows <- (max(delay_data$symptom_onset_date)-h_day_0)[[1]];
num_cols <- (max(delay_data$reported_date)-h_day_0)[[1]];

H_nr <- matrix(0, nrow=num_rows, ncol=num_cols);
for(i in 1:nrow(delay_data))
{
	onset <- (delay_data[i]$symptom_onset_date - h_day_0)[[1]];
	repor <- (delay_data[i]$reported_date - h_day_0)[[1]];
	H_nr[onset, repor] <- H_nr[onset, repor]+1;
}

# nlm - function from the state package - native R
MLE_res = nlm(
  f = negLL_Weibull_counts_matrix, # this function comes form the rightTruncation package - that's where they don't just use fitdistr from the MASS package
  p = c(2, 7), # estimates to start the fitting
  h_nr = H_nr # the matrix of time to report data
);

# the parameters for the distribution
shape_MLE <- MLE_res$estimate[1];
scale_MLE <- MLE_res$estimate[2];

# I reckon this is a good place to start ... we're at a pretty low number of cases (for the initial consitions of the ODE model), and not too biased based on previous case numbers

writeLines("\tRetrieving case data...\n");

download.file(
	"http://www.bccdc.ca/Health-Info-Site/Documents/BCCDC_COVID19_Regional_Summary_Data.csv",
	destfile = "case_data.csv"
);

Day0 <- "2020-07-01";

BC_cases <- data.table(read.csv("case_data.csv")) %>%
	separate(col=Date, into=c("date"), sep=' ') %>%
	mutate_at(vars(date), as.IDate) %>%
	subset(HA=="All" & Province=="BC") %>%
	subset(date>=Day0) %>%
	select(date, Cases_Reported) %>%
	rename(value=Cases_Reported);
BC_cases <- dplyr::as_tibble(BC_cases);

write.csv(BC_cases, file="specific_case_data.csv", row.names=FALSE);

start_time <- Sys.time();

writeLines("\tPreamble to fitting\n");

# fix this to account for Day0 and specific dates
# samp_frac <- rep(0.37, nrow(BC_cases));
samp_frac <- c(rep(0.14, 13), rep(0.21, 38));
samp_frac <- c(samp_frac, rep(0.37, nrow(BC_cases) - length(samp_frac)));

current_segment <- 1;
beta_means_f_param <- c();
beta_sds_f_param <- c();

f_seg <- c(0, rep(current_segment, nrow(BC_cases) - 1))
{
	current_segment <- current_segment + 1;
	beta_means_f_param <- c(beta_means_f_param, 0.4);
	beta_sds_f_param <- c(beta_sds_f_param, 0.2);
}
day_second_rise <- which(BC_cases$date == ymd("2020-10-01"));
if(length(day_second_rise))
{
	f_seg[seq(day_second_rise, length(f_seg))] <- current_segment;
	current_segment <- current_segment + 1;
	beta_means_f_param <- c(beta_means_f_param, 0.95);
	beta_sds_f_param <- c(beta_sds_f_param, 0.2);
}
day_they_got_smarter <- which(BC_cases$date == ymd("2020-11-15"));
if(length(day_they_got_smarter))
{
	f_seg[seq(day_they_got_smarter, length(f_seg))] <- current_segment;
	current_segment <- current_segment + 1;
	beta_means_f_param <- c(beta_means_f_param, 0.45);
	beta_sds_f_param <- c(beta_sds_f_param, 0.2);
}
day_going_up_again <- which(BC_cases$date == ymd("2020-12-15"));
if(length(day_going_up_again))
{
	f_seg[seq(day_going_up_again, length(f_seg))] <- current_segment;
	current_segment <- current_segment + 1;
	beta_means_f_param <- c(beta_means_f_param, 0.7);
	beta_sds_f_param <- c(beta_sds_f_param, 0.2);
}

writeLines("\tFitting the data...\n");

fit <- covidseir::fit_seir(
	daily_cases = BC_cases$value,
	samp_frac_fixed = samp_frac,
	f_seg = f_seg,
	R0_prior = c(log(2.6), 0.2),
	f_prior = cbind(beta_means_f_param, beta_sds_f_param),
	e_prior = c(0.8, 0.05),
	start_decline_prior = c(log(15), 0.1),
	end_decline_prior = c(log(22), 0.1),
	chains = 4,
	iter = 500,
	N_pop = 5.1e6,
	i0_prior = c(log(10), 1), # number of infected cases 30 days before day 1
	delay_shape = shape_MLE,
	delay_scale = scale_MLE,
	time_increment = 0.1,
	########################
	# save_state_predictions = TRUE,
	fit_type = "optimizing"
);

save(fit, file="RUN_fit.rda");
print(Sys.time() - start_time);

writeLines("\nFitting with covariates...\n")
start_time <- Sys.time();

dow <- data.frame(day_of_week = rep(gl(7, 1), 999)[-c(1:6)]);
dow <- dow[seq_len(nrow(BC_cases)), , drop = FALSE];
X <- model.matrix(~ 0 + day_of_week, dow);

fit_dow <- fit_seir(
	daily_cases = BC_cases$value,
	samp_frac_fixed = samp_frac,
	f_seg = f_seg,
	R0_prior = c(log(2.6), 0.2),
	f_prior =cbind(beta_means_f_param, beta_sds_f_param),
	e_prior = c(0.8, 0.05),
	start_decline_prior = c(log(15), 0.1),
	end_decline_prior = c(log(22), 0.1),
	iter = 500,
	N_pop = 5.1e6,
	i0_prior = c(log(10), 1),
	delay_shape = shape_MLE,
	delay_scale = scale_MLE,
	time_increment = 0.1,
	X = X, # <- the model matrix
	###########################
	fit_type = "optimizing",
);

save(fit_dow, file="RUN_fit_dow.rda");
print(Sys.time() - start_time);

writeLines("\n\tPreparing data for graphs...\n");

future::plan(future::multisession);

proj <- covidseir::project_seir(fit, iter = 1:50);
save(proj, file="RUN_proj.rda");

tidy_proj <- covidseir::tidy_seir(proj, resample_y_rep = 20);
first_day <- min(BC_cases$date, na.rm=TRUE);
last_day <- 300;
lut <- dplyr::tibble(
  day = seq_len(last_day),
  date = seq(first_day, first_day + length(day) - 1, by = "1 day")
);
tidy_proj <- dplyr::left_join(tidy_proj, lut, by = "day");
save(tidy_proj, file="RUN_tidy_proj.rda");

#> Finding the MAP estimate.
proj_cov <- covidseir::project_seir(fit_dow, iter = 1:40);
save(proj, file="RUN_proj_cov.rda");

tidy_proj_cov <- covidseir::tidy_seir(proj_cov);
tidy_proj_cov <- dplyr::left_join(tidy_proj_cov, lut, by = "day");
save(tidy_proj_cov, file="RUN_tidy_proj_cov.rda");

# writeLines("\n\tPlotting graphs...\n");
#
# text_size <- 15;
#
# delay_data_fit <- ggplot(delay_data) +
# 	geom_histogram(aes(time_to_report), binwidth=1, fill="grey", colour="black") +
# 	geom_line(data=the_dist, aes(weibull_points, V2), colour="black", size=1) +
# 	labs(x="Days from symptom onset to reporting", y='Frequency') +
# 	theme_bw() +
# 	theme(
# 		axis.title=element_text(size=text_size),
# 		axis.text=element_text(size=text_size-5),
# 		panel.grid.major = element_blank(),
# 		panel.grid.minor = element_blank()
# 	) +
# 	scale_x_continuous(expand=c(0,0)) +
# 	scale_y_continuous(expand=c(0,0));
# 	ggsave(plot=delay_data_fit, file="covidseir_reporting_delay.png", height=6, width=6, units="in", dpi=600);
#
# pl_time_series <- ggplot(BC_cases, aes(date, value)) +
# 	geom_line(colour="blue") +
# 	theme_bw() +
# 	labs(x="Date", y="Reported Cases") +
# 	xlim(c(Day0, max(BC_cases$date))) +
# 	theme(
# 		axis.title=element_text(size=text_size),
# 		axis.text=element_text(size=text_size-5)
# 	) +
# 	scale_y_continuous(expand=c(0,0)) +
# 	scale_x_date(date_breaks = "1 month");
# 	ggsave(plot=pl_time_series, file="covidseir_case_data.png", height=6, width=22, units="in", dpi=600);
#
# fit <- get(load("RUN_fit.rda"));
# tidy_proj <- get(load("RUN_tidy_proj.rda"));
#
# projection_viz <- covidseir::plot_projection(tidy_proj, obs_dat=BC_cases, value_column="value", date_column="date") +
# 	labs(x="Date", y="Reported Cases") +
# 	theme_bw() +
# 	theme(
# 		axis.title=element_text(size=text_size),
# 		axis.text=element_text(size=text_size-5),
# 		panel.grid.major = element_line(size=0),
# 		panel.grid.minor = element_line(size=0)
# 	) +
# 	scale_y_continuous(expand=c(0,0)) +
# 	scale_x_date(date_breaks = "1 month");
# 	ggsave(projection_viz, file="covidseir_data_fit.png", height=6, width=22, dpi=600);
#
# projection_viz <- covidseir::plot_projection(tidy_proj, obs_dat=BC_cases, value_column="value", date_column="date") +
# 	labs(x="Date", y="Reported Cases") +
# 	theme_bw() +
# 	theme(
# 		axis.title=element_text(size=text_size),
# 		axis.text=element_text(size=text_size-5),
# 		panel.grid.major = element_line(size=0),
# 		panel.grid.minor = element_line(size=0)
# 	) +
# 	scale_y_continuous(expand=c(0,0)) +
# 	scale_x_date(date_breaks = "1 month");
# 	ggsave(projection_viz, file="covidseir_data_fit.png", height=6, width=22, dpi=600);
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
# 	ggsave(projection_viz_cov, file="covidseir_data_fit_covariates.png", height=6, width=22, dpi=600);
#
# residuals_viz <- covidseir::plot_residuals(tidy_proj, obs_dat=BC_cases, obj=fit, date_column="date") +
# 	theme_bw() +
# 	labs(x="Date", y="Residual") +
# 	theme(
# 		axis.title=element_text(size=text_size),
# 		axis.text=element_text(size=text_size-5)
# 	) +
# 	scale_y_continuous(expand=c(0,0)) +
# 	scale_x_date(date_breaks = "1 month");
# 	ggsave(residuals_viz, file="covidseir_residual.png", height=6, width=9, dpi=600);
#
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
# 	ggsave(quantile_residuals_viz, file="covidseir_quantile_residual.png", height=6, width=9, dpi=600);
#
# set.seed(1)
# raw_residuals <- plot_residuals(tidy_proj, obs_dat=BC_cases, obj=fit, type="quantile", return_residuals=TRUE, date_column="date", value_column="value");
#
# pl_residuals <- ggplot(data.table(val=raw_residuals), aes(val), fill="grey", colour="black") +
# 	geom_histogram(binwidth=0.5) +
# 	theme_bw() +
# 	labs(x="Residual", y="Frequency") +
# 	theme(
# 		axis.title=element_text(size=text_size),
# 		axis.text=element_text(size=text_size-5),
# 		panel.grid.major = element_line(size=0),
# 		panel.grid.minor = element_line(size=0)
#
# 	);
# 	ggsave(pl_residuals, file="covidseir_quantile_residual.png", height=6, width=9, dpi=600);
#
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
# 	ggsave(pl_normal_qq, file="covidseir_normal_qq_plot.png", height=4, width=4, dpi=600);
#
# MEOECC <- purrr::map_dfr(1:7, ~tibble(dow=.x, b=fit_dow$post[[paste0("beta[", .x, "]")]]));
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
# 	ggsave(plot=pl_MEOECC, file="covidseir_violins.png", height=6, width=12, units="in", dpi=600);



writeLines("\tFinished\n");
