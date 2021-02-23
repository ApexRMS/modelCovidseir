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
);

ymd <- lubridate::ymd;

if(!require("rightTruncation")){ devtools::install_github("andrew-edwards/rightTruncation"); }
library(rightTruncation);

if(!require("covidseir")){ remotes::install_github("seananderson/covidseir", build_vignettes = TRUE); }
library(covidseir);

writeLines("\n\tGathering data...\n");

read.csv("specific_case_data.csv") %>% mutate_at(vars(date), as.IDate) %>% dplyr::as_tibble() -> BC_cases;


# I'm just using this to call the correct files output by the "fit_the_case_data" file
# since the optimising function runs fast (and we use it last file), let's read that set of covidseir objects into this file

TYPE_OF_FIT <- "optimizing"; # takes around 20 seconds
# TYPE_OF_FIT <- "NUTS"; # takes around 23 seconds

# text size for the output graphs
text_size <- 15;

# loading all the covidseir objects produced by the "fit_the_case_data" file
fit <- get(load(sprintf("RUN_%s_fit.rda", TYPE_OF_FIT)));
proj <- get(load(sprintf("RUN_%s_proj.rda", TYPE_OF_FIT)));
tidy_proj <- get(load(sprintf("RUN_%s_tidy_proj.rda", TYPE_OF_FIT)));

# setting up the date lut as in the last file
first_day <- min(BC_cases$date, na.rm=TRUE);
last_day <- 300;
lut <- dplyr::tibble(
  day = seq_len(last_day),
  date = seq(first_day, first_day+length(day)-1, by = "1 day")
);


# parallel processing
future::plan(future::multisession);


# number of days to project into the future
days_project <- 45;

writeLines(sprintf("\n\tMaking a %i day projection...\n", days_project));

# make and plot the prediction. under the assumption that f (the contribution of the
#   distanced category) will change in the future. Here, we're modelling the effects of
#   loosening restrictions during the predicted future
#
#   f_multi - vector of factors multiplying f (instead of setting the actual values) over time of the prediction

day_start_change <- 1;
proj2 <- covidseir::project_seir(
	fit, # the fit object loaded
	iter = 1:50, # MCMC iterations
	forecast_days = days_project, # number of days into the future to predict
	f_fixed_start = max(fit$days) + day_start_change, # the day at which to start changing f
	f_multi = rep(1.1, days_project - day_start_change + 1),
	f_multi_seg = 1, # which f segment to use
	parallel = TRUE # process in parallel
);

# smooth the prediction by resampling
tidy_proj2 <- covidseir::tidy_seir(proj2, resample_y_rep = 30);
tidy_proj2 <- dplyr::left_join(tidy_proj2, lut, by = "day");

save(proj2, file=sprintf("RUN_%s_proj2.rda", TYPE_OF_FIT));
save(tidy_proj2, file=sprintf("RUN_%s_tidy_proj2.rda", TYPE_OF_FIT));

# plot the prediction - the file name will give the name of the stan fit input and the projection length
projection_viz2 <- covidseir::plot_projection(tidy_proj2, obs_dat=BC_cases, value_column="value", date_column="date") +
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

ggsave(
	projection_viz2,
	file=sprintf("covidseir_%s_%i_day_projection.png", TYPE_OF_FIT, days_project),
	height=6, width=22, dpi=600
);


writeLines("\n\tGetting threshold values...\n");

# threshold scatter plot
png("covidseir_threshold_contact_rates.png", pointsize=text_size);
threshold <- covidseir::get_threshold(fit, iter=1:30, fs = seq(0.3, 0.8, length.out=5));
dev.off();

# threshold histogram
pl_histogram <- ggplot(data.table(val=threshold), aes(val))  +
  geom_histogram(binwidth=0.0125, fill="grey", colour="black") +
  theme_bw() +
  labs(x="Threshold", y="Frequency") +
  theme(
    axis.title=element_text(size=text_size),
    axis.text=element_text(size=text_size-5),
    panel.grid.major = element_line(size=0),
    panel.grid.minor = element_line(size=0)
  );
ggsave(pl_histogram, file=sprintf("covidseir_%s_threshold_histogram.png", TYPE_OF_FIT), height=6, width=9, dpi=600);

# threshold violin plot - ratios between posterior f divided by corresponding threshold values
f_ratios <- reshape2::melt(fit$post$f_s[seq_along(threshold), ]/threshold);

the_f_segs <- sort(as.character(unique(f_ratios$Var2)));
the_f_seg_numbers <- gsub("[^[:digit:].]","", as.character(unique(f_ratios$Var2)));

f_ratios$Var2 <- factor(
  f_ratios$Var2,
  levels=the_f_segs,
  labels=sprintf("f seg %s", the_f_seg_numbers)
);

pl_post_thresh_f_ratio_violin <- ggplot(f_ratios, aes(Var2, value)) +
  geom_violin(colour="black", fill="blue") +
  geom_hline(yintercept = 1, lty = 2) +
  labs(x="f segment", y=expression(atop("Contact threshold ratio", "(Posterior/Threshold)"))) +
  theme_bw() +
  theme(
    axis.title=element_text(size=text_size),
    axis.text=element_text(size=text_size-5),
  );
  ggsave(pl_post_thresh_f_ratio_violin, file=sprintf("covidseir_%s_threshold_ratio_violins.png", TYPE_OF_FIT), height=6, width=9, dpi=600);


writeLines("\n\tCalculating effective reproduction number...\n")

rt <- covidseir::get_rt(fit, iter = 1:50);

pl_Rt_viz <- ggplot(rt, aes(time, Rt, group = .iteration)) +
	geom_line(alpha = 0.2, na.rm = TRUE) +
	geom_hline(yintercept = 1, lty = 2) +
	labs(x="Time", y=expression("Effective Reproduction Number (R"['t']*')')) +
	theme_bw() +
	theme(
		axis.title=element_text(size=text_size),
		axis.text=element_text(size=text_size-5),
	) +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0));
	ggsave(plot=pl_Rt_viz, file="covidseir_Rt.png", height=6, width=6, units="in", dpi=600);

writeLines("\tDoubling time...\n");

# doubling time plot
png("covidseir_doubling_times.png", pointsize=text_size);
doubling_times <- covidseir::get_doubling_time(
  fit,
  iter = 1:50,
  # forecast_days = days_project
);
dev.off()

png("covidseir_doubling_times_histogram.png", pointsize=text_size);
hist(doubling_times);
dev.off();

# # doubling time histograms
# pl_double_time <- ggplot(data.table(val=doubling_times), aes(val)) +
#   geom_histogram(binwidth=0.5, fill="grey", colour="black") +
# 	labs(x="Time", y="log(Prevalence)") +
# 	theme_bw() +
# 	theme(
# 		axis.title=element_text(size=text_size),
# 		axis.text=element_text(size=text_size-5),
# 	) +
# 	scale_x_continuous(expand=c(0,0)) +
# 	scale_y_continuous(expand=c(0,0));
# 	ggsave(plot=pl_double_time, file="covidseir_Rt.png", height=6, width=6, units="in", dpi=600);
