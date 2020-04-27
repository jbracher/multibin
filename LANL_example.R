# Set working directory to the folder containing LANL_example.R
setwd("/home/johannes/Documents/CDClogScore/multibin")

# This script requires the repository https://github.com/FluSightNetwork/cdc-flusight-ensemble
# (at commit 1120 things work) to be downloaded to the same directory as the multibin repository.

# get some functions
source("functions_mbls.R")

# Load targets
targets <- read.csv("../cdc-flusight-ensemble/scores/target-multivals-20172018.csv")
targets <- subset(targets, Location == "US National")

tgs <- list(wk1 = "1 wk ahead", wk2 = "2 wk ahead", wk3 = "3 wk ahead",
            wk4 = "4 wk ahead", ot = "Season onset",
            pt = "Season peak week", pi = "Season peak percentage")

# create an additional folder containing just the LANL forecasts for 16/17
dir.create("Forecasts_LANL_1617")
setwd("Forecasts_LANL_1617")
files_to_copy <- (c(paste0("EW", 40:52, "-2016-LANL_DBMplus.csv"),
                   paste0("EW0", 1:9, "-2017-LANL_DBMplus.csv"),
                   paste0("EW", 10:20, "-2017-LANL_DBMplus.csv")))
file.copy(to = files_to_copy,
          from = paste0("../../cdc-flusight-ensemble/model-forecasts/component-models/LANL_DBMplus/",
                        files_to_copy),
          overwrite = TRUE)
setwd("..")

# get the forecasts from the different weeks (get_forecasts is a custom function
# accessing the different files and shaping things nicely)
F_1617 <- lapply(tgs,
                 get_forecasts,
                 path = "Forecasts_LANL_1617",
                 season = "2016/2017")
# pops one warning, can be ignored

# note: we simplifyingly treat the NA category in onset timing as if it was a regular one

# obtain optimized forecasts:
G_1617 <- list()

for(tg in names(F_1617)){
  print(paste("Starting", tg))

  d_temp <- ifelse(tg %in% c("ot", "pt"), 1, 5)

  G_1617[[tg]] <- optimize_forecasts(original_forecast = F_1617[[tg]], d = d_temp)
}

# warnigns come from cases where all mass is on one bin, this upsets the optimization algorithm


# get expected and observed multi-bin scores:
for(tg in names(F_1617)){
  print(paste("Starting", tg))

  d_temp <- ifelse(tg %in% c("ot", "pt"), 1, 5)

  # get expected scores:
  F_1617[[tg]]$e_mbls <- e_mbls_df(pred_df = F_1617[[tg]]$forecast,
                                   true_df = F_1617[[tg]]$forecast,
                                   d = d_temp)
  G_1617[[tg]]$e_mbls <- e_mbls_df(pred_df = G_1617[[tg]]$forecast,
                                   true_df = F_1617[[tg]]$forecast,
                                   d = d_temp)
  # get observed scores:
  observed_temp <- subset(targets,
                          Season == "2016/2017" &
                            Target == tgs[[tg]] &
                            Calendar.Week %in% c(40:52, 1:20))$Valid.Bin_start_incl

  F_1617[[tg]]$observed <- observed_temp

  F_1617[[tg]]$mbls <- evaluate_forecasts_mbls(forecast = F_1617[[tg]],
                                               observed_bin = observed_temp,
                                               d = d_temp)
  G_1617[[tg]]$mbls <- evaluate_forecasts_mbls(forecast = G_1617[[tg]],
                                               observed_bin = observed_temp,
                                               d = d_temp)
}


# store (needed in manuscript)
save(F_1617, G_1617, file = "forecasts_LANL_DBMplus_1617.rda")

# load("R/LANL_Example/forecasts_LANL_DBMplus_1617.rda")

sapply(F_1617, function(l) mean(l$mbls))
sapply(G_1617, function(l) mean(l$mbls))


# compare to results from CDC table:
dat_scores <- read.csv("../cdc-flusight-ensemble/scores/scores.csv")
table(dat_scores$Model)

scores_1617 <- subset(dat_scores,
                      Location == "US National" &
                        Season == "2016/2017" &
                        Model == "LANL-DBMplus")
aggregate(scores_1617$Multi.bin.score,
          by = list(scores_1617$Target),
          FUN = mean)
sapply(F_1617, function(l) mean(l$mbls))
sapply(G_1617, function(l) mean(l$mbls))

# compare visually:
plot(subset(scores_1617, Target == "1 wk ahead")$Multi.bin.score[c(21:33, 1:20)])
lines(F_1617$wk1$mbls)
lines(G_1617$wk1$mbls, col = "red")

plot(subset(scores_1617, Target == "2 wk ahead")$Multi.bin.score[c(21:33, 1:20)])
lines(F_1617$wk2$mbls)
lines(G_1617$wk2$mbls, col = "red")

plot(subset(scores_1617, Target == "3 wk ahead")$Multi.bin.score[c(21:33, 1:20)])
lines(F_1617$wk3$mbls)
lines(G_1617$wk3$mbls, col = "red")

plot(subset(scores_1617, Target == "4 wk ahead")$Multi.bin.score[c(21:33, 1:20)])
lines(F_1617$wk4$mbls)
lines(G_1617$wk4$mbls, col = "red")


plot(subset(scores_1617, Target == "Season onset")$Multi.bin.score[c(21:33, 1:20)])
lines(F_1617$ot$mbls)
lines(G_1617$ot$mbls, col = "red")

plot(subset(scores_1617, Target == "Season peak week")$Multi.bin.score[c(21:33, 1:20)])
lines(F_1617$pt$mbls)
lines(G_1617$pt$mbls, col = "red")

plot(subset(scores_1617, Target == "Season peak percentage")$Multi.bin.score[c(21:33, 1:20)])
lines(F_1617$pi$mbls)
lines(G_1617$pi$mbls, col = "red")


# restrict to evaluation periods from the original paper:
target_bounds <- read.csv("../cdc-flusight-ensemble/writing/comparison/data/all-target-bounds.csv")

subset(target_bounds, Location == "US National" & Season == "2016/2017")

inds_evaluated <- list(
  wk1 = paste0("iw", c(46:52, 1:17)),
  wk2 = paste0("iw", c(46:52, 1:17)),
  wk3 = paste0("iw", c(46:52, 1:17)),
  wk4 = paste0("iw", c(46:52, 1:17)),
  ot = paste0("iw", c(43:52, 1:4)),
  pt = paste0("iw", c(43:52, 1:17)),
  pi = paste0("iw", c(43:52, 1:17))
)

# results included into arxiv paper:
lapply(c("wk1", "wk2", "wk3", "wk4", "ot", "pt", "pi"),
       function(target){
         mean(F_1617[[target]]$mbls[inds_evaluated[[target]]])
       })

lapply(c("wk1", "wk2", "wk3", "wk4", "ot", "pt", "pi"),
       function(target){
         mean(G_1617[[target]]$mbls[inds_evaluated[[target]]])
       })
