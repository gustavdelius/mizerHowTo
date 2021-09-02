# Script used to add observed catch and observed biomass to species params
# data frame `nsParams`
library(mizerHowTo)
time_averaged_catches <- readr::read_csv("data-raw/time-averaged-catches.csv")
time_averaged_SSB <- readr::read_csv("data-raw/time-averaged-SSB.csv")

library(dplyr)
nsParams <- nsParams %>%
  left_join(time_averaged_catches, by = "species") %>%
  left_join(time_averaged_SSB, by = "species") %>%
  rename(yield_observed = Catch_1419_tonnes,
         biomass_observed = SSB_1419)
nsParams$cutoff_size <- sp$w_mat

usethis::use_data(nsParams, overwrite = TRUE)
