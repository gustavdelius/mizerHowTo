#' Interaction matrix for the North Sea
#'
#' ....
"interNS"

#' ICES fisheries mortality
#'
#' ....
"Fmat"

#' ICES fisheries mortality processed
#' On the ICES database, some species have multiple entries and some have none.
#' Fisheries mortality is weighted by the species abundance to be able to take into account multiple database entries
#' Gurnard data is extrapolated from whiting (same gear and physiology)
#' ....
"FmatWeightedInterpolated"

#' HTM1 final simulation saved (sim_guessed)
#' Used at the beginning of HTM2
#'
#' ....
"HTM1_sim"

#' HTM2 Fmsy plot file (step 2 - Fmsy1 plot_dat)
#' plotFmsy() runs slowly so saving the result for fast kniting
#'
#' ....
"HTM2_Fmsy1"

#' HTM2 Fmsy plot file (step 2 - Fmsy2 plot_dat)
#' plotFmsy() runs slowly so saving the result for fast knitting
#'
#' ....
"HTM2_Fmsy2"

#' HTM2 optim parallel output (step 1 - optim_result)
#' The Rmax vector output of optim parallel is saved for fast knitting
#' ....
"HTM2_optimParallel_Rmax1"

#' HTM2 optim parallel output (step 2 - optim_loop)
#' The Rmax vector output of optim parallel is saved for fast knitting
#' ....
"HTM2_optimParallel_Rmax2"

#' HTM2 saved simulation (step 1 - results)
#' Used as checkpoint/kick start if .rmd are not run from the beginning
#' ....
"HTM2_sim_optim1"

#' HTM2 saved simulation (step 2 - tweaked erepro)
#' Used as checkpoint/kick start if .rmd are not run from the beginning
#' ....
"HTM2_sim_optim2"

#' HTM2 saved simulation (step 2 - sim_loop)
#' Used as checkpoint/kick start if .rmd are not run from the beginning
#' ....
"HTM2_sim_optim3"

#' Species specific parameters for the North Sea
#'
#' ....
"nsParams"

#' averaged catch in the North Sea 2014 - 2019
#'
#' ....
"catchAvg"

