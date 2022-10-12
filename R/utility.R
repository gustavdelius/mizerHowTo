#' Open tutorials
#'
#' Opens the file of a tutorial in RStudio so that the user can read it or
#' play with the code.
#'
#' @param tutorial_name The name of the tutorial. Accepted values are "HTM1"
#' and "HTM2". Default is "HTM1".
#' @param extension The extension of the tutorial. Accepted values are "html"
#' "website" (direct to sizespectrum.org) and "Rmd". Default is "html".
#' @export
tutorial <- function(tutorial_name = "HTM1", extension = "html") {
    switch (tutorial_name,
            "HTM0" = {
                if(extension == "Rmd")
                {
                    rstudioapi::navigateToFile(
                        system.file("HTM0", "HTM0_whyUseMizer.Rmd",
                                    package = "mizerHowTo"))
                } else if (extension == "html") browseURL(system.file("HTM0", "HTM0_whyUseMizer.html",
                                                                      package = "mizerHowTo"))
                #else if (extension == "website") browseURL("https://sizespectrum.org/mizerHowTo/articles/HTM1_parametrisation.html") # is the URL ready?
            },
            "HTM1" = {
                if(extension == "Rmd")
                {
                    rstudioapi::navigateToFile(
                        system.file("HTM1", "HTM1_parametrisation.Rmd",
                                    package = "mizerHowTo"))
                } else if (extension == "html") browseURL(system.file("HTM1", "HTM1_parametrisation.html",
                                                                      package = "mizerHowTo"))
                else if (extension == "website") browseURL("https://sizespectrum.org/mizerHowTo/articles/HTM1_parametrisation.html")
            },
            "HTM2" = {
                if(extension == "Rmd")
                {
                    rstudioapi::navigateToFile(
                        system.file("HTM2", "HTM2_timeAveraged_calibration.Rmd",
                                    package = "mizerHowTo"))
                } else if (extension == "html") browseURL(system.file("HTM2", "HTM2_timeAveraged_calibration.html",
                                                                      package = "mizerHowTo"))
                else if (extension == "website") browseURL("https://sizespectrum.org/mizerHowTo/articles/HTM2_timeAveraged_calibration.html")
            },
            "HTM3" = {
                if(extension == "Rmd")
                {
                    rstudioapi::navigateToFile(
                        system.file("HTM3", "HTM3_timeSeries_calibration.Rmd",
                                    package = "mizerHowTo"))
                } else if (extension == "html") browseURL(system.file("HTM3", "HTM3_timeSeries_calibration.Rmd",
                                                                      package = "mizerHowTo"))
                # else if (extension == "website") browseURL("https://sizespectrum.org/mizerHowTo/articles/HTM2_timeAveraged_calibration.html") # is URL ready?
            },
            {print("Something went wrong")}
    )
}

#' Create ggplot legend object
#'
#' Function that extracts the legend from ggplot object to use it in grid plots
#'
#' @param a.gplot A ggplot object

g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}


#' Compares modelled and empirical output
#'
#' Function that takes as input a vector of Rmax values and empirical data of catch or SSB
#' It runs a mizer simulation with the Rmax values,calculate the difference of predicted
#' and observed yield (or biomass) and return the sum of squared errors of the difference.
#' @param vary Rmax vector, needs to be the same lenght as the number of species in params
#' @param params mizerParams object containing the species parameters
#' @param dat empirical data of catch or SSB
#' @param data_type type of the data given in `dat`. Either "catch" or "SSB". Default is "catch"
#' @param tol Default is 0.1
#' @param timetorun lenght of simulation. Default is 10
#'
#' @export
getError <- function(vary,params,dat,data_type="catch", tol = 0.1,timetorun=10)
{
    #env$params@species_params$R_max[]<-10^vary[1:12]
    params@species_params$R_max[]<-10^vary[1:length(params@species_params$R_max)]

    params <- setParams(params)
    # run to steady state and update params
    # env$params<- projectToSteady(env$params, distance_func = distanceSSLogN,
    #                 tol = tol, t_max = 200,return_sim = F)
    params<- projectToSteady(params, distance_func = distanceSSLogN,
                             tol = tol, t_max = 200,return_sim = F)

    # create sim object

    sim <- project(params, effort = 1, t_max = timetorun, progress_bar = F) #Change t_max to determine how many years the model runs for

    #
    # sim <- project(env$params, effort = 1, t_max = timetorun) #Change t_max to determine how many years the model runs for
    #
    # env$params <-sim@params
    #

    ## what kind of data and output do we have?
    if (data_type=="SSB") {
        output <-getSSB(sim)[timetorun,]   #could change to getBiomass if using survey, also check units.
    }
    if (data_type=="catch") {
        output <-getYield(sim)[timetorun,]/1e6
    }

    pred <- log(output)
    dat  <- log(dat)
    # sum of squared errors, here on log-scale of predictions and data (could change this or use other error or likelihood options)
    discrep <- pred - dat
    discrep <- (sum(discrep^2))

    # can use a strong penalty on the error to ensure we reach a minimum of 10% of the data (biomass or catch) for each species
    # if(any(pred < 0.1*dat)) discrep <- discrep + 1e10

    return(discrep)
}



#' Wrapper for optimParallel
#'
#' The function sets up the parallel environment to run optimParallel
#' on a params object
#'
#' @param params A mizer params object
#' @param vary Dataframe containing which parameters should vary and their upper and lower boundary
#' @param errorFun What error function should be used
#' @param errorArgs Arguments for getError if necessary as a list (could contain params)
#' @param observed_data Data getError is compared to
#'
#' @export

fastOptim <- function(params, vary, errorFun, errorArgs = NULL)
{
    # set up workers
    noCores <- parallel::detectCores() - 1 # keep some spare core
    cl <- parallel::makeCluster(noCores, setup_timeout = 0.5)
    setDefaultCluster(cl = cl)
    clusterExport(cl, varlist = "cl",envir=environment())
    clusterEvalQ(cl, {
        library(mizerExperimental)
        library(optimParallel)
    })

    optim_result <- optimParallel::optimParallel(par = vary$data,
                                                 fn = errorFun,
                                                 params = params,
                                                 errorArgs = errorArgs,
                                                 method   ="L-BFGS-B",
                                                 lower= vary$lower,
                                                 upper= vary$upper,
                                                 parallel=list(loginfo=TRUE, forward=TRUE))#,
    #plankton_forcing, therMizerEncounter, therMizerPredRate, therMizerEReproAndGrowth, scaled_temp_fun)
    stopCluster(cl)

    return(optim_result)
}


#' Improving getError
#'
#'
#'
#'




getErrorCustom <- function(vary, params, errorArgs, tol = 0.001,
                           timetorun = 10)
{
    # params@species_params$R_max[1:9]<-10^vary[1:9]
    # params@species_params$erepro[1:9]<-vary[10:18]
    # params@species_params$interaction_resource[1:9] <- vary[19:27]
    #
    # params <- setParams(params)
    #
    # interaction <- params@interaction
    # interaction[] <- matrix(vary[28:108],nrow = 9) # stop at 54 if looking only at 3 biggest species
    #
    # params <- setInteraction(params,interaction)

    params@species_params$erepro[1] <- vary$data

    params <- projectToSteady(params, distance_func = distanceSSLogN,
                              tol = tol, t_max = length(times), return_sim = F)

    sim <- project(params,t_max = length(times),t_start = 1950, progress_bar = F,initial_n =  params@initial_n, initial_n_pp = params@other_params$other$n_pp_array[1,])

    # get biomass through time
    biomass <- sweep(sim@n, 3, sim@params@w * sim@params@dw, "*")
    biomass <- biomass[which(dimnames(biomass)$time == "2003"):which(dimnames(biomass)$time == "2020"),,] # trageting fishing period

    #get yield through time from model:

    f_gear<-mizer::getFMortGear(params,effort)
    # f_gear <- f_gear[1:13,,,,drop=F]
    yield_species_gear <- apply(sweep(f_gear, c(1, 3, 4), biomass, "*"),
                                c(1, 2, 3), sum)
    # yield_species_gear

    yield_species <-apply(yield_species_gear, c(1, 3), sum)

    yield_frame <- melt(yield_species)

    # leave out spin up and change units to tonnes
    # y<-yield_frame[yield_frame$time >= 1947,]

    # disregard zeroes - these were NAs only filled in to run the model

    obs<-dat$catchT
    pred<-yield_frame$value[1:18] # only selecting D.ele for now

    # sum of squared errors, could use  log-scale of predictions and data (could change this or use other error or likelihood options)

    error <- sum((log(pred[1:13]) - log(obs[1:13]))^2,na.rm=T)

    plot_dat <- data.frame(obs,pred)
    plot_dat$Year <- 2003:2020

    p <- ggplot(plot_dat, aes(x = Year)) +
        geom_line(aes(y = pred)) +
        geom_line(aes(y = obs), color = "red")
    print(p)

    p2 <- mizer::plotBiomass(sim)
    print(p2)
    return(error)
}



#' Plot ....
#'
#' @export
getBiomassFrame2 <- function (sim, species = dimnames(sim@n)$sp[!is.na(sim@params@A)], min_w = NULL,
                              start_time = as.numeric(dimnames(sim@n)[[1]][1]), end_time = as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]]),
                              ylim = c(NA, NA), total = FALSE, ...)
{
    if(is.null(min_w)) b <- getBiomass(sim, ...)
    else {
        biom_per_size <- sim@n

        if(length(min_w) == 1) # can probably condense both cases in one
        {
            # find which size class is right after user-inputed w_min
            min_w_cell <- which(as.numeric(dimnames(biom_per_size)$w) >= min_w)[1]
            b <- apply(biom_per_size[,,min_w_cell:dim(biom_per_size)[3]],c(1,2),sum)
        } else if (length(min_w) == dim(biom_per_size)[2]){
            # find which size class is right after user-inputed w_min for each species
            min_w_cell <- NULL
            for(iW in min_w) min_w_cell <- c(min_w_cell,which(as.numeric(dimnames(biom_per_size)$w) >= iW)[1])
            # remove size before w_min
            for(iSpecies in 1:dim(biom_per_size)[2]) biom_per_size[,iSpecies,1: (min_w_cell[iSpecies]-1)] <- 0

            b <- apply(biom_per_size,c(1,2),sum)
        }

    }


    if (start_time >= end_time) {
        stop("start_time must be less than end_time")
    }
    b <- b[(as.numeric(dimnames(b)[[1]]) >= start_time) & (as.numeric(dimnames(b)[[1]]) <=
                                                               end_time), , drop = FALSE]
    b_total <- rowSums(b)
    if (total) {
        b <- cbind(b, Total = b_total)
        species <- c("Total", species)
    }
    bm <- mizer::melt(b)
    min_value <- 1e-20
    bm <- bm[bm$value >= min_value & (is.na(ylim[1]) | bm$value >=
                                          ylim[1]) & (is.na(ylim[2]) | bm$value <= ylim[2]), ]
    names(bm) <- c("Year", "Species", "Biomass")
    species_levels <- c(dimnames(sim@n)$sp, "Background", "Resource",
                        "Total")
    bm$Species <- factor(bm$Species, levels = species_levels)
    bm <- bm[bm$Species %in% species, ]
    return(bm)
}


#' Get the diet composition
#'
#' The diet \eqn{D_{ij}(w, w_p)} is the prey biomass density rate for a predator of
#' species \eqn{i} and weight \eqn{w}, resolved by prey species \eqn{j} and prey
#' size \eqn{w_p}. It is calculated from the predation kernel \eqn{\phi(w, w_p)},
#' the search volume \eqn{\gamma_i(w)}, the feeding level \eqn{f_i(w)}, the
#' species interaction matrix \eqn{\theta_{ij}} and the prey abundance density
#' \eqn{N_j(w)}:
#' \deqn{
#' D_{ij}(w, w_p) = (1-f_i(w)) \gamma_i(w) \theta_{ij} N_j(w_p)
#' \phi_i(w, w_p) w_p.
#' }
#' The prey index \eqn{j} can run over all species and the resource. The returned
#' values have units of 1/year.
#'
#' The total rate \eqn{D_{ij}(w)} at which a predator of species \eqn{i}
#' and size \eqn{w} consumes biomass from prey species \eqn{j} is
#' obtained by integrating over prey sizes:
#' \deqn{
#' D_{ij}(w) = \int D_{ij}(w, w_p) dw_p.
#' }
#' This aggregated diet can also be obtained directly from the `getDiet()` function.
#'
#' @param sim An object of class \linkS4class{MizerSim}
#' @return An array (predator species x predator size x
#'    (prey species + resource) x prey size)
#'
#' @export


getDietComp<- function(sim)
{
    # initialisation
    object <- sim@params
    feedinglevel=getFeedingLevel(object)
    pred_kernel <- getPredKernel(object)
    n = sim@n[dim(sim@n)[1],,]
    n_pp = sim@n_pp[dim(sim@n_pp)[1],]
    no_sp <- dim(object@species_params)[1]
    no_w <- length(object@w)
    no_w_full <- length(object@w_full)

    diet_comp<-array(0, c(no_sp, no_w, no_sp + 1, no_w_full),
                     dimnames=list( predator=as.character(object@species_params$species), pred_size = object@w,
                                    prey = c(as.character(object@species_params$species), "background"),
                                    prey_size = object@w_full))

    # Biomass by species
    n_total_in_size_bins<- sweep(n, 2, object@dw , "*")
    b_tot <- sweep(n_total_in_size_bins, 2, object@w , "*")

    # Index of predator size classes
    idx_sp<- object@w_full %in% object@w

    #  pred_kernel * interaction matrix
    for(iW in 1:no_w){
        for(iSpecies in 1:no_sp){
            diet_comp[iSpecies,iW,1:no_sp,idx_sp]<- sweep(sweep( b_tot, c(1), object@interaction[iSpecies, 1:no_sp], "*"), c(2),
                                                          pred_kernel[iSpecies,iW,idx_sp], "*")
        }
    }
    # Search rate *  feeding level * prey biomass
    diet_comp[,,1:no_sp,]<- sweep(sweep(sweep(diet_comp[,,1:no_sp,], c(1,2), object@search_vol,"*"),
                                        c(1,2),1-feedinglevel,"*"),
                                  c(1,2),b_tot,"*")  # Prey eaten: total g prey/ year  (given predator biomass density)

    # no interaction matrix for background spectrum
    b_background <- (sweep(pred_kernel[,,], c(3), object@dw_full*object@w_full*n_pp, "*"))
    #Search rate *  feeding level * predator biomass
    b_background<- sweep(b_background, c(1,2), object@search_vol,"*") #Scale up by search volume
    b_background<- sweep(b_background, c(1,2), feedinglevel,"*") # Scale according to feeding level. Prey eaten: g prey / year / g predator
    b_background_tot<-sweep(b_background,c(1,2), b_tot, "*") # Prey eaten: total g prey/ year  (given predator biomass density)

    # Store background eaten
    diet_comp[,,no_sp+1,]<- b_background_tot

    return(diet_comp)
}


getDietMizer <-
    function (params, n = initialN(params), n_pp = initialNResource(params),
              n_other = initialNOther(params), proportion = TRUE)
    {
        params <- validParams(params)
        species <- params@species_params$species
        no_sp <- length(species)
        no_w <- length(params@w)
        no_w_full <- length(params@w_full)
        no_other <- length(params@other_encounter)
        other_names <- names(params@other_encounter)
        # assert_that(identical(dim(n), c(no_sp, no_w)), length(n_pp) ==
        #               no_w_full)
        diet <- array(0, dim = c(no_sp, no_w, no_sp + 1 + no_other),
                      dimnames = list(predator = species, w = dimnames(params@initial_n)$w,
                                      prey = c(as.character(species), "Resource", other_names)))
        idx_sp <- (no_w_full - no_w + 1):no_w_full
        if (length(params@ft_pred_kernel_e) == 1) {
            ae <- matrix(params@pred_kernel[, , idx_sp, drop = FALSE],
                         ncol = no_w) %*% t(sweep(n, 2, params@w * params@dw,
                                                  "*"))
            diet[, , 1:no_sp] <- ae
            diet[, , no_sp + 1] <- rowSums(sweep(params@pred_kernel,
                                                 3, params@dw_full * params@w_full * n_pp, "*"), dims = 2)
        }
        else {
            prey <- matrix(0, nrow = no_sp + 1, ncol = no_w_full)
            prey[1:no_sp, idx_sp] <- sweep(n, 2, params@w * params@dw, "*")
            prey[no_sp + 1, ] <- n_pp * params@w_full * params@dw_full
            ft <- array(rep(params@ft_pred_kernel_e, times = no_sp + 1) * rep(mvfft(t(prey)), each = no_sp), dim = c(no_sp, no_w_full, no_sp + 1))
            ft <- matrix(aperm(ft, c(2, 1, 3)), nrow = no_w_full)
            ae <- array(Re(mvfft(ft, inverse = TRUE)/no_w_full),
                        dim = c(no_w_full, no_sp, no_sp + 1))
            ae <- ae[idx_sp, , , drop = FALSE]
            ae <- aperm(ae, c(2, 1, 3))
            ae[ae < 1e-18] <- 0
            diet[, , 1:(no_sp + 1)] <- ae
        }
        inter <- cbind(params@interaction, params@species_params$interaction_resource)
        diet[, , 1:(no_sp + 1)] <- sweep(sweep(diet[, , 1:(no_sp +
                                                               1), drop = FALSE], c(1, 3), inter, "*"), c(1, 2), params@search_vol,
                                         "*")
        for (i in seq_along(params@other_encounter)) {
            diet[, , no_sp + 1 + i] <- do.call(params@other_encounter[[i]],
                                               list(params = params, n = n, n_pp = n_pp, n_other = n_other,
                                                    component = names(params@other_encounter)[[i]]))
        }
        f <- getFeedingLevel(params, n, n_pp)
        fish_mask <- n > 0
        diet <- sweep(diet, c(1, 2), (1 - f) * fish_mask, "*")
        if (proportion) {
            total <- rowSums(diet, dims = 2)
            diet <- sweep(diet, c(1, 2), total, "/")
            diet[is.nan(diet)] <- 0
        }
        return(diet)
    }

