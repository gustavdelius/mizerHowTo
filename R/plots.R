# scripts grouping all the plot functions

cleanTheme <- function(font_size = 9){
    theme(text = element_text(size=font_size),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "gray"),
          panel.border = element_rect(colour = "gray", fill=NA, size=.5),
          legend.position = "none",
          legend.key = element_rect(fill = "white"))
}

noXTheme <- function(){
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
}

#' Summary plot displaying 8 differents plots useful to assess the state of
#' a mizerSim object.
#'
#' Left column (top to bottom) is: size spectrum, feeding level, predation
#' and fishing mortality on linear scale, predation and fishing mortality
#' on a log scale
#' Right column (top to bottom) is: predicted catches and biomass, RDD/RDI,
#' biomass, PPMR per size of predator (work in progress)
#' @param sim An object of class MizerSim
#' @param power The abundance is plotted as the number density times the weight raised to power.
#' The default power = 1 gives the biomass density, whereas power = 2 gives the biomass density
#' with respect to logarithmic size bins.
#' @param save_it Boleean value that determines whether to save the output of the function or not.
#' Default is FALSE
#' @param name_save Character string to give a specific name if saving the plot. Default is NULL
#' @param font_size Determines the font size of the labels and ticks. Default is 9.
#'
#' @export
plotSummary <- function (sim, power = 1, save_it = FALSE, name_save = NULL, font_size = 9)
{
    # need to display the legend at the bottom and only p1 has the background so using that one
    p1 <- mizer::plotSpectra(sim, power = power, wlim = c(0.001,NA), total = TRUE)
    p1 <- p1 +
        theme(legend.key = element_rect(fill = "white")) +
        guides(color = guide_legend(nrow=1))

    mylegend<-g_legend(p1) # save the legend
    p1 <- p1 +
        cleanTheme(font_size) +
        noXTheme()

    p2 <- plotFeedingLevel(sim, include_critical = T)
    p2 <- p2 +
        cleanTheme(font_size) +
        noXTheme()

    p3 <- plotPredMort(sim)
    p3 <- p3 +
        cleanTheme(font_size) +
        noXTheme()

    p4 <- plotFMort(sim)
    p4 <- p4 +
        cleanTheme(font_size)

    # yield and ssb |

    p5 <- plotBiomassVsCatch(sim)
    p5 <- p5 +
        cleanTheme(font_size) +
        noXTheme()

    #RDI / RDD

    p6 <- plotRdiVsRdd(sim)
    p6 <- p6 +
        cleanTheme(font_size)


    p7 <- plotBiomass(sim)
    p7 <- p7 +
        cleanTheme(font_size)

    # predator / prey mass comparison

    p8 <- plotPPMR(sim)
    p8 <- p8 +
        cleanTheme(font_size)


    plots_arranged <- plot_grid(p1,p2,p3,p4,p5, p6, p7,p8, byrow = F,
                                rel_widths = c(3,3), nrow = 4,
                                align = "v")

    p10 <- plot_grid(plots_arranged, mylegend,
                     rel_heights = c(10,1),
                     ncol = 1)

    if(save_it & !is.null(name_save)) ggsave(p10, filename = paste(name_save,".png",sep=""), units = "cm", width = 21, height = 29)
    else if (save_it & is.null(name_save)) ggsave(p10, filename = "tempSummary.png", units = "cm", width = 21, height = 29)

    return(p10)
}


#' @title Plot showing emergent PPMR
#'
#' @description Yet another plot
#'
#' @param sim An object of class MizerSim
#' @param return_data A boolean value that determines whether the formatted data
#'   used for the plot is returned instead of the plot itself. Default value is
#'   FALSE
#'
#' @export

plotPPMR <- function(sim, return_data = FALSE){

    diet_dat <- getDietComp(sim)
    SpIdx <- sim@params@species_params$species
    tempSimDf <- NULL

    for(iSpecies in SpIdx) # for each species
    {
        diet_dat_sp <- diet_dat[iSpecies,,,]
        diet_dat_sp<- apply(diet_dat_sp,c(1,3),sum) # sum prey identity, keep size class
        speciesPPMR <- NULL
        size_name_vec <- NULL
        size_preferred <- NULL
        for(iSize in dimnames(diet_dat_sp)$pred_size) # for each size class need PPMR value
        {
            if(sum(diet_dat_sp[iSize,])) # if there is at least one diet data
            {
                size_name_vec <- c(size_name_vec,iSize)
                sizeDat <- diet_dat_sp[iSize,] # select the size
                densityDat <- sizeDat / as.numeric(as.character(names(sizeDat)))# adjust biomass > density

                # calculating realised PPMR
                PreferredSizeClass <- which(densityDat == max(densityDat)) # which size class is most feed upon
                sizePPMR <- as.numeric(iSize)/as.numeric(names(PreferredSizeClass)) # calculate PPMR
                speciesPPMR <- c(speciesPPMR,sizePPMR)

                # what's the favorite mass? (taking the name of the size class)
                size_preferred <- c(size_preferred,as.numeric(names(PreferredSizeClass)))

                # what's the mean mass? converting from biomass in bin to mass I guess
                # how to calculate the average from a set of discrete values? I would need to duplicate the discrete classes by the biomass number (or individual whatever)
                # and then calculate the mean from that, is it legit?
                # for now, simple soluttion, mean mass is most eaten mass (assuming normal distribution)
                # temp <- sizeDat / sim@params@w / sim@params@dw
                # f1n <- MASS::fitdistr(sizeDat,"normal")
                # mean(temp)
                # c <- hist(sizeDat)
                # size_mean <- c(size_mean,mean(sizeDat[sizeDat != 0]))

            }
        }
        tempSpeciesDf <- data.frame("species" = rep(iSpecies,length(speciesPPMR)), "w" = as.numeric(size_name_vec), "rPPMR" = speciesPPMR, "prey_mass" = size_preferred)
        tempSimDf <- rbind(tempSimDf,tempSpeciesDf) # create a df of species
    }

    # plottin the data
    plot_dat <- tempSimDf

    legend_col <- intersect(sim@params@species_params$species,
                            names(sim@params@linecolour))

    p <- ggplot(plot_dat) +
        geom_line(aes(x = w, y = prey_mass, color = species)) +
        scale_x_continuous(name = "Predator mass (g)", trans = "log10") +
        scale_y_continuous(name = "Normalised PPMR", trans = "log10") +
        scale_colour_manual(values = sim@params@linecolour[legend_col])

    if(return_data) return(plot_dat) else return(p)
}




#' @title Plot displaying biomass versus catch per species
#'
#' @description Yet another plot
#'
#' @inheritParams plotPPMR
#'
#' @export

plotBiomassVsCatch <- function(sim, return_data = FALSE){

    bm <- getBiomassFrame2(sim)
    plot_dat <- filter(bm, Year == max(unique(bm$Year)))
    yieldDat <- getYield(sim)
    plot_dat$yield <- yieldDat[dim(yieldDat)[1],]
    plot_dat$Year <- NULL
    colnames(plot_dat) <- c("Species", "average biomass", "average catch")
    plot_dat$Species <- factor(as.character(plot_dat$Species),
                               levels = c(as.character(sim@params@species_params$species)))
    plot_dat <- melt(plot_dat,"Species")
    plot_dat$w_inf <- rep(sim@params@species_params$w_inf,2)

    legend_col <- intersect(sim@params@species_params$species,
                            names(sim@params@linecolour))

    p <- ggplot(plot_dat)+
        geom_point(aes(x = w_inf, y = value, color = Species, shape = variable), size = 6, alpha = .8) +
        ggrepel::geom_text_repel(data = filter(plot_dat,variable == "average biomass"),
                                 aes(x = w_inf, y = value, label = Species), hjust = 0, nudge_x = 0.05)+
        geom_line(aes(x = w_inf, y = value, color = Species)) +
        scale_y_continuous(name = "Catch and Biomass", limits = c(NA,NA), trans = "log10") +
        scale_x_continuous(name = "Asymptotic size (g)", trans = "log10") +
        scale_colour_manual(values = sim@params@linecolour[legend_col]) +
        scale_shape_manual(name = "Data", values = c(16,17))


    if(return_data) return(plot_dat) else return(p)

}

#' @title Plot density independent versus density dependent reproduction rate
#'
#' @description yet another plot
#'
#' @param object An object of class \linkS4class{MizerSim} or
#'   \linkS4class{MizerParams}.
#' @param return_data A boolean value that determines whether the formatted data
#'   used for the plot is returned instead of the plot itself. Default value is
#'   FALSE
#'
#' @export

plotRdiVsRdd <- function(object, return_data = FALSE){

    if (is(object, "MizerSim")) {
        params <- object@params
        params <- setInitialValues(params, object)
    } else if (is(object, "MizerParams")) {
        params <- validParams(object)
    }

    plot_dat <- as.data.frame(getRDI(params)/getRDD(params))
    plot_dat$species <- factor(rownames(plot_dat),params@species_params$species)
    colnames(plot_dat)[1] <- "ratio"
    plot_dat$w_inf <- params@species_params$w_inf

    legend_col <- intersect(params@species_params$species,
                            names(params@linecolour))

    p <- ggplot(plot_dat)+
        geom_point(aes(x = w_inf, y = ratio, color = species), size = 6, alpha = .8) +
        ggrepel::geom_text_repel(aes(x = w_inf, y = ratio, label = species), hjust = 0, nudge_x = 0.05)+
        # geom_line(aes(x = w_inf, y = value, color = Species)) +
        scale_y_continuous(name = "Density-independent / density-dependent reproduction rate", trans = "log10") +
        scale_x_continuous(name = "Asymptotic size (g)", trans = "log10") +
        scale_colour_manual(values = params@linecolour[legend_col])

    if(return_data) return(plot_dat) else return(p)

}

#' A function that plots the predicted yield versus observed yield.
#'
#' @inheritParams plotPPMR
#' @param dat A dataframe containing the observed yield values
#'
#' @export
plotPredObsYield <-function(sim, dat, returnData = FALSE){
    ## check obs vs. predicted yield
    plot_dat <-melt(getYield(sim)[100,]/1e6)
    plot_dat$obs <- log10(dat)
    plot_dat$value <- log10(plot_dat$value)
    plot_dat$Species <-row.names(plot_dat)

    w_inf <- log10(sim@params@species_params$w_inf)
    names(w_inf) <- sim@params@species_params$species

    # window size
    winLim <- c(min(plot_dat$obs,plot_dat$value), max(plot_dat$obs,plot_dat$value))
    winLim <- c(0,max(plot_dat$obs,plot_dat$value)) # abline doesn't show anymore 18/06/2021
    p <- ggplot(plot_dat) + # plot predicted and observed yields
        geom_point(aes(x = value, y = obs, color = Species, size = Species)) +
        ggrepel::geom_text_repel(aes(x = value, y = obs, label = Species), hjust = 0, nudge_x = 0.05)+
        scale_size_manual(values = w_inf) +
        scale_color_manual(values = sim@params@linecolour) +
        geom_abline(color = "black", slope = 1, intercept = 0, linetype = "dashed", alpha = .5) +
        scale_x_continuous(name = "log10 Predicted Yield in t/year", limits = winLim) +
        scale_y_continuous(name = "log10 Observed Yield in t/year", limits = winLim) +
        theme(legend.position = "none", legend.key = element_rect(fill = "white"),
              panel.background = element_blank(), panel.grid.minor = element_line(color = "gray"))

    if(returnData) return(plot_dat) else return(p)
}

#' Function showing the diet proportion of predators
#'
#' @param sim An object if class MizerSim
#' @param species A character string of the species name. If NULL, all species
#' are shwon as facets. Default is NULL
#' @param xlim A numeric vector of length two providing lower and upper limits for the x axis.
#' Use NA to refer to the existing minimum or maximum. Default is c(1,NA).
#' @param returnData A boolean value that determines whether to return the plot
#' or the data itself. Default is FALSE
#'
#' @export
plotDiet2 <- function (sim, species = NULL, xlim = c(1,NA), returnData = F)
{
    params <- sim@params
    # if (is.integer(species)) {
    #     species <- params@species_params$species[species]
    # }


    # diet <- getDiet(params)[params@species_params$species ==
    #     species, , ]
    # prey <- dimnames(diet)$prey
    # prey <- factor(prey, levels = rev(prey))
    # plot_dat <- data.frame(Proportion = c(diet), w = params@w,
    #     Prey = rep(prey, each = length(params@w)))
    # plot_dat <- plot_dat[plot_dat$Proportion > 0, ]
    #
    # ggplot(plot_dat) + geom_area(aes(x = w, y = Proportion, fill = Prey)) +
    #     scale_x_log10(limits = xlim) + labs(x = "Size [g]") +
    #   scale_fill_manual(values = sim@params@linecolour) +
    #   ggtitle(species)


    diet <- getDiet(params)
    plot_dat <- melt(diet)
    plot_dat <- plot_dat[plot_dat$value > 0, ]
    colnames(plot_dat) <- c("Predator", "size", "Prey", "Proportion")

    if(is.null(species)) p <- ggplot(plot_dat) + facet_wrap(.~Predator, scales = "free") else p <- ggplot(filter(plot_dat, Predator == species))

    p <- p +
        geom_area(aes(x = size, y = Proportion, fill = Prey))+
        scale_x_continuous(limits = c(1,NA), name = "Size [g]", trans = "log10") +
        scale_fill_manual(values = sim@params@linecolour)+
        theme(legend.position = "right", legend.key = element_rect(fill = "white"),
              panel.background = element_blank(), panel.grid.minor = element_line(color = "gray"),
              strip.background = element_blank())

    if(returnData) return(plot_dat) else return(p)

}

#' A function plotting the fisheries effort versus yield of a species.
#'
#' @param params An object of class MizerParams
#' @param effortRes A numeric value determinin the number of simulation
#' to be run per species. A high number offers a better resolution of
#' the plot albeit for a longer computing time. Default is 20
#' @param returnData A boolean value that determines whether to return the plot
#' or the data itself. Default is FALSE
#' @param speciesData A list of size 2. If provided, the function uses the
#' first slot (character string of a species name) to compute only that particular
#' species and update it in the second slot (a plotFmsy(returnData = TRUE) output).
#' Used in shiny apps to decrease the run time. Default is NULL
#'
#' @export
plotFmsy <- function(params, effortRes = 20, returnData = F, speciesData = NULL)
{
    # make one gear per species so we can vary the effort per species
    gear <- gear_params(params)
    gear$gear <- params@species_params$species
    gear_params(params) <- gear

    catchability <- params@species_params$catchability

    xlim <- 1.5 # maximum effort* catchability / xaxis limit

    # we want to vary effort value so we get a scale from 0 to 1 of effort * catchability per species

    # the "species" arg allows to run the function for only one species, which should be faster but it means "species" must also contain the result of every other species (so it's a two object list)

    if(!is.null(speciesData))
    {
        speciesName <- speciesData[[1]] # which species are we changing?
        plot_dat <- speciesData[[2]] # plot_dat of all species
        plot_dat <- filter(plot_dat, species!= speciesName) # remove previous result of the concerned species
        iSpecies <- which(params@species_params$species == speciesName)
        counter = 0 # sim counter
        # determine effort range
        effortMax <- round(xlim/catchability[iSpecies],1)+.1
        SpDat <- NULL

        effortSeq <- exp(seq(0,log(effortMax+1), length.out =  effortRes)) -1
        effortSeq <- effortSeq[effortSeq<effortMax] # creating an exponentially increasing effort sequence
        for(iEffort in effortSeq)
        {
            effort_vec <- rep(1,dim(params@species_params)[1]) # all effort set to one
            effort_vec[iSpecies] <- iEffort # except that one which varies

            if(!counter )
            {
                tempSim <- project(params, effort = effort_vec, t_max = 20)
                counter <- 1
            } else  tempSim <- project(params, effort = effort_vec, t_max = 10, initial_n = tempSim@n[dim(tempSim@n)[1],,],
                                       initial_npp = tempSim@n_pp[dim(tempSim@n_pp)[1],])
            #catch
            yieldDat <- getYield(tempSim)
            SpDat <- rbind(SpDat,c(yieldDat[dim(yieldDat)[1],iSpecies],iEffort))
        }
        SpDat <- as.data.frame(SpDat)
        SpDat$species <- params@species_params$species[iSpecies]
        SpDat$V2 <- SpDat$V2*catchability[iSpecies] # so V2 is effort * catchability
        colnames(SpDat) <- c("yield","effort","species")
        plot_dat <- rbind(plot_dat,SpDat)



    } else {

        plot_dat <- NULL
        for(iSpecies in 1:dim(params@species_params)[1])
        {
            counter = 0 # sim counter
            # determine effort range
            effortMax <- round(xlim/catchability[iSpecies],1)+.1
            SpDat <- NULL

            effortSeq <- exp(seq(0,log(effortMax+1), length.out =  effortRes)) -1 # every .1 takes 2 min to run, evry .2 takes 1 min but lesser resolution
            effortSeq <- effortSeq[effortSeq<effortMax] # creating an exponentially increasing effort sequence
            for(iEffort in effortSeq)
            {
                effort_vec <- rep(1,dim(params@species_params)[1]) # all effort set to one
                effort_vec[iSpecies] <- iEffort # except that one which varies

                if(!counter )
                {
                    tempSim <- project(params, effort = effort_vec, t_max = 20)
                    counter <- 1
                } else  tempSim <- project(params, effort = effort_vec, t_max = 10, initial_n = tempSim@n[dim(tempSim@n)[1],,],
                                           initial_npp = tempSim@n_pp[dim(tempSim@n_pp)[1],])
                #catch
                yieldDat <- getYield(tempSim)
                SpDat <- rbind(SpDat,c(yieldDat[dim(yieldDat)[1],iSpecies],iEffort))
            }
            SpDat <- as.data.frame(SpDat)
            SpDat$species <- params@species_params$species[iSpecies]
            SpDat$V2 <- SpDat$V2*catchability[iSpecies] # so V2 is effort * catchability
            colnames(SpDat) <- c("yield","effort","species")
            plot_dat <- rbind(plot_dat,SpDat)

        }
    }

    plot_dat$species <- factor(plot_dat$species, levels = params@species_params$species)
    # colnames(plot_dat) <- c("yield","effort","species")
    if(!is.null(speciesData)) p <- ggplot(filter(plot_dat, species == speciesName)) else p <- ggplot(plot_dat)

    p <- p + geom_line(aes(x = effort , y = yield, color = species))+
        facet_wrap(species~., scales = "free") +
        scale_x_continuous(limits= c(0,xlim),name = "fishing mortality rate")+#, limits = c(1e10,NA))+
        scale_y_continuous(trans = "log10") +
        scale_color_manual(name = "Species", values = params@linecolour) +
        theme(legend.position = "none", legend.key = element_rect(fill = "white"),
              panel.background = element_blank(), panel.grid.minor = element_line(color = "gray"),
              strip.background = element_blank())

    if(returnData) return(plot_dat) else return(p)

}


#' Plot growth curves giving weight as a function of age
#'
#' When the growth curve for only a single species is plotted, horizontal lines are included
#' that indicate the maturity size and the maximum size for that species. If furthermore the
#' species parameters contain the variables a and b for length to weight conversion and the
#' von Bertalanffy parameter k_vb (and optionally t0), then the von Bertalanffy growth curve
#' is superimposed in black.
#' @param object An object of class \linkS4class{MizerSim} or
#'   \linkS4class{MizerParams}.
#' @param return_data A boolean value that determines whether the formated data
#' used for the plot is returned instead of the plot itself. Default value is FALSE
#' @param species The species to be selected. Optional. By default all target species
#' are selected. A vector of species names, or a numeric vector with the species indices,
#' or a logical vector indicating for each species whether it is to be selected (TRUE) or not.
#' @param max_age The age up to which to run the growth curve. Default is 20.
#' @param percentage Boolean value. If TRUE, the size is given as a percentage of the maximal size.
#' @param species_panel  If TRUE, display all species with their Von Bertalanffy curves as facets
#' (need species and percentage to be set to default). Default FALSE.
#' @param highlight Name or vector of names of the species to be highlighted.
#' @param A boolean value that determines whether the formated data used for the plot is returned
#' instead of the plot itself. Default value is FALSE
#' @export
plotGrowthCurves2 <- function (object,
                               species = NULL,
                               max_age = 20,
                               percentage = FALSE,
                               species_panel = FALSE,
                               highlight = NULL,
                               returnData = F)
{
    if (is(object, "MizerSim")) {
        params <- object@params
        t <- dim(object@n)[1]
        params@initial_n[] <- object@n[t, , ]
        params@initial_n_pp <- object@n_pp[t, ]
    }
    else if (is(object, "MizerParams")) {
        params <- validParams(object)
    }
    species <- valid_species_arg(params, species)
    ws <- getGrowthCurves(params, species, max_age, percentage)
    plot_dat <- melt(ws)
    plot_dat$Species <- factor(plot_dat$Species, params@species_params$species)
    plot_dat$legend <- "model"
    if (all(c("a", "b", "k_vb") %in% names(params@species_params))) {
        if ("t0" %in% names(params@species_params)) {
            t0 <- params@species_params$t0
        }
        else {
            t0 <- 0
        }
        VBdf <- data.frame(species = params@species_params$species,
                           w_inf = params@species_params$w_inf, a = params@species_params$a,
                           b = params@species_params$b, k_vb = params@species_params$k_vb,
                           t0 = t0)
        VBdf$L_inf <- (VBdf$w_inf/VBdf$a)^(1/VBdf$b)
        plot_dat2 <- plot_dat
        plot_dat2$value <- apply(plot_dat, 1, function(x) {
            sel <- VBdf$species == x[1]
            length <- VBdf$L_inf[sel] * (1 - exp(-VBdf$k_vb[sel] *
                                                     (as.numeric(x[2]) - VBdf$t0[sel])))
            VBdf$a[sel] * length^VBdf$b[sel]
        })
        plot_dat2$legend <- "von Bertalanffy"
        plot_dat <- rbind(plot_dat, plot_dat2)
    }
    p <- ggplot(filter(plot_dat, legend == "model")) + geom_line(aes(x = Age,
                                                                     y = value, colour = Species, linetype = Species, size = Species))
    y_label <- if (percentage) "Percent of maximum size" else "Size [g]"
    linesize <- rep(0.8, length(params@linetype))
    names(linesize) <- names(params@linetype)
    linesize[highlight] <- 1.6
    p <- p + scale_x_continuous(name = "Age [Years]") + scale_y_continuous(name = y_label) +
        scale_colour_manual(values = params@linecolour) + scale_linetype_manual(values = params@linetype) +
        scale_size_manual(values = linesize)
    if (!percentage) {
        if (length(species) == 1) {
            idx <- which(params@species_params$species == species)
            w_inf <- params@species_params$w_inf[idx]
            p <- p + geom_hline(yintercept = w_inf, colour = "grey") +
                annotate("text", 0, w_inf, vjust = -1, label = "Maximum")
            w_mat <- params@species_params$w_mat[idx]
            p <- p + geom_hline(yintercept = w_mat, linetype = "dashed",
                                colour = "grey") + annotate("text", 0, w_mat,
                                                            vjust = -1, label = "Maturity")
            if ("von Bertalanffy" %in% plot_dat$legend)
                p <- p + geom_line(data = filter(plot_dat, legend ==
                                                     "von Bertalanffy"), aes(x = Age, y = value))
        }
        else if (species_panel) {
            p <- ggplot(plot_dat) +
                geom_line(aes(x = Age, y = value, colour = legend)) +
                scale_x_continuous(name = "Age [years]") +
                scale_y_continuous(name = "Size [g]") +
                facet_wrap(.~Species, scales = "free") +
                geom_hline(aes(yintercept = w_mat),
                           data = tibble(Species = as.factor(object@params@species_params$species[]),
                                         w_mat = object@params@species_params$w_mat[]),
                           linetype = "dashed", colour = "grey") +
                geom_hline(aes(yintercept = w_inf),
                           data = tibble(Species = as.factor(object@params@species_params$species[]),
                                         w_inf = object@params@species_params$w_inf[]),
                           linetype = "solid", colour = "grey") +
                theme(panel.background = element_blank(), panel.grid.minor = element_line(color = "gray"),
                      strip.background = element_blank(), legend.key = element_blank())+
                scale_color_discrete(name = "Growth", labels = c("Modelled","von Bertalanffy"))
        }
    }
    if(returnData) return(plot_dat) else return(p)
}


#' @title Summary function bundling helpful plots during calibration.
#'
#' @description  Stage 1 shows size spectra, biomass and predicted versus observed yield
#' Stage 2 shows the pannel of fisheries effort versus yield
#' Stage 3 shows the pannel of growth curves
#' @inheritParams plotSummary
#' @param stage Numeric values which determines the function's output. Range from 1 to 3.
#' Default is 1.
#' @param effortRes A numeric value determinin the number of simulation
#' to be run per species. A high number offers a better resolution of
#' the plot albeit for a longer computing time. Default is 10.
#' @param catch_dat A dataframe containing the observed yield values. If not
#' included, plotPredObsYield will not be displayed. Default is NULL.
#' @export
plotCalibration <- function(sim, catch_dat = NULL, stage = 1, wlim = c(.1,NA), power = 1, effortRes = 10, font_size = 9)
{
    xlim = c(NA,10^log10(max(sim@params@species_params$w_inf)))

    switch (stage,
            "1" = {
                p1 <- plotSpectra(sim, power = power, wlim = wlim, total = TRUE)

                p1 <- p1 +
                    theme(legend.position = "right",
                          legend.key = element_rect(fill = "white"))

                mylegend<-g_legend(p1) # save the legend

                p1 <- p1 +
                    cleanTheme(font_size)

                p2 <- plotBiomass(sim)
                p2 <- p2 +
                    cleanTheme(font_size)

                if(is.null(catch_dat))
                {
                    leftCol <- plot_grid(p1,p2,
                                         ncol = 1, align = "v")
                    p <- plot_grid(leftCol, mylegend,
                                   rel_widths = c(6,1),
                                   ncol = 2)
                } else {

                    p3 <- plotPredObsYield(sim,catch_dat)
                    p3 <- p3 +
                        cleanTheme(font_size)


                    leftCol <- plot_grid(p1,p2,p3,
                                         ncol = 1, align = "v")
                    p <- plot_grid(leftCol, mylegend,
                                   rel_widths = c(6,1),
                                   ncol = 2)
                }
            },
            "3" = {
                p <- plotFmsy(sim@params,effortRes = effortRes)
            },
            "2" = {
                p <- plotGrowthCurves2(sim, species_panel = T)

            },
            {print("Unknow stage selected.")
                p <- NULL}
    )
    return(p)
}





#' Plot ....
#'
#' @export

plot_relative_biomass = function(sim0, sim1, ratio = FALSE) {

    # assume sim0 is steady state, sim1 is some kind of variation such as fishing
    # mean of last five years
    fish_sim <- apply(N(sim1)[(dim(N(sim1))[1]-5):(dim(N(sim1))[1]),,],2:3,mean)
    unfish_sim <- apply(N(sim0)[(dim(N(sim0))[1]-5):(dim(N(sim0))[1]),,],2:3,mean)

    if (ratio) {
        relative_n <- melt(fish_sim/unfish_sim) # Julia's original calculation
    } else {
        relative_n <- melt((fish_sim - unfish_sim) / (fish_sim + unfish_sim)) # Gustav's suggested calculation
    }
    colnames(relative_n)[1] <- "Species"
    legend_levels <- intersect(names(sim0@params@linecolour), relative_n$Species)
    p <- ggplot(relative_n) +
        geom_line(aes(x = w, y = value, colour = Species), size = 1) +
        scale_x_continuous(trans = "log10", name = "Weight [g]") +
        scale_color_manual(values = sim0@params@linecolour[legend_levels]) +
        theme(legend.key = element_rect(fill = "white"))

    if (ratio == T) {
        p = p + scale_y_continuous(trans = "log10") +
            geom_hline(yintercept = 1, linetype = 1, colour="dark grey", size=0.75) +
            labs(y="Relative abundance")
    } else {
        p = p + geom_hline(yintercept = 0, linetype = 1, colour="dark grey", size=0.75) +
            labs(y="Relative difference")
    }

    print(p)

}

