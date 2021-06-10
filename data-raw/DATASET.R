# # # # Get most recent annual ICES information (assumes have already downloaded and converted and saved file from http://www.ices.dk/datacentre/StdGraphDB.asp, called "Fishdata.txt")
# # # # 1) Computes time-averaged values (here, over 1985-1995) for SSB, Landings and Fishing mortality rates for calibrating the size-based model to a reference period
# # # # 2) Creates species-specific and time-varying matrices for: SSB, Landings and Fishing mortality rates for all years available, for simulating dynamical changes through time
# # # # 3) Also calculates species-specific CVs for recruitment, for use in stochastic simulations.



##### JLB UPDATED 3/9/2020

# get most recent ICES stock assessment summary outputs for North Sea model species

# install.packages("icesSAG") # archived
# library(devtools)
# install_version("icesVocab","1.1.8")
# install_version("icesSAG","1.3-6")

library(icesSAG)
library(tidyverse)
# ?icesSAG


### get summary data for each species by assessmentkey:

#Sprat: https://standardgraphs.ices.dk/ViewSourceData.aspx?key=13322

Sprat_sumtab <- getSummaryTable(13322)

# Sandeel : https://standardgraphs.ices.dk/ViewCharts.aspx?key=13303
# https://standardgraphs.ices.dk/ViewCharts.aspx?key=13301
# https://standardgraphs.ices.dk/ViewCharts.aspx?key=13298
# https://standardgraphs.ices.dk/ViewCharts.aspx?key=13304

Sandeel_sumtab <- getSummaryTable(c(13303,13301,13298,13304))

# N.pout - https://standardgraphs.ices.dk/ViewCharts.aspx?key=13166

N.pout_sumtab <- getSummaryTable(c(13166))


# Herring - https://standardgraphs.ices.dk/ViewCharts.aspx?key=13422

Herring_sumtab <- getSummaryTable(c(13422))

# Dab - https://standardgraphs.ices.dk/ViewCharts.aspx?key=13184

Dab_sumtab <- getSummaryTable(c(13184))

# Whiting - https://standardgraphs.ices.dk/ViewCharts.aspx?key=13525
# https://standardgraphs.ices.dk/ViewCharts.aspx?key=13506

Whiting_sumtab <- getSummaryTable(c(13525))

# Sole - https://standardgraphs.ices.dk/ViewCharts.aspx?key=13743
# https://standardgraphs.ices.dk/ViewCharts.aspx?key=13495
# https://standardgraphs.ices.dk/ViewCharts.aspx?key=13828

Sole_sumtab <- getSummaryTable(c(13743))

# Gurnard - https://standardgraphs.ices.dk/ViewCharts.aspx?key=13493

Gurnard_sumtab <- getSummaryTable(c(13493))

# Plaice - https://standardgraphs.ices.dk/ViewCharts.aspx?key=13484
# https://standardgraphs.ices.dk/ViewCharts.aspx?key=13744

Plaice_sumtab <- getSummaryTable(c(13484))

# Haddock - https://standardgraphs.ices.dk/ViewCharts.aspx?key=13204

Haddock_sumtab <- getSummaryTable(c(13204))

# Cod - https://standardgraphs.ices.dk/ViewCharts.aspx?key=13740
# https://standardgraphs.ices.dk/ViewCharts.aspx?key=13838

Cod_sumtab <- getSummaryTable(c(13740,13838))

# Saithe - https://standardgraphs.ices.dk/ViewCharts.aspx?key=13511

Saithe_sumtab <- getSummaryTable(c(13511))



############ Need to use the data above and create new versions of the following matrices (rows: year, cols: species):

SpIdx <- c("Sprat","Sandeel","N.pout","Herring","Dab","Whiting","Sole","Gurnard","Plaice","Haddock","Cod","Saithe")
# dataIdx <- c("Sprat","Sandeel","N.pout","Herring","Dab","Whiting","Sole","Gurnard","Plaice","Haddock","Cod","Saithe")

# Which data set is the longest?
dim_sum <- NULL
for(iSpecies in SpIdx)
{
  switch (iSpecies,
          Sprat = {myVar = Sprat_sumtab[[1]]},
          Sandeel = {myVar = Sandeel_sumtab[[1]]},
          N.pout = {myVar = N.pout_sumtab[[1]]},
          Herring = {myVar = Herring_sumtab[[1]]},
          Dab = {myVar = Dab_sumtab[[1]]},
          Whiting = {myVar = Whiting_sumtab[[1]]},
          Sole = {myVar = Sole_sumtab[[1]]},
          Gurnard = {myVar = Gurnard_sumtab[[1]]},
          Plaice = {myVar = Plaice_sumtab[[1]]},
          Haddock = {myVar = Haddock_sumtab[[1]]},
          Cod = {myVar = Cod_sumtab[[1]]},
          Saithe = {myVar = Saithe_sumtab[[1]]},
          {}
  )
  dim_sum <- c(dim_sum,dim(myVar)[1])
  print(max(dim_sum))
}
# Herring it is
timePeriod <- Herring_sumtab[[1]]$Year

# Extract "F" column and create time series of fishing mortalty rates: fmat.csv

fmat <- matrix(NA,nrow = length(timePeriod), ncol = length(SpIdx), dimnames = list("year" = timePeriod, "species" = SpIdx))

for(iSpecies in SpIdx)
{
  switch (iSpecies,
          Sprat = {myVar = Sprat_sumtab[[1]]$F},
          Sandeel = {
            myMat <- matrix(NA, nrow = dim(Sandeel_sumtab[[1]])[1], ncol = 4)
            mySSB <- matrix(NA, nrow = dim(Sandeel_sumtab[[1]])[1], ncol = 4)

            for(icol in 1:4)
            {
              tempF <- Sandeel_sumtab[[icol]]$F*Sandeel_sumtab[[icol]]$SSB
              while(length(tempF)<dim(Sandeel_sumtab[[1]])[1]) tempF <- c(NA,tempF)
              myMat[,icol] <- tempF

              tempSSB <- Sandeel_sumtab[[icol]]$SSB
              while(length(tempSSB)<dim(Sandeel_sumtab[[1]])[1]) tempSSB <- c(NA,tempSSB)
              mySSB[,icol] <- tempSSB

            }
            ssbSum <- apply(mySSB,1,sum,na.rm=T)
            myVar <- apply(myMat,1,mean,na.rm=T)/ssbSum
            myVar[is.nan(myVar)] <- NA # in case there was a mean of NAs only, keeping the NA in place
          },
          N.pout = {myVar = N.pout_sumtab[[1]]$F},
          Herring = {myVar = Herring_sumtab[[1]]$F},
          Dab = {myVar = Dab_sumtab[[1]]$F},
          Whiting = {myVar = Whiting_sumtab[[1]]$F},
          Sole = {myVar = Sole_sumtab[[1]]$F},
          Gurnard = {myVar = Gurnard_sumtab[[1]]$F},
          Plaice = {myVar = Plaice_sumtab[[1]]$F},
          Haddock = {myVar = Haddock_sumtab[[1]]$F},
          Cod = {
            myMat <- matrix(NA, nrow = dim(Cod_sumtab[[1]])[1], ncol = 2)
            mySSB <- matrix(NA, nrow = dim(Cod_sumtab[[1]])[1], ncol = 2)

            for(icol in 1:2)
            {
              tempF <- Cod_sumtab[[icol]]$F*Cod_sumtab[[icol]]$SSB
              while(length(tempF)<dim(Cod_sumtab[[1]])[1]) tempF <- c(NA,tempF)
              myMat[,icol] <- tempF

              tempSSB <- Cod_sumtab[[icol]]$SSB
              while(length(tempSSB)<dim(Cod_sumtab[[1]])[1]) tempSSB <- c(NA,tempSSB)
              mySSB[,icol] <- tempSSB

            }
            ssbSum <- apply(mySSB,1,sum,na.rm=T)
            myVar <- apply(myMat,1,mean,na.rm=T)/ssbSum
            myVar[is.nan(myVar)] <- NA # in case there was a mean of NAs only, keeping the NA in place
          },
          Saithe = {myVar = Saithe_sumtab[[1]]$F},
          {}
  )
  while (length(myVar)<length(timePeriod)) myVar <- c(NA,myVar)
  fmat[,iSpecies] <- myVar
}
write.csv(fmat, file = "data/fmatWeighted.csv")

#nseaparams.csv , "catachabliity column" - here puttinhg in the time-averaged F as a baseline reference value

fAvg <- fmat[which(rownames(fmat) == "2014"):which(rownames(fmat) == "2019"),]
fAvg <- apply(fAvg,2,mean,na.rm=T)
fAvg[is.nan(fAvg)] <- NA

nsparams <- readr::read_csv("data/old/nsparams.csv")
GurnardCatch <- nsparams$catchability[8]
nsparams$catchability <- fAvg
nsparams$catchability[8] <- GurnardCatch # there is no data so using catchability from old df to not have NAs
write.csv(nsparams, file = "data/nsparams.csv")
# Extract "catches" column and create time series of total catches (including discards) | RF catches column is already total catches

catchesMat <- matrix(NA,nrow = length(timePeriod), ncol = length(SpIdx), dimnames = list("year" = timePeriod, "species" = SpIdx))

for(iSpecies in SpIdx)
{
  switch (iSpecies,
          Sprat = {myVar = Sprat_sumtab[[1]]$catches},
          Sandeel = {
            myMat <- matrix(NA, nrow = dim(Sandeel_sumtab[[1]])[1], ncol = 4)
            for(icol in 1:4)
            {
              tempF <- Sandeel_sumtab[[icol]]$catches
              while(length(tempF)<dim(Sandeel_sumtab[[1]])[1]) tempF <- c(NA,tempF)
              myMat[,icol] <- tempF
            }
            myVar <- apply(myMat,1,sum,na.rm=T)
            myVar[which(myVar == 0)] <- NA # summing NA will do this, keeping the NA in place
          },
          N.pout = {myVar = N.pout_sumtab[[1]]$landings}, # N.pout does not have catches data, just landings (and no discards)
          Herring = {myVar = Herring_sumtab[[1]]$catches},
          Dab = {myVar = Dab_sumtab[[1]]$catches},
          Whiting = {myVar = Whiting_sumtab[[1]]$catches},
          Sole = {myVar = Sole_sumtab[[1]]$catches},
          Gurnard = {myVar = Gurnard_sumtab[[1]]$catches},
          Plaice = {myVar = Plaice_sumtab[[1]]$catches},
          Haddock = {myVar = Haddock_sumtab[[1]]$catches},
          Cod = {
            myMat <- matrix(NA, nrow = dim(Cod_sumtab[[1]])[1], ncol = 2)
            for(icol in 1:2)
            {
              tempF <- Cod_sumtab[[icol]]$catches
              while(length(tempF)<dim(Cod_sumtab[[1]])[1]) tempF <- c(NA,tempF)
              myMat[,icol] <- tempF
            }
            myVar <- apply(myMat,1,sum,na.rm=T)
            myVar[which(myVar == 0)] <- NA # summing NA will do this, keeping the NA in place
          },
          Saithe = {myVar = Saithe_sumtab[[1]]$catches},
          {}
  )
  while (length(myVar)<length(timePeriod)) myVar <- c(NA,myVar)
  catchesMat[,iSpecies] <- myVar
}
write.csv(catchesMat, file = "data/catchesMat.csv") # loading from csv won't work with below code due to losing "years" as rownames
#time-averaged-catches.csv
# averaged subset
catchAvg <- catchesMat[which(rownames(catchesMat) == "2014"):which(rownames(catchesMat) == "2019"),]
catchAvg <- apply(catchAvg,2,mean,na.rm=T)
catchAvg[is.nan(catchAvg)] <- NA
catchAvg <- data.frame("species" = SpIdx, "Catch_1419_tonnes" = catchAvg,row.names = NULL)

write.csv(catchAvg, file = "data/time-averaged-catches.csv",row.names = F)

# Extract "SSB" column and create time series of SSB

SSBmat <- matrix(NA,nrow = length(timePeriod), ncol = length(SpIdx), dimnames = list("year" = timePeriod, "species" = SpIdx))

for(iSpecies in SpIdx)
{
  switch (iSpecies,
          Sprat = {myVar = Sprat_sumtab[[1]]$SSB},
          Sandeel = {
            myMat <- matrix(NA, nrow = dim(Sandeel_sumtab[[1]])[1], ncol = 4)
            for(icol in 1:4)
            {
              tempF <- Sandeel_sumtab[[icol]]$SSB
              while(length(tempF)<dim(Sandeel_sumtab[[1]])[1]) tempF <- c(NA,tempF)
              myMat[,icol] <- tempF
            }
            myVar <- apply(myMat,1,sum,na.rm=T)
            myVar[which(myVar == 0)] <- NA # summing NA will do this, keeping the NA in place
          },
          N.pout = {myVar = N.pout_sumtab[[1]]$SSB}, # N.pout does not have catches data, just landings (and no discards)
          Herring = {myVar = Herring_sumtab[[1]]$SSB},
          Dab = {myVar = Dab_sumtab[[1]]$SSB},
          Whiting = {myVar = Whiting_sumtab[[1]]$SSB},
          Sole = {myVar = Sole_sumtab[[1]]$SSB},
          Gurnard = {myVar = Gurnard_sumtab[[1]]$SSB},
          Plaice = {myVar = Plaice_sumtab[[1]]$SSB},
          Haddock = {myVar = Haddock_sumtab[[1]]$SSB},
          Cod = {
            myMat <- matrix(NA, nrow = dim(Cod_sumtab[[1]])[1], ncol = 2)
            for(icol in 1:2)
            {
              tempF <- Cod_sumtab[[icol]]$SSB
              while(length(tempF)<dim(Cod_sumtab[[1]])[1]) tempF <- c(NA,tempF)
              myMat[,icol] <- tempF
            }
            myVar <- apply(myMat,1,sum,na.rm=T)
            myVar[which(myVar == 0)] <- NA # summing NA will do this, keeping the NA in place
          },
          Saithe = {myVar = Saithe_sumtab[[1]]$SSB},
          {}
  )
  while (length(myVar)<length(timePeriod)) myVar <- c(NA,myVar)
  SSBmat[,iSpecies] <- myVar
}
write.csv(SSBmat, file = "data/SSBmat.csv")
#time-averaged-SSB.csv
# averaged subset
SSBavg <- SSBmat[which(rownames(SSBmat) == "2014"):which(rownames(SSBmat) == "2019"),]
SSBavg <- apply(SSBavg,2,mean,na.rm=T)
SSBavg[is.nan(SSBavg)] <- NA
SSBavg <- data.frame("species" = SpIdx, "SSB_1419" = SSBavg,row.names = NULL)
write.csv(SSBavg, file = "data/time-averaged-SSB.csv", row.names = F)



#### Then apply these new files to the toyexample1.R (time-averaged calibration) and toyexample2.R (time-series evaluation/exploration)
#### Then apply these using Mike's code to do time-series fitting.





#' how to get more F per year/species
#' Fmat has data from 1983 up to 2010 for all species / should we use it?
#' comparing old and new code over 2000-2010 period
#'

newF <- fMatW[54:64,]
oldF <- Fmat[,44:54]
oldF <-t(oldF)
oldF <- cbind(2000:2010,oldF)

diff <- newF - oldF
diff$X <- 2000:2010

plot_dat <- reshape2::melt(diff,"X")

ggplot(plot_dat)+
  geom_line(aes(x = X, y = value, color = variable ))+
  ggtitle("newF - oldF")




# usethis::use_data(DATASET, overwrite = TRUE)
