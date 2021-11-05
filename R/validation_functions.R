library(raster)
library(hydroGOF)
library(reshape2)
library(RColorBrewer)
library(leaps)
source("R/basic_functions.R")

### Data =======================================================================
r0 <- readRDS("data/results/obs_hydro.RDS")

outliers <- c(11, 564, 255, 190, 403, 393, 530)
# 11 is on the border with costa rice
### 564 just seems poorly delineateed- assigned to the wrong stream
### 403 is in the middle of indigenous territory - bad pcp
### 255 only has two points of observation?
# 94 is downstream of the Bayano dam -- FIXED, it's sub 981 now
### 190 is affected by lake Gatun, the station is downstream of a dam
# 530 also on the border with costa rica
## Three hashes means they're the most downstream observation pt of that river
## Explanation for those: border, poorly delineated, downstream of dam

wshed <- shapefile("data/wshed_regionalized.shp")
upst <- readRDS("data/upstream.rds")
downst <- readRDS("data/downstream.rds")

alphas_x <- readRDS("data/optim_alpha_wet_month.RDS")
# interpolation kernel decay parameter values for wet/dry
alphas_q <- readRDS("data/opt_qidw_elev.RDS") 
# interpolation kernel decay parameter values for quantity

met_st <- readRDS("data/met_stations.RDS") # meteo. stations locations
hyd_st <- readRDS("data/hyd_stations.RDS") # hydro. stations locations

pcp90 <- readRDS("data/all_days_records.RDS")

### Functions ==================================================================


# Read reach output file stored at filepath
# full = T for general use, full = F for nutrients - only sampled time periods
read_rch <- function(filepath, full = TRUE, daily = FALSE,
                     cols = "default",
                     antebayanohack=TRUE,
                     simyears = 2005:2020) {
  # filepath = paste("~/panama_swat/Scenarios/", scen, "/TxtInOut", sep = '')
  # setwd(filepath)
  print(paste("Reading output.rch..."))
  
  if(cols == "default") {
    rch <- read.fortran(filepath, skip = 9,
                        c("6X", "I4", "X9", "F6", "2X12", "F12" ,"44X12"))
    names(rch) <- c("SUB", "MON", "FLOW_OUT")
  } else {
    rch <- read.fortran(filepath, skip = 9,
                        c("6X", "I4", "X9", "F6", "47F12"))
    names(rch) <- c("SUB", "MON", "AREA", "FLOW_IN", "FLOW_OUT", "EVAP", "TLOSS",
                    "SED_IN", "SED_OUT", "SEDCONC", "ORGN_IN", "ORGN_OUT", 
                    "ORGP_IN", "ORGP_OUT", "NO3_IN", "NO3_OUT", "NH4_IN", 
                    "NH4_OUT", "NO2_IN", "NO2_OUT", "MINP_IN", "MINP_OUT", 
                    "CHLA_IN", "CHLA_OUT", "CBOD_IN", "CBOD_OUT", "DISOX_IN",
                    "SOLPST_IN", "SOLPST_OUT", "SORPST_IN", "SORPST_OUT", 
                    "REACTPST", "VOLPST", "SETTLPST", "RESUSP_PST", "DIFFUSEPST", 
                    "REACBEDPST", "BURYPST", "BED_PST", "BACTP_OUT", "BACTLP_OUT",
                    "CMETAL1", "CMETAL2", "CMETAL3", "TOTN", "TOTP", "NO3CONC")
    rch <- rch[, cols]
  }
  
  print("Processing...")
  
  if(antebayanohack) {
    # delineation problem, we will be merging the outputs of two subs 95 & 97
    # into a virtual "SUB 981" for the ante bayano station
    Qantebayano <- rch$FLOW_OUT[rch$SUB==95] + rch$FLOW_OUT[rch$SUB==97]
    
    nsub <- length(unique(rch$SUB))
    nsubnew <- nsub + 1
    nday <- nrow(rch)/nsub
    rch3 <- data.frame(SUB = rep(1:nsubnew, nday))
    indicesall <- 1:nrow(rch3)
    indicesold <- indicesall[-which(indicesall%%nsubnew==0)]
    indicesnew <- which(indicesall%%nsubnew==0)
    rch3$MON <- NA
    rch3$MON[indicesold] <- rch$MON
    rch3$FLOW_OUT <- NA
    rch3$FLOW_OUT[indicesold] <- rch$FLOW_OUT
    rch3$MON[indicesnew] <- rch3$MON[indicesnew-1]
    rch3$FLOW_OUT[indicesnew] <- Qantebayano
    rch <- rch3 ## imputing virtual subbasin 981 values to the dataframe
  }
  
  ## Adding year and month
  nsub <- length(unique(rch$SUB))
  yeardays <- ifelse(is_leap(simyears), 366, 365)
  rch$YEAR <- rep(simyears, yeardays*nsub)[1:nrow(rch)]
  names(rch)[2] <- "DAY"
  rch <- rch[order(rch$SUB, rch$YEAR, rch$DAY),]
  rch$MON <- rep(unlist(sapply(simyears, function(y){
    if(is_leap(y)) {
      out <- rep(1:12, month_days_leap)
    } else {
      out <- rep(1:12, month_days)
    }
    
    if(y == max(simyears)) {
      end <- length(rch$YEAR[rch$YEAR == y])/981
      return(out[1:end])
    } 
    return(out)
  })), nsub)
  
  return(rch)
}

# Make vector if not
vectorize <- function(x) {
  if(!is.vector(x)) {
    x <- as.vector(x)
  }
  return(x)
}

# R^2 of linear model
rsq <- function(x, y) {
  x <- vectorize(x)
  y <- vectorize(y)
  summary(lm(x ~ y))$r.squared
}




match_day <- function(rch) {
  obs <- r0
  obs <- obs[-which(obs$sub %in% outliers), ]
  obs <- obs[which(obs$year %in% simyears), ]
  
  rch <- rch[which(rch$SUB %in% unique(obs$sub)), ]
  nsub <- length(unique(rch$SUB))
  
  qpred <- rch$FLOW_OUT
  names(qpred) <- paste0(rch$SUB, rch$YEAR, rch$DAY)
  obs$name <- paste0(obs$sub, obs$year, obs$day0)
  obs$pred <- qpred[obs$name]
  results <- data.frame(sub = obs$sub, mon = obs$mon, day = obs$day0,
                        year = obs$year, Qs = obs$pred, Qo = obs$q, 
                        res = obs$q - obs$pred)
  
  return(results)
}

validate <- function(v, nfl = 3, by = 'subm', type='mean', ...) { 
  
  ## nfl is number of maxima for flood analysis
  ## valid 'by': 'sub', 'subm', 'suby'
  if(any(v$sub %in% outliers)) {
    v <- v[-which(v$sub %in% outliers),]
  }
  
  v <- v[which(v$sub %in% names(upst)), ]
  # only most downstream point of each watershed used in analysis
  
  tapplygrp <- switch(by, 'sub' = v$sub, 'subm' = list(v$sub, v$mon),
                      'submy' = list(v$sub, v$mon, v$year))
  
  qobs <- tapply(v$Qo, tapplygrp, function(x) {
    if (type == 'mean') {
      return(mean(x, na.rm=T))
    } else if (type == 'sd') {
      return(sd(x, na.rm=T))
    } else if (type == 'fl') {
      return(sort(x, decreasing=T)[1:nfl])
    }
  })
  
  qsim <- tapply(v$Qs, tapplygrp, function(x)  {
    if (type == 'mean') {
      return(mean(x, na.rm=T))
    } else if (type == 'sd') {
      return(sd(x, na.rm=T))
    } else if (type == 'fl') {
      return(sort(x, decreasing=T)[1:nfl])
    }
  })
  
  if(type == 'fl') {
    qo <- unlist(qobs); qs <- unlist(qsim)
    sub <- rep(sapply(row.names(qobs), function(n) rep(n, nfl)), 12)
    mon <- unlist(lapply(1:12, function(n) rep(n, nfl*nrow(qobs))))
    q2 <- data.frame(sub = sub, mon = mon, qo = qo, qs = qs)
  } else {
    qo <- setNames(melt(qobs), c('sub', 'mon', 'qo'))
    qs <- setNames(melt(qsim), c('sub', 'mon', 'qs'))
    q2 <- merge(qo, qs)
  }
  
  
  print(rsq(q2$qo, q2$qs))
  print(NSE(q2$qs, q2$qo))
  print(paste("Pbias", pbias(q2$qs, q2$qo)))
  
  plot(q2$qs, q2$qo, pch=19, ...)
  abline(0, 1)
  return(q2)
  
}








