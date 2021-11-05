
## This script should set up functions and data necessary for precipitation
## interpolation and forecasting.
## Sourcing other scripts ======================================================

source("R/basic_functions.R")

## Some base data ==============================================================

# stat <- readRDS("../commondata/stations_regionalized.RDS")
stat <- readRDS("data/precipitation_stations.RDS")

cent <- readRDS("data/cent980.RDS")
# 980 subbasins, includes the ones with the hydro stations as their own

# pcp00 <- readRDS("../weather/data/hm_00s/gwr_daily/gwr_lin/pmodeldata_lin_full.rds")
pcp00 <- readRDS("data/2000s_precipitation.RDS")
pcp00$index0 <- 1:nrow(pcp00)
# index1 corresponds to days

wetmodels <- readRDS("data/logistic_to_wetdry.RDS")
## models for IDW occurrence -> wet/dry probability

## Distance matrices
dm1 <- as.matrix(dist(cbind(stat$Longitude, stat$Latitude), 
                      diag = T, upper = T)) 
# stations to stations
allpoints <- rbind(cbind(cent$Long_, cent$Lat), cbind(stat$Longitude, stat$Latitude))
dm2 <- as.matrix(dist(allpoints))[1:nrow(cent), (nrow(cent)+1):nrow(allpoints)]
# stations to centroids

scalarsfile1 <- readRDS("data/scalars_combinedmodel.RDS")

day_indices <- unlist(sapply(2005:2020, function(y) {
  v1 <- sapply(1:ifelse(is_leap(y), 366, 365), function(n) {
    formatC(n, width=3, format="d", flag="0")
  })
  return(paste0(y, v1))
}))

## Functions ===================================================================

defaultlog <- "log.txt"
Log <- function(text, tofile=TRUE, filenam=defaultlog, ...) { # Log function 
  msg <- sprintf(paste0(as.character(Sys.time()), ": ", text), ...)
  cat(msg, "\n")
  if(tofile) {
    write(msg, defaultlog, append=TRUE)
  }
}

plot_obs <- function(index1, dat, col="obs") {
  dat2 <- dat[dat$index1 == index1, ]
  # column <- as.symbol(col)
  p1 <- ggplot(data = pan_shp, aes(x=long, y=lat)) + 
    geom_polygon(fill="grey80", aes(group=group)) +
    geom_point(data=dat2, aes(x=long, y=lat, color=get(col)), size=1.5)
  plot(p1)
}

########### Interpolation based on regions ###########
interpolate_reg <- function(par, pcp_out, pcp_avail, dis, pcpmodel="lin", 
                            alg, outtype=1, scalarsfile = scalarsfile1) {

  if(missing(pcp_avail) & missing(pcp_out)) {
    stop("need at least one precipitation data frame")
  } else if(missing(pcp_avail)) {
    pcp_avail <- pcp_out
  } else if(missing(pcp_out)) {
    pcp_out <- pcp_avail
  } ## If only one frame supplied, make it fit on itself
  
  fit_to_same <- identical(pcp_out, pcp_avail)
  
  Log(paste0("Parameters: ", paste0(as.character(par), collapse=" ")))
  interp <- vector(length = nrow(pcp_out))
  
  alpha <- par # parameter is decay of distance exponential
  
  ## Getting the right indices to fill in with foreach
  # Log("getting indices...")
  alldays = unique(pcp_out$index1)
  indices <- sapply(alldays, function(x) {
    return(which(pcp_out$index1 == x))
  })
  
  ndays <- length(alldays)
  failcount = 0
  
  if(ndays%%ncores != 0) {
    extra_rows <- ndays%%ncores
    # Log(paste0("extra rows for last core: ", extra_rows))
    i_subset <- floor(ndays/ncores)
    do_extra <- TRUE
  } else {
    i_subset <- ndays/ncores
    do_extra <- FALSE
  }
  i_subset0 <- i_subset
  
  interp[indices] <- foreach(n = 1:ncores, .combine='c') %dopar% {
    if(do_extra) {
      if(n == ncores) {
        i_subset <- i_subset + extra_rows
      }
    }
    
    unlist(lapply(1:i_subset, function(x) {
      i <- alldays[((n-1)*i_subset0 + x)]
      p_obs <- if(alg %in% c("wetfc", "wetprev")) {
        pcp_avail[which(pcp_avail$index1 == (i-1)), ]
      } else {
        pcp_avail[which(pcp_avail$index1 == i), ]
      }
      pcp_avail[which(pcp_avail$index1 == i), ]
      
      if(nrow(p_obs) == 0) return(NA)
      
      p_out <- pcp_out[which(pcp_out$index1 == i), ]
      # if(i %% 80 == 0) Log(paste0("core: ", n, " index: ", x/i_subset))
      # All precipitation records for day i
      out_tmp <- unlist(lapply(1:nrow(p_out), function(r) {
        # Kernel with values weighted by closeness to the point
        # Weight drops off with d as an exponential function
        # r1 <- ptmp$region[r]
        # alpha_r <- par[r1]
        
        stindex <- if(fit_to_same) {
          which(stat$Number == p_out$id[r]) 
        } else {
          p_out$sub[r]
        }
        
        pos <- if(fit_to_same) {
          # cat("a")
          which(p_obs$index0 != p_out$index0[r]) 
        } else {
          1:nrow(p_obs)
        }
        
        if(alg == "quantidw") {
          pos_wet <- if(fit_to_same) {
            which(p_obs$index0 != p_out$index0[r] & p_obs$wet_obs == 1)
          } else {
            which(p_obs$wet_obs == 1)}
        }
        
        n <- length(pos)
        if(n == 0){
          ## The no neighbor case, return a 'default' value 
          ## failcount tracks total number of such cases (should be negligible)
          failcount = failcount + 1
          return(switch(alg, 
                        # "eps"=ptmp$epsbar[r], 
                        # "eps2"=ptmp$epsbar[r],
                        "eps3"=0, "wet"=0.5, "wet2"=0.5, "wet3"=0.5,
                        "quantidw"=0.5))
        }
        
        # ids <- switch(alg, "wet"=p_obs$id[pos], "quantidw"=p_obs$id[pos_wet])
        ids <- p_obs$id[pos]
        
        # if(test3) {
        #   # testing regional averages, 3 and 4 are w/d and quantity
        #   v_interp <- as.numeric(runif(1) < mean(p_obs$wet_obs[pos]))
        #   return(v_interp)
        # } else if (test4) {
        #   v_interp <- mean(p_obs$obs[pos])
        #   return(v_interp)
        # }
        
        values <- switch(alg, ## Values to interpolate with
                         "quantidw"=p_obs$obs[pos_wet],
                         "wet"=p_obs$wet_obs[pos])
        
        d <- unlist(lapply(1:n, function(x) {
          return(dis[stindex, which(stat$Number == ids[x])])
        }))
        
        if(alg == "wet") {
          v_interp <- (sum(values*exp(-alpha * d)) / (sum(exp(-alpha * d)))) 
        } else if(alg == "quantidw") {
          dnum <- unlist(lapply(1:length(pos_wet), function(x) {
            return(dis[stindex, which(stat$Number == ids[x])])
          }))
          v_interp <- (sum(values*exp(-alpha*dnum)) / sum(exp(-alpha * d)))
        }
        
        return(v_interp)
      }))
      return(out_tmp)
      # e_interp[ptmp$index0] <- out_tmp
    }))
  }
  
  pcp2 <- pcp_out
  
  if(alg == "quantidw") {
    
    pcp2$predidw <- interp
    
    if(fit_to_same) {
      sumofsq <- sum((pcp2$obs - pcp2$predidw)^2)
      Log(paste0("Sum of squares: ", sumofsq))
      Log(paste0("NSE: ", NSE(pcp2$predidw, pcp2$obs), 
                 ", Region: ", unique(pcp2$region)))
      optmetric <- sumofsq
    }
    
  } else if(alg == "wet") {
    pcp2$w_interp2 <- interp
    if(fit_to_same) {
      if(any(is.na(pcp2$w_interp2))) {
        pcp2 <- pcp2[-which(is.na(pcp2$w_interp2)), ]
      }
      fit1 <- glm(wet_obs ~ w_interp2, pcp2, family="binomial")
      probs <- predict(fit1, pcp2, type='response')
      error_rate <- mean(unlist(lapply(1:100, function(x) {
        rng <- runif(length(probs)) 
        preds <- rng < probs
        return(mean(preds != pcp2$wet_obs))
      })))
      Log(paste0("Error rate: ", error_rate))
      optmetric <- error_rate
    } else {
      cat(pcp2$region[1], pcp2$mon[1], "\n")
      wetmodel1 <- wetmodels[[pcp2$region[1]]][[pcp2$mon[1]]] ### make less janky
      pcp2$w_interp3 <- predict(wetmodel1, pcp2, type='response')
    }
  } 
  
  # Log(paste0("Failcount: ", failcount))
  return(switch(outtype, optmetric, pcp2, fit1))
}


########### Interpolate all 6 regions ###########
all_interp <- function(pcpraw, pcpout, step, parset=0, nregions=6, algorithm, ...) {
  
  all_interp_mainfunction <- function(step, algorithm, ...) {
    if(step == "fit") {
      opttmp <- optimize(f=interpolate_reg, pcp_avail=pcp0, dis=dm1,
                         interval=c(0.001, 30), alg=algorithm, outtype=1, ...)
      return(opttmp$minimum)
    } else if (step == "gen") {
      par0 <- parset[[which(reg_mon$region == region & reg_mon$month == month)]]
      return(interpolate_reg(par=par0, pcp_avail=pcp0, pcp_out=pcp2, outtype=2, 
                             alg=algorithm, ...))
    } else if (step == "getcoeff") {
      par0 <- parset[[region]]
      return(interpolate_reg(par0, pcp_avail=pcp0, pcp_out=pcp2, outtype=3, 
                             alg=algorithm, ...))
    }
  }
  
  if(missing(pcpout)) pcpout <- pcpraw
  reg_mon <- expand.grid(1:6, 1:12); names(reg_mon) <- c("region", "month")
  
  for(region in 1:nregions) {
    for(month in 1:12) {
      assign(paste0("pcp_", month, region), 
             pcpraw[which(pcpraw$region == region & pcpraw$mon == month), ])
      assign(paste0("pcpo_", month, region), 
             pcpout[which(pcpout$region == region & pcpout$mon == month), ])
    }
  }
  
  out <- vector("list", length=nrow(reg_mon))
  for(region in 1:nregions) {
    for(month in 1:12) {
      cat("Region:", region, "Month:", month, "\n")
      rmindex <- which(reg_mon$region==region & reg_mon$month==month)
      pcp0 <- get(paste0("pcp_", month, region))
      # if(test2) {
      #   # remove what wasnt there originally
      #   pos0 <- unlist(lapply(1:nrow(pcp0), function(x) {
      #     if(length(intersect(which(pcp00$id == pcp0$id[x]),
      #                         which(pcp00$index1 == pcp0$index1[x]))) == 0) {
      #       return(x)
      #     }}))
      #   if(length(pos0) > 0) pcp0 <- pcp0[-pos0,]
      # }
      pcp2 <- get(paste0("pcpo_", month, region))
      # 
      out[[rmindex]] <- all_interp_mainfunction(step, algorithm, ...)
    }
  }
  if(step=="gen") {
    return(do.call(rbind, out))
  }
  return(out)
}

########### Preparing output df and retrieving GWR values ###########
build_outdf <- function(nxy=980, nday=4017, startyr=2005, pmeans) {
  # pmeans can be "wclim", "wclim50", or "gwr"
  
  nyear <- round(nday/365.25)
  leaps <- unlist(lapply(0:(nyear-1), function(x) {
    is_leap(startyr + x)
  }))
  
  mon0 <- unlist(lapply(leaps, function(leapx) {
    if(leapx) {
      return(day_months_leap)
    } else {
      return(day_months)
    }
  }))
  
  out <- data.frame(sub = unlist(lapply(1:nxy, function(x) rep(x, nday))),
                    index1 = rep(1:nday, nxy),
                    mon = rep(mon0[1:nday], nxy))
  
  p0 <- data.frame(sub = rep(1:nxy, 12), 
                   mon = unlist(lapply(1:12, function(x) {
                     rep(x, nxy)})))
  
  p0$pred0 <- switch(pmeans, 
                     "default" = unlist(lapply(1:12, function(x) {
                       return(NA)
                     })),
                     "gwr" = unlist(lapply(1:12, function(x) {
                       gwrdf[[x]]$pred
                     })),
                     "wclim" = unlist(lapply(1:12, function(x) {
                       cat("Month", x, "\n")
                       longlat <- cbind(gwrdf[[x]]$Long_, gwrdf[[x]]$Lat)
                       return(extract(wcpresent[[x]], longlat) / month_days2[x])
                     })),
                     "wclim50" = unlist(lapply(1:12, function(x) {
                       cat("Month", x, "\n")
                       longlat <- cbind(gwrdf[[x]]$Long_, gwrdf[[x]]$Lat)
                       return(extract(wcfuture[[x]], longlat) / month_days2[x])
                     })))
  
  out <- merge(out, p0)
  out$region <- cent$region[out$sub]
  out <- out[order(c(out$sub, out$mon)), ]
  
  if(any(is.na(out$sub))) {
    out <- out[-which(is.na(out$sub)), ]
  }
  
  return(out)
}

########### Precipitation generators & SAC #############

generate_weights <- function(coords, alpha) {
  W <- 1/(exp(alpha * as.matrix(dist(coords, diag=T, upper=T))))
  diag(W) = 0
  W <- W/rowSums(W) #normalizing
  return(W)
}

gamma_fit <- function(x, region, month) {
  dat <- x[[region, month]]
  
  if(any(is.na(dat))) {
    dat <- dat[-which(is.na(dat))]
  }
  
  offset <- 0
  if(min(dat) < 0) {
    offset <- -min(dat)
    dat <- dat + offset
  }
  
  fit1 <- fitdist(dat, distr="gamma", method="mme")
  plot(fit1)
  return(list(offset, fit1))
}

# generate_V <- function(gamma, coords, alpha, u) {
#   # u is a vector of random numbers (0, 1) same length as coords
#   n <- if(length(dim(coords)) > 1) {
#     nrow(coords)
#   } else {
#     length(coords)
#   }
#   if(missing(u)) {
#     u <- runif(n) # random numbers for each station
#   }
# 
#   W <- generate_weights(coords, alpha)
# 
#   # cdf <- unlist(lapply(1:10000, function(x) {
#   # return(gamma * W %*% u + u)
#   # }))
# 
#   # V <- gamma * W %*% u + u # %*% is matrix multiplication
# 
#   test1 <- sort(unlist(lapply(1:5000, function(x) {
#     u <- runif(n)
#     gamma * W %*% as.vector(u) + u
#   })))
#   if(any(test1 < 0)) {
#     testtmp <- test1 - min(test1)
#     test1_scaled <- testtmp/max(testtmp)
#   } else {
#     test1_scaled <- test1/max(test1)
#   }
# 
#   V0 <- gamma * W %*% u + u
#   V <- unlist(lapply(V0, function(x) {
#     pos<-which.min(abs(test1-x))
#     return(test1_scaled[pos])
#   }))
#   # V <- gamma * W %*% u + (1 - gamma) * u
#   return(V)
# }

ac_day <- function(day, df, maxdist=100, var="obs", metric="moran", verbose = F,
                   zeroes = F, nscale = T)  {
  if(verbose) cat("Day:", day, "\n")
  if(!(day %in% df$index1)) return(NA)
  
  ptmp <- df[df$index1 == day, ]
  
  if(zeroes == F & var == "obs") {
    ptmp <- ptmp[which(ptmp$obs > 0), ]
    if(nrow(ptmp) < 4) {
      return(NA)
    }
  }
  ### Plotting semivariogram
  # p.vgm <- variogram(pcp ~ 1, ptmp)
  # p.fit <- fit.variogram(p.vgm, model = vgm("Sph"))
  # plot(p.vgm, p.fit, main = day)
  
  n <- nrow(ptmp)
  coords <- cbind(ptmp$long, ptmp$lat)
  W <- 1/as.matrix(dist(coords, diag=T, upper=T)); diag(W) = 0
  W <- W/rowSums(W)
  if(any(is.na(W) | is.infinite(W))) {
    W[is.na(W) | is.infinite(W)] <- 0
  }
  wlist <- mat2listw(W)
  
  out <- switch(metric,
                "geary" = tryCatch({
                  gearyc <- geary.test(ptmp[[var]], wlist, zero.policy = FALSE)
                  return(gearyc$estimate[1])
                }, warning = function(war) {
                  print(paste("Warning:", war))
                  return(gearyc$estimate[1])
                }, error = function(err) {
                  print(paste("Error:", err))
                  return(NA)
                }),
                "moran" = moran(ptmp[[var]], wlist, n, n)$I)
  
  if(nscale) out <- out + 1/(n-1)
  
  return(out)
}

####### Writing pcp in ##############

postprocess_pcp <- function(p, predname="predidw") {
  # Assumes you've already run the interpolation for wet/dry and quantity
  # need to be named wet_obs
  pred <- p[, predname]
  pred[which(is.na(pred))] <- 0
  pred[p$wet_obs == 0] <- 0
  if(any(pred < 0.5 & p$wet_obs == 1)) {
    pred[pred < 0.5 & p$wet_obs == 1] <- 0.5
  }
  p$pcp <- pred
  return(p)
}

postprocess_and_scale <- function(p_out, scale) {
  # scale can be 'scalars', 'redist', or 'n'
  # test1 = F
  
  scalars <- scalarsfile1
  if(scale=='redist') {
    means <- tapply(p_out$predidw, list(p_out$mon, p_out$sub), function(x) {
      mean(x, na.rm=T) })
  }
  
  p_out$wet_obs <- as.numeric(runif(nrow(p_out)) < p_out$w_interp3)
  p_out$pcp <- p_out$predidw
  
  # if(!test1){
  #   p_out$pcp[p_out$wet_obs == 0] <- 0
  #   p_out$pcp[p_out$wet_obs == 1 & p_out$pcp < 0.5] <- 0.5
  # }
  p_out$pcp[is.na(p_out$pcp)] <- 0
  
  p_scaled <- p_out
  if(scale=='redist') {
    newmeans <- tapply(p_out$pcp, list(p_out$mon, p_out$sub), function(x) {
      mean(x, na.rm=T)
    })
    p_scaled$pcp <- unlist(lapply(1:nrow(p_out), function(x) {
      if(x %% 1000 == 0) cat(x, "\n")
      i = p_out$mon[x]; j = p_out$region[x]
      return(p_out$pcp[x] * means[i,j]/newmeans[i,j])
    }))
  } else if(scale=='n') {
    
  } else if(scale=='scalars'){
    p_scaled$pcp <- unlist(lapply(1:nrow(p_out), function(x) {
      if(x %% 1000 == 0) cat(x, "\n")
      return(p_out$pcp[x] * scalars[p_out$mon[x], p_out$region[x]])
    }))
  }
  
  p_scaled2 <- write_swatpcp(dat=p_scaled, output_type=".pcp")
  return(p_scaled2)
}

## MAKE SURE THERE IS A BACKUP FIRST
write_swatpcp <- function(dat, output_type='.pcp', model_type='interpolated',
                          sim="data", 
                          startdate="20050101", cent=cent) {
  if(output_type == '.pcp' & model_type == 'interpolated') {
    pcp1 <- readLines("data/pcpnew.pcp")
    sampleline <- pcp1[[5]]
    pcpnew <- pcp1[1:4]
    
    dat <- dat[order(dat$sub, dat$index1), ] # IMPORTANT
    ndays <- max(dat$index1)
    
    for(d in 1:ndays) {
      
      cat(d, "/", ndays, "\n")  
      
      outd <- dat[dat$index1 == d, ]
      outd$pcp <- round(outd$pcp, 1)
      outd$pcp <- as.character(outd$pcp)
      pcpin <- outd$pcp
      # pcpin <- outd$pcp[nbr] ## IMPORTANT
      
      ## 915 forecasting seed stations but 980 stations in the model
      
      pcpin <- unlist(lapply(pcpin, function(s) {
        if (length(grep("\\.", s)) == 0) 
          s <- paste0(s, ".0")
        L <- nchar(s)
        if (L < 5)  
          s <- paste0(strrep("0", 5-L), s)
        return(s)
      }))
      
      ### substr 1:8 is the date
      newline <- paste0(day_indices[d], paste0(pcpin, collapse=''))
      pcpnew[[d+4]] <- newline
    }
    return(pcpnew)
  }
  if(output_type == 'p-files') {
    # Works more reliably than creating new pcp1.pcp for new watershed delination
    # For interpolated data
    if(model_type == 'interpolated') {
      # Input interpolated data
      if(length(grep("pfiles", getwd())) == 0) setwd("./pfiles/")
      
      # names(dat)[2] <- 'index1'
      ndays <- max(dat$index1)
      
      for(i in 1:max(dat$sub)) {
        print(paste0("Working on station: ", i))
        dat1 <- dat[dat$sub == i, ]
        pcpi <- unlist(lapply(0:ndays, function(d) {
          if(d == 0) return(startdate)
          else {
            val <- dat1$pcp[dat1$index1 == d]
            return(as.character(round(val, 1)))
          }
        }))
        writeLines(pcpi, paste0("p", i, ".txt"))
      }
      
      pcps <- data.frame(ID = 1:max(dat$sub),
                         NAME = paste0("p", 1:max(dat$sub)),
                         LAT = cent$Lat - 0.001,
                         LONG = cent$Long_,
                         ELEVATION = cent$Elev)
      write.csv(pcps, "pcp.txt", quote = F, row.names = F)
      
    } else if(model_type == 'observed') {
      # For 'null' model, i.e. measured data
      
      dir.create("./pfiles_null")
      setwd("./pfiles_null/")
      
      startdate <- "20050101"
      statns <- unique(dat$id)
      
      iterations <- length(statns)
      
      if(iterations%%ncores != 0) {
        extra_rows <- iterations%%ncores
        Log(paste0("extra rows: ", extra_rows))
        subset <- floor(iterations/ncores)
        subset0 <- subset
        do_extra <- TRUE
      } else {
        subset <- iterations/ncores  
        Log(paste("iterations:", iterations, "per core:", subset))
        do_extra <- FALSE
      }
      
      df0 <- foreach(n = 1:ncores) %dopar% {
        
        if(do_extra & n == ncores) {
          subset <- subset + extra_rows
        }
        
        l1 <- lapply(1:subset, function(x) {
          i <- statns[(n-1)*subset0 + x]
          if(x%%5 == 0) {
            Log((paste("Core:", n, "Progress:", x/subset)))
          }
          pcpi <- unlist(lapply(0:4017, function(d) {
            if(d == 0) return(startdate)
            else {
              val <- dat$obs[dat$id == i & dat$index1 == d]
              if(length(val) == 0) return("-99.0") # SWAT will interpolate -99s
              return(as.character(round(val, 1)))
            }
          }))
        })
      }
      
      nstn = 0
      for(i in 1:length(df0)) {
        df00 <- df0[[i]]
        for(j in 1:length(df00)) {
          nstn <- nstn + 1
          cat(nstn, '\n')
          writeLines(df00[[j]], paste0("p", nstn, ".txt"))
        }
      }
      stations <- readRDS("../weather/data/stations.rds")
      #### NEED A DIFFERENT OPTION FOR CENTROIDS
      pcps <- data.frame(ID = 1:n,
                         NAME = paste0("p", 1:n),
                         LAT = unlist(lapply(statns, function(s) {
                           return(stations$Latitude[stations$Number == s] - 0.001)
                         })),
                         LONG = unlist(lapply(statns, function(s) {
                           return(stations$Longitude[stations$Number == s])
                         })),
                         ELEVATION = unlist(lapply(statns, function(s) {
                           return(stations$Elevation.m[stations$Number == s])
                         })))
      write.csv(pcps, "pcp.txt", quote = F, row.names = F)
    }
  }
}

interpolate_full <- function(pcp, pred1, pmeans='default', 
                             pars_wetdryfile = "data/optim_alpha_wet_month.RDS",
                             pars_quantfile = "data/opt_qidw_elev.RDS",
                             scalarsfile0,
                             simname = "~/panama_swat3/Scenarios/sim1_p2step/",
                             pname,
                             tofile = TRUE,
                             distmat = dm2,
                             nday = 5843,
                             ...) {
  
  pars_wetdry <- readRDS(pars_wetdryfile)
  pars_quant <- readRDS(pars_quantfile)
  
  cat("Building data frame...", "\n")
  pcp_out00 <- build_outdf(pmeans=pmeans, nday=nday)
  
  pcp$wet_obs <- as.numeric(pcp$obs >= 0.5)
  if(any(is.na(pcp$wet_obs))) {
    posna <- which(is.na(pcp$wet_obs))
    pcp$wet_obs[posna] <- sample(0:1, length(posna), replace=T)
  }
  
  obspos <- which(colnames(pcp) == pred1)
  names(pcp)[obspos] <- "obs" ## need to do this to pass correctly to quantidw
  
  cat("Calculating wet days...", "\n")
  pcp_outwd <- all_interp(pcpraw=pcp, pcpout=pcp_out00, step="gen", dis=distmat,
                          parset=pars_wetdry, algorithm="wet")
  pcp_outwd$wet_obs <- unlist(lapply(pcp_outwd$w_interp2, function(p) {
    return(as.numeric(runif(1) < p))
  }))
  cat("Calculating quantities...", "\n")
  pcp_outq <- all_interp(pcpraw=pcp, pcpout=pcp_outwd, step="gen", dis=distmat,
                         parset=pars_quant, algorithm="quantidw")
  
  if(tofile) {
    cat("Writing to file...", "\n")
    pcpfile <- write_swatpcp(pcp_outq, sim=simname, ...)
    outpath <- paste0(simname, "TxtInOut/", pname)
    writeLines(pcpfile, outpath)
  } else {
    return(pcp_outq)
  }
}

####### Writing temp in ##############

write_temp <- function(type='future', rcp, model, year,
                       scen="~/panama_swat3/Scenarios/sim1_p2step/TxtInOut/") {
  oldwd <- getwd()
  setwd(scen)
  
  
  if(type == 'future') {
    tmin <- raster::getData('CMIP5', var='tmin', rcp=rcp, model=model, year=year,
                            res=2.5)
    tmax <- raster::getData('CMIP5', var='tmax', rcp=rcp, model=model, year=year,
                            res=2.5)
  } else {
    # wc present rasters
    tmin <- raster::getData('worldclim', var='tmin', res=0.5, lon=-80, lat=8)
    tmax <- raster::getData('worldclim', var='tmax', res=0.5, lon=-80, lat=8)
  }
  
  longlat <- cbind(wgen915$WLONGITUDE, wgen915$WLATITUDE)
  
  tmins <- lapply(1:12, function(m) {
    extract(tmin[[m]], longlat)
  })
  tmins <- do.call(cbind, tmins)
  tmins <- tmins/10
  tmins <- format(round(tmins, 2), nsmall=2)
  
  tmaxs <- lapply(1:12, function(m) {
    extract(tmax[[m]], longlat)
  })
  tmaxs <- do.call(cbind, tmaxs)
  tmaxs <- tmaxs/10
  tmaxs <- format(round(tmaxs, 2), nsmall=2)
  
  wgn <- list.files()[grep(".wgn", list.files())]
  
  for(w in wgn) {
    cat("Editing", w, "...", "\n")
    wtmp <- readLines(w)
    
    w1 <- wtmp[[1]]
    pos1 <- regexpr("NAME:", w1)[1] + 5
    pos2 <- pos1 + regexpr(" ", substr(w1, pos1, stop=1000)) - 1
    stationid <- as.numeric(substr(w1, pos1, pos2))
    
    wtmp[5] <- paste0(" ", paste(tmaxs[stationid,], collapse=" "))
    wtmp[6] <- paste0(" ", paste(tmins[stationid,], collapse=" "))
    
    writeLines(wtmp, w)
  }
  
  #cleanup
  setwd(oldwd)
}


## writing null
# setwd("../weather/data/hm_00s/")
# dat <- readRDS("pcpdaily.rds")
# names(dat)[10] <- "obs"
# ncores=15
# registerDoParallel(cores = ncores)
# write_swatpcp(dat=dat, output_type="p-files", model_type="observed")

# 980 subbasin centroids =

# shp1 <- shapefile("../commondata/full_watershed/wshednew.shp")
# regions <- readRDS("../commondata/regionalized_watershed.rds")
# 
# cent980 <- as.data.frame(shp1)
# coordinates(cent980) <- ~Long_ + Lat
# crs(cent980) <- crs(regions)
# 
# cent980$region <- over(cent980, regions)$region
# saveRDS(cent980, "../commondata/cent980.RDS")

## Preparing distance matrix & neighbor lists =
# 
# dm1 <- as.matrix(dist(rbind(cbind(cent$Long_, cent$Lat), 
#                             cbind(stat$Longitude, stat$Latitude)), 
#                       method = "euclidean"))
# dm1 <- dm1[1:nrow(cent), (nrow(cent)+1):(nrow(cent)+nrow(stat))]
# 
# idmap <- data.frame(id1 = unique(p$id),
#                     id2 = unique(as.numeric(as.factor(p$id))))
# p$id <- as.numeric(as.factor(p$id))
# 
# if(interp_type == 1) {
#   neighbors1 <- lapply(1:nrow(cent), function(x) {
#     neighb_indices <- which(dm1[x, ] < bw)
#     return(idmap$id2[idmap$id1 %in% stat$Number[neighb_indices]])
#   })
# }

## REgionalizing observed data ===
# reg <- readRDS("../../../../commondata/stations_regionalized.RDS")
# pcpraw$region <- unlist(lapply(1:nrow(pcpraw), function(x) {
#   reg$region[which(reg$Number == pcpraw$id[x])]
# }))

## Switches ==

# nregions <- 6 # threshold if 0?
# algorithm <- "wet" # choices are "eps" "eps2" "eps3" "wet" "wet2"
# parset <- readRDS("p_wet_fit.RDS") # a list with parameter values by region
# postprocess <- FALSE
# outname <- "p_wet_fit_coeffs.RDS"
# 
# step = "getcoeff"
# if (step == "fit") {
#   par_start <- c(0.01, 30)
# }
# "fit" to fit using optim or optimize
# "gen" to generate interpolated data using parameters
