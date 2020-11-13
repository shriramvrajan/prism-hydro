### Validation of mean, sd, extremes

# source("R/general_stuff.R")

### Data =======================================================================

r0 <- readRDS("results/obs_hydro.RDS")
r1 <- readRDS("results/SWATnaive_model.RDS")
r2 <- readRDS("results/SWATsampling_model.RDS")
r3 <- readRDS("results/SWATfinal_model.RDS")

# r4 <- readRDS("../panama_scripts/discharge/results/")

outliers <- c(11, 564, 255, 190, 403, 393, 530)
# 11 is on the border with costa rice
# 564 just seems poorly delineateed- assigned to the wrong stream
# 403 is in the middle of indigenous territory - bad pcp
# 255 only has two points of observation?
# 94 is downstream of the Bayano dam -- FIXED, it's sub 981 now
# 190 is affected by lake Gatun, the station is downstream of a dam
# 393 also on the border with costa rica
# 530 also on the border with costa rica

wshed <- shapefile("data/wshed_regionalized.shp")
upst <- readRDS("data/upstream.rds")

alphas_x <- readRDS("data/optim_alpha_wet_month.RDS")
# alpha values for wet/dry interpolation
alphas_q <- readRDS("data/opt_qidw_noelev.RDS") 
# alpha values for quantity

pcp90 <- readRDS("~/panama_scripts/weather/data/hm_90s/interp2/pmodel_winterp.RDS") 

## Functions ===================================================================

match_day <- function(rch, simyears=2006:2015) {
  obs <- r0
  obs <- obs[-which(obs$sub %in% outliers), ]
  obs <- obs[which(obs$year %in% simyears), ]
  
  rch <- rch[which(rch$SUB %in% unique(obs$sub)), ]
  nsub <- length(unique(rch$SUB))
  
  qpred <- rch$FLOW_OUT
  names(qpred) <- paste0(rch$SUB, rch$YEAR, rch$DAY)
  obs$name <- paste0(obs$sub, obs$year, obs$day0)
  obs$pred <- qpred[obs$name]
  results <- data.frame(sub = obs$sub, mon = obs$mon, 
                        year = obs$year, Qs = obs$pred, Qo = obs$q, 
                        res = obs$q - obs$pred)
  
  return(results)
}

validate <- function(v, nfl = 1, by = 'subm', type='mean', ...) { 
  
  ## nfl is number of maxima for flood analysis
  ## valid 'by'   : 'sub', 'subm', 'suby'
  ## valid 'type' : 'mean', 'sd', 'fl'
  
  if(any(v$sub %in% outliers)) {
    v <- v[-which(v$sub %in% outliers),]
  }
  
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
    qobs <- unlist(as.list(as.data.frame(qobs)))
    qsim <- unlist(as.list(as.data.frame(qsim)))
  }
  
  print(rsq(as.vector(qobs), as.vector(qsim)))
  print(NSE(as.vector(qsim), as.vector(qobs)))
  plot(qobs, qsim, pch=19, ...)
  abline(0,1)
  
  return(data.frame(qo = qobs, qs = qsim))
}


### Mean SD max ================================================================
md1 <- match_day(r3)
v1 <- validate(md1)

md1 <- match_day(r1)
md2 <- match_day(r2)
md3 <- match_day(r3)

{ # means
  par(mfrow = c(1,3))
  v1 <- validate(md1, cex=0.7, main='Mean discharge,\n NN model')
  v2 <- validate(md2, cex=0.7, main='Mean discharge,\n IWGEN model')
  v3 <- validate(md3, cex=0.7, main='Mean discharge,\n RDW model')
}

{ #sd and extreme
  par(mfrow = c(2,3))
  v1s <- validate(md1, cex=0.7, type='sd', main='SD discharge,\n NN model')
  v2s <- validate(md2, cex=0.7, type='sd', main='SD discharge,\n IWGEN model')
  v3s <- validate(md3, cex=0.7, type='sd', main='SD discharge,\n RDW model')
  
  v1e <- validate(md1, cex=0.7, type='fl', 
                  main='Max. discharge,\n NN model')
  v2e <- validate(md2, cex=0.7, type='fl', 
                  main='Max. discharge,\n IWGEN model')
  v3e <- validate(md3, cex=0.7, type='fl', 
                  main='Max. discharge,\n RDW model')
}

### Drought ====================================================================
{ # drought
  
  rankyears <- function(v=v1, type='s') {
    if(any(v$sub %in% outliers)) v <- v[-which(v$sub %in% outliers), ]
    
    if(type == 'o') {
      sums <- tapply(v$Qo, list(v$sub, v$year), function(x) sum(x, na.rm=T))
    } else {
      sums <- tapply(v$Qs, list(v$sub, v$year), function(x) sum(x, na.rm=T))
    }
    ranks <- apply(sums, 1, function(x) {
      # NAs <- length(which(is.na(x)))
      out <- order(x, decreasing=FALSE, na.last=NA)
    })
  }
  
  rnk0 <- rankyears(md1, type='o')
  rnk1 <- rankyears(md1)
  rnk2 <- rankyears(md2)
  rnk3 <- rankyears(md3)
  # rnk4 <- rankyears(v4)
  
  whichyear <- function(sim, obs) {
    lowestobs <- obs[1]
    return(which(sim == lowestobs))
  }
  
  matchyears <- function(rnk) {
    return(mapply(whichyear, rnk, rnk0))
  }
  
  m1 <- matchyears(rnk1)
  m2 <- matchyears(rnk2)
  m3 <- matchyears(rnk3)
  # m4 <- matchyears(rnk4)
  
  par(mfrow = c(1, 1))
  
  hist(m3, col=rgb(0,0,1,0.4), breaks=10, border=NA, main = 'Simulated ranks of driest year by location')
  hist(m1, col=rgb(1,0,0,0.4), breaks=10, border=NA, add=T)
  legend(x=5, y=40, legend=c("Naive model", "Final model"), 
         fill=c(rgb(1,0,0,0.4), rgb(0,0,1,0.4)),
         border=NA,
         bty="n")
  
}

### Alphas =====================================================================
{ # alphas
  adf <- data.frame(region = unlist(lapply(1:6, function(x) rep(x, 12))),
                    month = rep(1:12, 6),
                    x = unlist(alphas_x), 
                    q = unlist(alphas_q))   # data frame of alphas
  idf <- adf[, 1:2]
  idf$I <- unlist(lapply(1:nrow(idf), function(x) {
    cat(idf$region[x], idf$month[x], "\n")
    p90 <- pcp90[pcp90$region == idf$region[x] & pcp90$mon == idf$month[x], ]
    return(sd(p90$obs))
  }))
  idf <- tapply(idf$I, list(idf$region, idf$mon), function(x) x)
  
  pal <- brewer.pal(6, 'Accent')
  par(mfrow=c(3,1))
  for(i in 1:6) {
    if(i == 1) {
      plot(1:12, adf$x[adf$region == i], type='l', col=pal[i], 
           lwd = 2, ylim = c(0, 16))
    } else {
      lines(1:12, adf$x[adf$region == i], col=pal[i], lwd=2) 
    }
  }
  for(i in 1:6) {
    if(i == 1) {
      plot(1:12, adf$q[adf$region == i], type='l', col=pal[i], 
           lwd = 2, ylim = c(0, 11))
    } else {
      lines(1:12, adf$q[adf$region == i], col=pal[i], lwd=2) 
    }
  }
  
  for(i in 1:6) {
    if(i == 1) {
      plot(idf[i,], type='l', col=pal[i], lwd=2, 
           ylim=c(0, max(idf, na.rm=T) + 0.05))
    } else {
      lines(idf[i,], col=pal[i], lwd=2)
    }
  }
}




### Appendix ===================================================================
# 
# 
# { # more outliers
#   par(mfrow=c(1,1))
#   plot(wshed, col='gray', border=NA)
#   plot(wshed[wshed$Subbasin %in% c(370, 375, 397), ], col='red', add=T)
# }
# {
#   ## 90s?
#   r9 <- read.csv("results/8500_STRI_LU_HM_GWR_rch.csv")
#   r9 <- r9[, c("SUB", "YEAR", "MON", "FLOW_OUTcms")]
#   names(r9)[4] <- "FLOW_OUT"
#   
#   r9 <- r9[order(r9$SUB, r9$YEAR, r9$MON),]
#   
#   lookup <- readRDS("data/commondata/980to915nbr.RDS")
#   r9$SUB <- sapply(r9$SUB, function(x) which(lookup == x))
#   
#   r00 <- r0[r0$year %in% r9$YEAR,]
#   r00 <- tapply(r00$q, list(r00$sub, r00$mon), function(x) mean(x, na.rm=T))
#   r99 <- r9[r9$SUB %in% r0$sub, ]
#   r99 <- tapply(r9$FLOW_OUT, list(r9$SUB, r9$YEAR), mean)
#   
#   pstat <- read.csv("data/commondata/available_pstations_00s.csv")
#   pstat2 <- SpatialPoints(pstat[, c(2,3)], proj4string=wgs84)
#   pstat2 <- spTransform(pstat2, crs(wshed))
#   p2 <- over(pstat2, wshed)
#   p2 <- p2[, c("Subbasin", "region")]
#   pstat <- cbind(pstat, p2)
#   pstat <- pstat[-which(is.na(pstat$Subbasin)), ]
#   saveRDS(pstat, "data/commondata/active_pstations_00.RDS")
# }


# n1 <- sum(tapply(md1$res, md1$sub, function(x) !all(is.na(x))))
# n2 <- sum(tapply(md2$res, md2$sub, function(x) !all(is.na(x))))
# 
# # length(names(upst))
# # r0 <- r0[r0$year %in% 2006:2015, ]
# # upst <- upst[sort(unique(r0$sub))]
# # ind <- unlist(lapply(upst, function(u) {
# #   try <- unlist(lapply(upst, function(v) is_subset(u, v)))
# #   if(sum(try) == 1) return(1)
# # }))
# # upst2 <- upst[names(ind)]
# # plot(wshed, col='gray', border=NA)
# # plot(wshed[wshed$Subbasin %in% aa,], col='blue', border=NA, add=T) # total basins
# # plot(wshed[wshed$Subbasin %in% bb,], col='red', border=NA, add=T) # terminuses
# # plot(wshed[wshed$Subbasin %in% cc,], col='yellow', add=T)
# # 
# # rch3 <- r3[r3$SUB %in% as.numeric(names(upst2)), ]
# # md <- match_day(rch3)
# # md2 <- match_day(r3)
# # 
# # vv <- validate(md, type='fl', nfl=3)
# # 
# # vv2 <- validate(match_day(r3), type='fl', nfl=3)

# n <- tapply(v$Qo, tapplygrp, length)
# 
# # col1 <- switch(cols, 
# #                'none' = rep('black', length(qobs)),
# #                'mon' = unlist(lapply(1:12, function(n) {
# #                  rep(brewer.pal(12, 'Paired')[n], nrow(qobs))
# #                })))


# leaps <- is_leap(simyears[1:nskip])
# rch <- rch[((nskip * nsub * ifelse(leaps, 366, 365)) + 1):nrow(rch), ]
# rch$YEAR <- unlist(lapply(simyears[(nskip + 1):length(simyears)], 
#                           function(y) {
#                             nday <- ifelse(is_leap(y), 366, 365)
#                             rep(y, nsub * nday)
#                           }))[1:nrow(rch)]


# 
# if(any(is.infinite(qsim))) {
#   qsim[is.infinite(qsim)] <- NA  
# }

# par(mfrow=c(1,1))
# if(group=='sub') {
#   qobs <- rowMeans(qobs); qsim <- rowMeans(qsim)
# }

# add_months <- function(rch, simyears=2005:2015, nskip=1) {
#   nsub <- length(unique(rch$SUB))
#   leaps <- is_leap(simyears[1:nskip])
#   rch <- rch[((nskip * nsub * ifelse(leaps, 366, 365)) + 1):nrow(rch),]
# }

# 
# switch(cols,
#        'none' = print("no legend"),
#        'mon' = legend("bottomright", fill=brewer.pal(12, 'Paired'), bty="n",
#                       legend=month_names),
#        'sub' = legend())