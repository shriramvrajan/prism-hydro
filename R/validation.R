library(raster)
library(hydroGOF)
library(reshape2)
library(RColorBrewer)

### Data =======================================================================

r0 <- readRDS("data/results/obs_hydro.RDS")
r1 <- readRDS("data/results/SWATnaive_model.RDS")
r2 <- readRDS("data/results/SWATsampling_model.RDS")
r3 <- readRDS("data/results/SWATfinal_model.RDS")

# outliers <- c(999999,999399)
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
# alpha values for wet/dry interpolation
alphas_q <- readRDS("data/opt_qidw_noelev.RDS") 
# alpha values for quantity

pcp90 <- readRDS("~/panama_scripts/weather/data/hm_90s/interp2/pmodel_winterp.RDS") 

## Functions ===================================================================

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

validate <- function(v, nfl = 3, by = 'subm', type='mean', ...) { 
  
  ## nfl is number of maxima for flood analysis
  ## valid 'by'   : 'sub', 'subm', 'suby'
  if(any(v$sub %in% outliers)) {
    v <- v[-which(v$sub %in% outliers),]
  }
  
  v <- v[which(v$sub %in% names(upst)), ]
  
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
  plot(q2$qo, q2$qs, pch=19, ...)
  return(q2)
  
}

### Mean SD max ================================================================

md1 <- match_day(r1)
md2 <- match_day(r2)
md3 <- match_day(r3)

{
  png(filename="plots/fig3.png", width=800, height=800, 
      pointsize=20, units="px")
  par(mfrow = c(3,3))
  cex1 = 0.5
  
  # mean
  v1 <- validate(md1, cex=cex1, main='Mean discharge,\n NN model')
  v2 <- validate(md2, cex=cex1, main='Mean discharge,\n IWGEN model')
  v3 <- validate(md3, cex=cex1, main='Mean discharge,\n RDW model')
  
  # sd
  v1s <- validate(md1, cex=cex1, type='sd', main='SD discharge,\n NN model')
  v2s <- validate(md2, cex=cex1, type='sd', main='SD discharge,\n IWGEN model')
  v3s <- validate(md3, cex=cex1, type='sd', main='SD discharge,\n RDW model')
  
  # max
  nfl1 = 3 # top n max days
  v1e <- validate(md1, cex=cex1, type='fl', nfl=nfl1,
                  main='Max. discharge,\n NN model')
  v2e <- validate(md2, cex=cex1, type='fl', nfl=nfl1,
                  main='Max. discharge,\n IWGEN model')
  v3e <- validate(md3, cex=cex1, type='fl', nfl=nfl1,
                  main='Max. discharge,\n RDW model')
  dev.off()
}


## Alphas ======================================================================

{
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
  
  png(filename="plots/fig5.png", width=450, height=800,
      pointsize=15, units="px")
  par(mfrow=c(2,1))
  for(i in 1:6) {
    if(i == 1) {
      plot(1:12, adf$x[adf$region == i], type='l', col=pal[i], 
           lwd = 2, ylim = c(0, 16), xlab = "Month", 
           main = "Strength of interpolation parameters\n by region and month",
           ylab = expression(paste(alpha, " for interpolating occurrence")))
    } else {
      lines(1:12, adf$x[adf$region == i], col=pal[i], lwd=2) 
    }
  }
  abline(v = c(4, 10))
  for(i in 1:6) {
    if(i == 1) {
      plot(1:12, adf$q[adf$region == i], type='l', col=pal[i], 
           lwd = 2, ylim = c(0, 11), xlab = "Month",
           ylab = expression(paste(alpha, " for interpolating quantity")))
    } else {
      lines(1:12, adf$q[adf$region == i], col=pal[i], lwd=2) 
    }
  }
  abline(v = c(4, 10))
  
  dev.off()
}

### Within-basin analysis ======================================================

b_nse <- function(v) {
  return(by(v, v$sub, function(x) {
    return(NSE(x$qs, x$qo))
  }))
}

bn1 <- b_nse(v1)
bn2 <- b_nse(v2)
bn3 <- b_nse(v3)

### Predictors of failure ======================================================

# GLM: NSE ~ upstream * elev * nsub (3 outliers)
basin_nse <- by(v3, v3$sub, function(v) {
  return(NSE(v$qs, v$qo))
})
basin_rsq <- by(v3, v3$sub, function(v) {
  return(rsq(v$qs, v$qo))
})

# getting distances to nearest station & n stations
met_st <- readRDS("data/met_stations.RDS")
hyd_st <- readRDS("data/hyd_stations.RDS")

pof <- data.frame(sub = names(basin_nse),
                  nse = as.vector(basin_nse),
                  rsq = as.vector(basin_rsq))

# adding elevation
pof$elev <- sapply(pof$sub, function(sub) {
  return(hyd_st$elev[which(hyd_st$sub == sub)])
})

# adding dist. to nearest p station
lonlat <- rbind(cbind(met_st$lon, met_st$lat), 
                cbind(hyd_st$lon, hyd_st$lat))

dm1 <- as.matrix(dist(lonlat, diag=T, upper=T))[204:262, 1:203]
# columns are meteorological stations, rows are hydrological stations
nearest <- sapply(1:nrow(dm1), function(h) {
  min(dm1[h,])
})
names(nearest) <- hyd_st$sub
pof$near <- nearest[pof$sub]

# adding number of pcp stations by region
pof$region <- sapply(pof$sub, function(s) {
  return(wshed$region[which(wshed$Subbasin == s)])
})
ngauge <- tapply(met_st$region, met_st$region, length)
pof$ngauge <- ngauge[pof$region]

# adding number of upstream/downstream subs 
pof$up <- sapply(pof$sub, function(s) {
  length(upst[[as.character(s)]])
})
pof$down <- sapply(pof$sub, function(s) {
  length(downst[[as.character(s)]])
})
pof$nsub = pof$up + pof$down # watershed size

subs_per_region <- tapply(wshed$region, wshed$region, length)
pof$propg <- pof$ngauge/subs_per_region[pof$region]

# nse ~ down * elev * ngauge
lm1 <- lm(nse ~ down * elev * propg, data=pof)
plot(pof$nse, predict.lm(lm1), xlab = "obs", ylab = "sim")
rsq(pof$nse, predict.lm(lm1))
summary(lm1)

