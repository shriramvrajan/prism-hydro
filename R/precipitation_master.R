# ideally would like this script to do the following
# 1a. fitting parameters for interpolation 
# 1b. forecasting based on means / other switches
# 2. interpolate based on forecasted seeds and/or fitted parameters
# 3. change file.cio - precipitation record name
# 4. change temperature file
# 5. run SWAT

### SWITCHES ===================================================================

source("R/precipitation_functions.R")

datecurr <- as.character(Sys.Date())

actual_interp <- TRUE ## Only F for forecasting

int_parallel = T                     ## use multiple CPU cores?
ncores <- ifelse(int_parallel, 3, 1) ## how many CPU cores  
int_outname = paste0("pcp_out_", datecurr, ".RDS")
int_pcpname = "pcp2020rd.pcp" 
int_scale = 'n'                 ## scaling to fix 2-step quantity loss

val_fnam <- paste0("0514_daily_", datecurr, ".rch")

do_pfitting = T
do_interp = F
do_scaling = F
do_swat = F
do_val = F

### 1a. Fitting =================================================================

if(do_pfitting) {
  Log("Fitting rainfall parameters...")
  opt2 <- all_interp(step="fit", algorithm="1step", pcpraw = pcp90) 
  p_3 <- all_interp(pcpraw=pcp00, parset=opt2, algorithm="1step", 
                    step="gen", dis=dm1)
  ## troubleshoot HERE
}


### 1b. Forecasting =============================================================

# if(do_forecast) {
#   setwd("~/panama_scripts/weather/data/forecasting")
#   Log("Forecasting seeds..")
#   out1 <- generate_fullpcp(startdate='2005-01-01', enddate='2015-01-31', XSAC=T,
#                            means=fc_means, Qmodel='gamma')
#   saveRDS(out1, fc_outname)
#   plot_pcp(out1)
# }

### 2. Interpolating ===========================================================

if(do_interp) {
  source("R/precipitation_functions.R")
  
  Log("Interpolating to all locations...")
  if(int_parallel == TRUE) {
    registerDoParallel(cores = ncores)
    Log(paste0("number of workers: ", getDoParWorkers()))
  }
  if(exists("out1")) {
    pcp01 <- out1 
  } else if(fc_outname %in% list.files() & !actual_interp) {
    pcp01 <- readRDS(fc_outname)
  } else if(actual_interp) {
    pcp01 <- pcp00
  }
  # names(pcp01)[7] <- "obs"
  pcp01$obs[which(is.na(pcp01$obs))] <- 0
  pcp01$wet_obs <- as.numeric(pcp01$obs >= 0.5)
  
  p_out <- interpolate_full(pcp=pcp01, pred1="obs", pname="NA", alg="1step",
                            tofile=F)
  saveRDS(p_out, paste0("output/", int_outname))
}

### 3. Scaling and writing .pcp file ===========================================
if(do_scaling) {
  Log("Processing and scaling interpolated precipitation...")
  p_out <- readRDS(paste0("output/", int_outname))
  p_sc <- postprocess_and_scale(p_out, scale=int_scale)

  writeLines(p_sc, paste0("output/", int_pcpname))
} else {
  ## no need to scale for 1-step algorithm
  writeLines(p_out, paste0("output/", int_pcpname))
}

### 4. Changing temperature based on WC ========================================
if(do_temp) {
  Log("Changing temperature values...")
  write_temp(type=wc_type, rcp=wc_rcp, model=wc_model, year=wc_year,
             scen="~/panama_swat3/Scenarios/sim1_p2step/TxtInOut/")
}

### 5. Editing file.cio and running SWAT =======================================
fpath <- paste0("../../rchfiles/", val_fnam)
if(do_swat) {
  Log("Setting precipitation input for SWAT run...")
  edit_pcpin(pname = int_pcpname)
  
  run_swat("~/panama_swat3/Scenarios/sim1_p2step/TxtInOut/")
  
  system2("mv", args=c("output.rch", fpath))
  
  # file.copy("~/panama_swat3/Scenarios/sim1_p2step/TxtInOut/output.rch",
  #            "~/panama_swat3/Scenarios/rchfiles/") ##
  # setwd("~/panama_swat3/Scenarios/rchfiles/")
  # file.rename("output.rch", val_fnam)
  ## move it because read_rch messes up in TxtInOut for some reason
  Log("SWAT run finished.")
}

### 6. Validating SWAT output ==================================================
if(do_val) {
  Log("Validating...")
  
  r1 <- read_rch(fpath)
  v1 <- validate_output(r1, simyears=2005:2009, nskip=1, daily=T)
  p1 <- plot1(v1, out='df')
  Log(paste0("NSE:", NSE(p1$qsim, p1$qobs)))
  Log(paste0("R^2:", rsq(p1$qsim, p1$qobs)))
  
  saveRDS(r1, paste0("~/panama_scripts/discharge/data/r_", "wcmay5", ".RDS"))
}

### APPENDIX ===================================================================

### 1. Notes on speed ==========================================================

# Precipitation generation : 07 min
# Wet/dry interpolation    : 07 min
# Quantity interpolation   : 11 min
# Generating .pcp file     : 04 min
# SWAT run (10 years)      : 50 min
#----------------------------------
# Total                    : 79 min

# STill need to figure out SD
# why doesnt z score work?
# if(FALSE) {
#   coords <- merge(sample(1:100, 30), sample(1:100, 30)); n <- nrow(coords)
#   # mu_mon <- 1:nrow(coords)
#   mu_mon <- runif(n, min=2, max=10)
#   k <- 0.6741434 # empirically derived, see gamma section in appendix
#   theta <- mu_mon / k
#   sd_mon <- mu_mon/(sqrt(k))
#   # mean(rgamma(10000, shape=k, scale=theta[100]))
#   
#   type='gamma'
#   k <- rep(k, n)
#   
#   W <- generate_weights(coords, 1)
#   
#   gamma <- 0.8
#   
#   VV <- lapply(1:2000, function(i) {
#     cat(i,"\n")
#     # test1 <- lapply(1:5000, function(x) {
#     #   u <- switch(type,
#     #               'skewnormal'=Rd(type=type, mu=mu_mon, sd=sd_mon, sk=skew_mon),
#     #               'exp'=Rd(type=type, rate=rates),
#     #               'gamma'=Rd(type=type, k=k, theta=theta))
#     #   # u <- Rd(mu_mon, sd_mon, skew_mon)
#     #   # return(gamma * W %*% as.vector(u) + u)
#     #   out <- gamma * W %*% as.vector(u) + (1 - gamma) * u
#     #   return(out)
#     # })
#     # test1 <- do.call(cbind, test1)
#     # # mu_gen <- rowMeans(test1)
#     # # sd_gen <- unlist(lapply(1:nrow(test1), function(x) {
#     # #   sd(test1[x,])
#     # # }))pww*pw + pwd*(1-pw)
#     # mu_gen <- rowMeans(test1)
#     # sd_gen <- unlist(lapply(1:n, function(r) {
#     #   sd(test1[r,])
#     # }))
#     mu_gen <- gamma * W %*% mu_mon + (1-gamma) * mu_mon
#     
#     c <- mu_mon/mu_gen
#     
#     u0 <- Rd(type=type, k=k, theta=theta)
#     V <- gamma * W %*% as.vector(u0) + (1-gamma) * u0
#     V1 <- V*c
#     ## Rescale by z-score + mean + sigma instead here
#     # z_v <- (V - mu_gen)/sd_gen
#     # V1 <- mu_mon + sd_mon * z_v
#     return(V1)
#   })
#   VB <- do.call(cbind, VV)
#   
#   gammafits <- lapply(1:nrow(VB), function(x) {
#     cat(x, "\n")
#     fitdist(VB[x,], distr='gamma')
#   })
#   gammafits2 <- lapply(1:nrow(VB), function(x) {
#     cat(x, "\n")
#     fitdist(VB[x,]/c[x], distr='gamma')
#   })
#   k2 <- unlist(lapply(gammafits2, function(x) x$estimate[1]))
#   theta2 <- unlist(lapply(gammafits2, function(x) 1/x$estimate[2]))
#   # plot(1:nrow(coords), rowMeans(VB)); abline(0,1)
#   mu_gen1 <- as.vector(gamma * W %*% mu_mon + (1-gamma) * mu_mon)
#   mu_gen2 <- rowMeans(VB)
#   c <- mu_mon/mu_gen1
#   par(mfrow=c(1,2))
#   plot(mu_mon, mu_gen2);abline(0,1)
#   plot(sd_mon, sd_gen);abline(0,1)
#   
#   # sd_gen <- unlist(lapply(1:nrow(VB), function(x) sd(VB[x,]/c, na.rm=T)))
#   # sk_gen <- unlist(lapply(1:nrow(VB), function(x) skewness(VB[x,]/c, na.rm=T)))
#   sd_gen <- unlist(lapply(1:nrow(VB), function(x) sd(VB[x,]/c[x], na.rm=T)))
#   # sd_sim <- gamma * W %*% sd_mon + (1 - gamma) * sd_mon # can't just do this for sd
#   
#   
#   sd_test1 <- unlist(lapply(1:nrow(W), function(x) {
#     sum(W[x, ] * sd_mon)
#   }))
#   
#   eig <- eigen(W)
#   
#   ### checking some stuff about the pcp files
#   
#   setwd("~/panama_swat3/Scenarios/sim1_p2step/")
#   pcplist <- list.files()[grep(".pcp", list.files())]
#   
#   
#   pcp1 <- readLines("bst00.pcp")
#   pcp2 <- readLines("emp8.pcp")
#   pcp3 <- readLines("pexp.pcp")
#   pcp3 <- readLines("p_wc00.pcp")
#   pcp4 <- readLines("p_wc00sc.pcp")
#   
#   pcp_to_df <- function(pcp, sums=T) {
#     df1 <- lapply(5:length(pcp), function(t) {
#       cat(t, "\n")
#       linet <- pcp[[t]]
#       year <- as.numeric(substr(linet, 1, 4))
#       day <- as.numeric(substr(linet, 5, 7))
#       val <- as.numeric(sapply(seq(from=8, to=nchar(linet), by=5), 
#                                function(i) substr(linet, i, i+4)))
#       return(c(year, day, val))
#     })
#     df1 <- as.data.frame(do.call(rbind, df1))
#     names(df1) <- c("y", "d", paste0("s", 1:980))
#     if(sums) return(colSums(df1)[3:982])
#     return(df1)
#   }
#   
#   sum1 <- pcp_to_df(pcp1, sums=F)
#   sum2 <- pcp_to_df(pcp2, sums=F)
#   sum3 <- pcp_to_df(pcp3, sums=F)
#   
#   par(mfrow=c(3,4))
#   teststations <- sample(3:982, 12)
#   for(i in teststations) {
#     plot(m_a(sum1[, i], 60), type='l')
#     lines(m_a(sum2[, i], 60), col='red')
#     lines(m_a(sum3[, i], 60), col='blue')
#   }
#   
#   plot(sum1, sum2);abline(0,1)
# 
#   
#   ### SKEWNESS
#   pcpwet <- pcp00[pcp00$wet_obs == 1,]
#   skew <- tapply(pcpwet$obs, list(pcpwet$id, pcpwet$month), skewness)
#   kurt <- tapply(pcpwet$obs, list(pcpwet$id, pcpwet$month), kurtosis)
# 
#   hist(skew); abline(v=2/sqrt(k[1]), col='red')
#   hist(kurt); abline(v=6/k[1], col='red')    
#   
#   ## also try turning things off
#   
#   
#   ### simulation
#   
#   xy <- data.frame(x = runif(30, 0, 10), y = runif(30, 0, 10))
#   xymeans <- runif(30, 0, 40)
#   xyk <- runif(30, 0, 10)
#   xyk2 <- rep(0.6, 30) 
#   xytheta <- xymeans/xyk
#   xytheta2 <- xymeans/xyk2
#   
#   vals <- data.frame(t = 1:10000)
#   for(n in 1:30) {
#     tsk1 <- rgamma(nrow(vals), shape=xyk[n], scale=xytheta[n])
#     tsk2 <- rgamma(nrow(vals), shape=xyk2[n], scale=xytheta2[n])
#     vals <- cbind(vals, tsk1, tsk2)
#   }
#   vals <- vals[, c(seq(2, ncol(vals), 2), seq(3, ncol(vals), 2))]
#   names(vals)[1:ncol(vals)] <- c(paste0("tk1", 1:30), paste0("tk2", 1:30))
#   row.names(vals) <- 1:nrow(vals)
# # names(vals)[2:62] <- paste0("tsk", 1:30)
#   # vals <- vals[,2:31]; row.names(vals) <- 1:1000
#   par(mfrow=c(2,2))
#   plot(colMeans(vals)[1:30], xymeans); abline(0,1)
#   points(colMeans(vals)[31:60], xymeans, col='red')
#   plot(lapply(vals[,1:30], sd), sqrt(xyk)*xytheta); abline(0,1)
#   points(lapply(vals[,31:60], sd), sqrt(xyk)*xytheta, col='red')
#   plot(lapply(vals[,1:30], skewness), 2/sqrt(xyk)); abline(0,1)
#   points(lapply(vals[,31:60], skewness), 2/sqrt(xyk), col='red')
#   plot(lapply(vals[,1:30], kurtosis), (6/xyk)+3); abline(0,1)
#   points(lapply(vals[,31:60], kurtosis), (6/xyk)+3, col='red')
#   
#   p <- data.frame(x=c(2,7,4), y=c(2,4,6))
#   plot(xy$x, xy$y, pch=17, cex=0.8); points(p$x, p$y, col='red', pch=19)
#   d <- as.matrix(dist(rbind(xy, p)))[31:33,1:30]
#   
#   intp <- matrix(nrow=nrow(vals), ncol=nrow(p)*2)
#   
#   for(j in 1:nrow(p)) {
#     for(i in 1:nrow(vals)) {
#       cat(j, i, "\n")
#       neighbors1 <- vals[i,1:30]
#       neighbors2 <- vals[i,31:60]
#       weights <- exp(-d[j,])
#       intp[i,j] <- sum(neighbors1*weights) / sum(weights)
#       intp[i,j+3] <- sum(neighbors2*weights) / sum(weights)
#       # intp[i,j] <- 
#     }
#   }
#   
#   intp0 <- intp
#   
#   
# }
# plot(sd_mon, sd_gen)
# abline(0,1)

# a <- sd_mon/sd_gen
# b <- mu_mon/mu_gen
# # plot(a,b)
# par(mfrow=c(1,1))
# plot(a, col='blue', pch=19, cex=0.8, ylim=range(c(a,b)))
# points(b, col='red', pch=19, cex=0.8)

## temperature scratch 

# setwd("~/panama_scripts/weather/data")
# wgen <- read.csv("wgen_hm_gwr_centroids.csv")





# p11$wet_obs <- as.numeric(runif(nrow(p11)) < p11$w_interp3)
# p11$pcp <- p11$predidw
# p11$pcp[p11$wet_obs == 0] <- 0
# p11$pcp[p11$wet_obs == 1 & p11$pcp < 0.5] <- 0.5
# p11$pcp[is.na(p11$pcp)] <- 0
# 
# scalars <- readRDS("scalars_combinedmodel.RDS")
# p_scaled <- p11
# p_scaled$pcp <- unlist(lapply(1:nrow(p11), function(x) {
#   if(x %% 1000 == 0) cat(x, "\n")
#   return(p11$pcp[x] * scalars[p11$mon[x], p11$region[x]])
# }))
# p_sc <- write_swatpcp(dat=p_scaled, output_type=".pcp")
# p_final <- write_swatpcp(dat=p11, output_type=".pcp")
# setwd("~/panama_swat3/Scenarios/sim1_p2step/TxtInOut/")
# writeLines(p_final, int_pcpname)
# 
# 
# source("~/panama_scripts/weather/rscripts/interpolation/fitting_and_imputing.R")
# 

# rm(list = ls())
# setwd("~/panama_scripts/weather/data/hm_90s/interp2")
# 
# p0 <- readRDS("pcp_out_FINALMODELFORSWAT00.RDS")
# p1 <- readRDS("pcp_out_expmodel_00smeans.RDS")
# p2 <- readRDS("pcp_out_interp_00s_wclimgamma.RDS")
# p3 <- readRDS("pcp_out_interp_00sgamma.RDS")
# # 
# p0$pcp <- p0$predidw * p0$w_interp2
# # 
# m0 <- tapply(p0$pcp, list(p0$sub, p0$mon), function(x) mean(x, na.rm=T))
# m1 <- tapply(p1$pcp, list(p1$sub, p1$mon), function(x) mean(x, na.rm=T))
# m2 <- tapply(p2$pcp, list(p2$sub, p2$mon), function(x) mean(x, na.rm=T))
# m3 <- tapply(p3$pcp, list(p3$sub, p3$mon), function(x) mean(x, na.rm=T))
# # 
# par(mfrow=c(1,3))
# plot(m0, m1, main='exp'); abline(0,1)
# plot(m0, m2, main='gamma fc'); abline(0, 1)
# plot(m0, m3, main='gamma int'); abline(0, 1)
# # 
# sd0 <- tapply(p0$pcp, list(p0$sub, p0$mon), function(x) sd(x, na.rm=T))
# sd1 <- tapply(p1$pcp, list(p1$sub, p1$mon), function(x) sd(x, na.rm=T))
# sd2 <- tapply(p2$pcp, list(p2$sub, p2$mon), function(x) sd(x, na.rm=T))
# sd3 <- tapply(p3$pcp, list(p3$sub, p3$mon), function(x) sd(x, na.rm=T))
# # 
# plot(sd0, sd1, main='exp var'); abline(0,1)
# plot(sd0, sd2, main='gamma var');abline(0,1)
# plot(sd0, sd3, main='gamma int var');abline(0,1)

# 
# sk0
# 
# 
# s0 <- readRDS("../../hm_00s/gwr_daily/gwr_lin/pmodeldata_lin_full.rds")
# s1 <- readRDS("../../forecasting/TEST_OUT_EXP_XMODELCOUNTS_00.RDS")
# s2 <- readRDS("../../forecasting/TEST_OUT_WCLIMPRESENTGAMMA.RDS")
# 
# m00 <- tapply(s0$obs, list(s0$id, s0$mon), function(x) mean(x, na.rm=T))
# m11 <- tapply(s1$pcpmm, list(s1$id, s1$month), function(x) mean(x, na.rm=T))
# m22 <- tapply(s2$pcpmm, list(s2$id, s2$month), function(x) mean(x, na.rm=T))
# 
# ids <- intersect(intersect(rownames(m00), rownames(m11)), rownames(m22))
# m00 <- m00[ids,]; m11 <- m11[ids,]; m22 <- m22[ids,]












##### SWAT thing
# library(raster)
# setwd("~/panama_swat4/Watershed/")
# wshd <- shapefile("watershed.shp")
# rch <- shapefile("reach.shp")
# 
# wshd <- wshd[, c("OBJECTID", "Area", "Subbasin")]
# names(wshd)[1] <- "PolygonID"
# writeOGR(wshd, ".", "watershed2", driver="ESRI Shapefile")
# 
# setwd("~/panama_swat3/Watershed/Shapes")
# a <- shapefile("riv1.shp")
