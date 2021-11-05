## Generally useful things to have for the type of stuff I do

## Libraries ===================================================================
library(raster)
library(hydroGOF)
library(spgwr)
library(RColorBrewer)
library(foreach)
library(doParallel)
library(ggplot2)
library(fitdistrplus)
library(mgcv)
library(lhs)
library(scales)
library(ape)
library(moments)

## CRS definitions =============================================================
wgs84 <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
utm17 <- CRS(paste0("+proj=utm +zone=17 +datum=WGS84 +units=m", 
                    " +no_defs +ellps=WGS84 +towgs84=0,0,0"))

## Functions ===================================================================

# Moving average
m_a <- function(x, n = 5, sides = 1) { 
  filter(x, rep(1/n, n), sides = sides)
}

# Moving sum
m_s <- function(x, n = 5, sides = 1) { 
  filter(x, rep(1, n), sides = sides)
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

inv_dist_weight <- function(x, y) {
  
}

'%nin%' <- Negate('%in%')

extractNA <- function(x, y) {
  # set.seed(2)
  # 
  # # create a 10x10 raster
  # r <- raster(ncol=10,nrow=10, xmn=0, xmx=10, ymn=0,ymx=10)
  # r[] <- 1:10
  # r[sample(1:ncell(r), size = 25)] <- NA
  # 
  # # plot the raster
  # plot(r, axes=F, box=F)
  # segments(x0 = 0, y0 = 0:10, x1 = 10, y1 = 0:10, lty=2)
  # segments(y0 = 0, x0 = 0:10, y1 = 10, x1 = 0:10, lty=2)
  # 
  # # create sample points and add them to the plot
  # xy = data.frame(x=runif(10,1,10), y=runif(10,1,10))
  # points(xy, pch=3)
  # text(x = xy$x, y = xy$y, labels = as.character(1:nrow(xy)), pos=4, cex=0.7, xpd=NA)
  # 
  # # use normal extract function to show that NAs are extracted for some points
  # extracted = extract(x = r, y = xy)
  
  # then take the raster value with lowest distance to point AND non-NA value in the raster
  sampled = apply(X = y, MARGIN = 1, 
                  FUN = function(y) {
                    r@data@values[which.min(replace(distanceFromPoints(r, xy), 
                                                    is.na(r), NA))]})
  
  # show output of both procedures
  print(data.frame(xy, extracted, sampled))
}


## Time ========================================================================

month_indices <- list(1:31, 32:59, 60:90, 91:120, 121:151, 152:181, 182:212,
                213:243, 244:273, 274:304, 305:334, 335:365)
month_indices_leap <- month_indices
month_indices_leap[[2]] <- 32:60
for(m in 3:12) {
  month_indices_leap[[m]] <- month_indices_leap[[m]] + 1
}
rm(m)

month_days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
month_days_leap <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
month_days2 <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

day_months <- unlist(lapply(1:12, function(x) { 
  rep(x, month_days[x]) }))
day_months_leap <- unlist(lapply(1:12, function(x) { 
  rep(x, month_days_leap[x]) }))

is_leap <- function(y) {
  unlist(lapply(y, function(y) {
  if(y %% 4 == 0 & (y %% 100 != 0 | y %% 400 == 0)) {
    return(TRUE)
  } else {
    return(FALSE)
  }}))
}

month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

p_factorize <- function(n) {
  
  # series <- 1:(n/2)
  # factors <- vector()
  # 
  # for(i in series) {
  #   n <- n/i
  #   if(n > 2)
  # }
  # 
  # return(factors)
}

## Data ========================================================================
# 
# pan_shp <- shapefile("~/panama_scripts/commondata/pan_shp/PAN_adm0.shp")
# 
# nbr <- readRDS("~/panama_scripts/commondata/980to915nbr.RDS")
