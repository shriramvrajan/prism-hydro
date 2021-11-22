rm(list = ls())
library('rWind')
library('plotrix')

## Wind data from https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdQAwindmday.html
## The wind data will be in vectors (u and v for x and y).  
## I found the R package 'rWind' to be useful in transforming that data to speed and direction.
## gcWind (Gulf of Chiriqui), gpWind(Gulf of Panama), apWind (Azuero Peninsula). 
## Data subset corresponding to 1x1 degree grids ~30-km south from shore.
pcp90 <- readRDS("data/1990s_precipitation.RDS")
adf <- readRDS("data/alpha_values.RDS")

wgc <- read.csv("../data/gcWind.csv"); wgc$loc <- "gc"
wgp <- read.csv("../data/gpWind.csv"); wgp$loc <- "gp"
wap <- read.csv("../data/apWind.csv"); wap$loc <- "ap"

wnd <- rbind(wgc, wgp, wap)
wnd <- wnd[-which(is.na(wnd$u) | is.na(wnd$v)), ]
wnd <- cbind(wnd, uv2ds(wnd$u, wnd$v)) ## adding direction and speed
wnd$mon <- sapply(wnd$date, function(s) {
  return(as.numeric(substr(s, 6, 7)))
})

wdir <- tapply(wnd$dir, list(wnd$loc, wnd$mon), mean)
wspd <- tapply(wnd$speed, list(wnd$loc, wnd$mon), mean)


## Alphas/
pal <- brewer.pal(6, 'Accent')
par(mfrow=c(3,1))
for(i in 1:6) {
  if(i == 1) {
    plot(1:12, adf$x[adf$region == i], type='l', col=pal[i], 
         lwd = 2, ylim = c(0, 16), xlab = "Month", 
         # main = "Strength of interpolation parameters\n by region and month",
         ylab = expression(paste(alpha, ", decay parameter for occurrence")))
  } else {
    lines(1:12, adf$x[adf$region == i], col=pal[i], lwd=2) 
  }
}
abline(v = c(4, 10))

plot(wdir[2,], type='l', col='red', ylab='wind direction')
lines(wdir[1,], col='blue'); lines(wdir[3,], col='green'); abline(v=c(4,10))

plot(wspd[3,], type='l', col='green', ylab='wind speed', ylim=c(0,7.5))
lines(wspd[2,], col='red'); lines(wspd[1,], col='blue'); abline(v=c(4,10))

for(i in 1:6) {
  if(i == 1) {
    plot(1:12, adf$q[adf$region == i], type='l', col=pal[i], 
         lwd = 2, ylim = c(0, 11), xlab = "Month",
         ylab = expression(paste(alpha, ", decay parameter for quantity")))
  } else {
    lines(1:12, adf$q[adf$region == i], col=pal[i], lwd=2) 
  }
}
abline(v = c(4, 10))



