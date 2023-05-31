rm(list = ls())
source("R/validation_functions.R")

### Data =======================================================================
r1 <- readRDS("data/results/SWATnaive_model.RDS")
r2 <- readRDS("data/results/SWAT1step_model.RDS")
r3 <- readRDS("data/results/SWAT2step_model.RDS")

simyears1 <- 2005:2015
redo_figures <- FALSE

### Paper figures :
### Mean SD max ================================================================

md1 <- match_day(r1, simyears = simyears1)
md2 <- match_day(r2, simyears = simyears1)
md3 <- match_day(r3, simyears = simyears1)

if (redo_figures) {
  png(
    filename = "figures/fig3.png", width = 700, height = 900,
    pointsize = 20, units = "px"
  )
  cex1 <- 0.5
  nfl1 <- 3 # top n max days
  ylab <- expression(paste("Observed, m"^"3" * "/s"))
  xlab <- expression(paste("Simulated, m"^"3" * "/s"))

  meansdmax <- function(md, name, ...) {
    vm <- validate(md, cex = cex1, xlab = xlab, ylab = ylab)
    title(main = paste0("Mean discharge,\n ", name), line = 2)
    vs <- validate(md, cex = cex1, type = "sd", xlab = xlab, ylab = ylab)
    title(main = paste0("SD discharge,\n ", name), line = 2)
    ve <- validate(md, cex = cex1, type = "fl", xlab = xlab, ylab = ylab)
    title(main = paste0("Max. discharge,\n ", name), line = 2)
    return(list(vm, vs, ve))
  }

  par(mfrow = c(3, 3), mar = c(5, 5, 4, 2))
  results_nn <- meansdmax(md1, "default model")
  results_rdw1 <- meansdmax(md2, "RDW1 model")
  results_rdw2 <- meansdmax(md3, "RDW2 model")
  dev.off()
}


## Alphas ======================================================================

adf <- data.frame(
  region = unlist(lapply(1:6, function(x) rep(x, 12))),
  month = rep(1:12, 6),
  x = unlist(alphas_x),
  q = unlist(alphas_q),
  q0 = unlist(alphas_q0)
) # data frame of alphas
idf <- adf[, 1:2]
idf$I <- unlist(lapply(1:nrow(idf), function(x) {
  cat(idf$region[x], idf$month[x], "\n")
  p90 <- pcp90[pcp90$region == idf$region[x] & pcp90$mon == idf$month[x], ]
  return(sd(p90$obs))
}))
idf <- tapply(idf$I, list(idf$region, idf$mon), function(x) x)
pal <- c(
  "#0a2463", "#58a4b0", "#896a67",
  "#0c7c59", "#f6bd60", "#5a0066"
)

if (redo_figures) {
  png(
    filename = "figures/fig3.png", width = 450, height = 700,
    pointsize = 15, units = "px"
  )
  par(mfrow = c(3, 1), mar = c(4, 5, 3, 2))
  for (i in 1:6) {
    if (i == 1) {
      plot(1:12, adf$x[adf$region == i],
        type = "l", col = pal[i],
        lwd = 2, ylim = c(0, 16), xlab = "",
        ylab = expression(paste(alpha, " for occurrence, RDW2"))
      )
    } else {
      lines(1:12, adf$x[adf$region == i], col = pal[i], lwd = 2)
    }
  }
  abline(v = c(4, 10))
  for (i in 1:6) {
    if (i == 1) {
      plot(1:12, adf$q[adf$region == i],
        type = "l", col = pal[i],
        lwd = 2, ylim = c(0, 11), xlab = "",
        ylab = expression(paste(alpha, " for quantity, RDW2"))
      )
    } else {
      lines(1:12, adf$q[adf$region == i], col = pal[i], lwd = 2)
    }
  }
  abline(v = c(4, 10))
  legend(
    x = 5.5, y = 11, legend = c(
      "Region 1", "Region 2", "Region 3",
      "Region 4", "Region 5", "Region 6"
    ),
    lty = 1, col = pal, lwd = 2, box.lty = 0
  )
  for (i in 1:6) {
    if (i == 1) {
      plot(1:12, adf$q0[adf$region == i],
        type = "l", col = pal[i],
        lwd = 2, ylim = c(0, 20), xlab = "Month",
        ylab = expression(paste(alpha, " for quantity, RDW1"))
      )
    } else {
      lines(1:12, adf$q0[adf$region == i], col = pal[i], lwd = 2)
    }
  }
  abline(v = c(4, 10)) ### whole thing should just be a function

  dev.off()
}

### Within-basin analysis ======================================================

b_nse <- function(md) {
  md <- md[md$sub %in% names(upst), ]
  by(md, md$sub, function(x) {
    qobs <- tapply(x$Qo, list(x$mon, x$year), mean)
    qsim <- tapply(x$Qs, list(x$mon, x$year), mean)
    return(NSE(as.vector(qsim), as.vector(qobs)))
  })
}

bn1 <- b_nse(md1) # default model
bn2 <- b_nse(md2) # rdw1
bn3 <- b_nse(md3) # rdw2

### Predictors of failure ======================================================

pof <- data.frame(
  sub = names(bn3),
  nse = as.vector(bn3)
)

# adding elevation
pof$elev <- sapply(pof$sub, function(sub) {
  return(hyd_st$elev[which(hyd_st$sub == sub)])
})

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
nsub <- readRDS("data/nsubs.RDS")
pof$nsub <- unlist(sapply(pof$sub, function(s) {
  if (s == 981) s <- 94
  nsub$nsubs[which(nsub$Subbasin == s)]
})) # watershed size

# adding density of stations per region (no. gauges/no. subs)
subs_per_region <- tapply(wshed$region, wshed$region, length)
pof$dens <- pof$ngauge / subs_per_region[pof$region]

# adding SWAT generated SD
pof$sd <- sapply(pof$sub, function(s) {
  sd(md3$Qs[which(md3$sub == s)], na.rm = T)
})

# SWAT generated mean
pof$mean <- sapply(pof$sub, function(s) {
  mean(md2$Qs[which(md2$sub == s)], na.rm = T)
})

# adding whether a station exists in that sub 1/0
pof$pstat <- sapply(pof$sub, function(s) {
  return(as.numeric(s %in% met_st$sub))
})

# model of all variables
pof0 <- pof ## backup
xvars <- c("elev", "ngauge", "down", "nsub", "sd", "pstat")
for (s in xvars) {
  pof[, s] <- as.vector(scale(pof[, s]))
}
form1 <- paste("nse ~ ", paste(xvars, collapse = "+"))
model1 <- lm(form1, data = pof)
summary(model1)

k <- step(model1, scope = ~ . + .^2 + .^3, direction = "forward", trace = 1)
summary(k)
AIC(k)

pof_model <- summary(k)$coefficients
row.names(pof_model) <- c(
  "Intercept", "1. Elevation",
  "2. Total number of gauges in region",
  "3. Number of downstream subbasins",
  "4. Size of watershed (# subbasins)",
  "5. Simulated standard deviation of discharge",
  "6. Presence of rain gauge within subbasin",
  "Interaction term 3*4",
  "Interaction term 4*5",
  "Interaction term 1*3",
  "Interaction term 2*5"
)

write.csv(pof_model, "figures/table3.csv")
