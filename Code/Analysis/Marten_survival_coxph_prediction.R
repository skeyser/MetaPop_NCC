## ----------------------------------------------------------
##
## Script name: Survival Model Validation
##
## Script purpose:
##
## Author: Spencer R Keyser & Matthew M. Smith
##
## Date Created: 2024-01-01
##
## Email: skeyser@wisc.edu
##
## Github: https://github.com/skeyser
##
## -----------------------------------------------------------
##
## Notes:
##
##
## -----------------------------------------------------------

## Defaults
options(scipen = 6, digits = 4)

## -----------------------------------------------------------

## Package Loading
library("pec")
library("blockCV")
library("survAUC")
library("SurvMetrics")
library("ggplot2")
library("dplyr")
library("survival")
library("sf")
library("mapview")
library("here")

#############################################################################################################################################
#################### Read in data and create sf object for blockCV ##########################################################################
#############################################################################################################################################
data.off <- read.csv(here("Data/Marten_Survival/Snow_off_coxph_input.csv"))
data.on <- read.csv(here("Data/Marten_Survival/Snow_on_coxph_input.csv"))

data.off.sf <- st_as_sf(data.off, coords = c("easting", "northing"), crs = 26915) ## NAD83, UTM Zone 15N
mapview(data.off.sf) # check locations

data.on.sf <- st_as_sf(data.on, coords = c("easting", "northing"), crs = 26915) ## NAD83, UTM Zone 15N
mapview(data.on.sf) # check locations

#############################################################################################################################################
#################### Block CV ###############################################################################################################
#############################################################################################################################################
folds <- 5

sc.off <- cv_cluster(x = data.off.sf, k = folds)
# now plot the created folds
cv_plot(cv = sc.off, # a blockCV object
        x = data.off.sf, # sample points
        points_alpha = 0.5,
        nrow = 2)

sc.on <- cv_cluster(x = data.on.sf, k = folds)
# now plot the created folds
cv_plot(cv = sc.on, # a blockCV object
        x = data.on.sf, # sample points
        points_alpha = 0.5,
        nrow = 2)

## List for the partitions
sp.cv.off <- vector(mode = "list", length = folds)

for(i in 1:folds){
  sp.off.train <- data.off.sf[sc.off$folds_list[[i]][[1]],]
  sp.off.train$Type <- "Train"
  sp.off.test <- data.off.sf[sc.off$folds_list[[i]][[2]],]
  sp.off.test$Type <- "Test"

  sp.cv.off[[i]] <- list(sp.off.train, sp.off.test)

}

sp.cv.on <- vector(mode = "list", length = folds)

for(i in 1:folds){
  sp.on.train <- data.on.sf[sc.on$folds_list[[i]][[1]],]
  sp.on.train$Type <- "Train"
  sp.on.test <- data.on.sf[sc.on$folds_list[[i]][[2]],]
  sp.on.test$Type <- "Test"

  sp.cv.on[[i]] <- list(sp.on.train, sp.on.test)

}




#############################################################################################################################################
#################### 80/20 testing ##########################################################################################################
#############################################################################################################################################
#Create training and testing set
on.train <- data.on %>% sample_frac(0.80) #Create test set
on.test  <- anti_join(data.on, on.train, by = 'ID')
#Create training and testing set
off.train <- data.off %>% sample_frac(0.80) #Create test set
off.test  <- anti_join(data.off, off.train, by = 'ID')

### AUC calculations for snow on ###
on.train.fit <- survival::coxph(survival::Surv(exit,event) ~ snow, x=TRUE, data=on.train)
summary(on.train.fit)

on.lpnew <- predict(on.train.fit, newdata = on.test)
on.lp <- predict(on.train.fit)

on.surv.rsp <- survival::Surv(on.train$exit, on.train$event)
on.surv.rsp.new <- survival::Surv(on.test$exit, on.test$event)

on.times <- seq(from = 1, to = max(data.on$exit), by = 3)

# Concordance index
on.Cstat <- UnoC(on.surv.rsp, on.surv.rsp.new, on.lpnew)
on.Cstat

# AUC estimators
on.AUC.Uno <- AUC.uno(on.surv.rsp, on.surv.rsp.new, on.lpnew, on.times)
on.AUC.Uno
names(on.AUC.Uno)
on.AUC.Uno$iauc
plot(on.AUC.Uno)

on.AUC.sh <- AUC.sh(on.surv.rsp, on.surv.rsp.new, on.lp, on.lpnew, on.times)
on.AUC.sh
names(on.AUC.sh)
on.AUC.sh$iauc
plot(on.AUC.sh)

on.AUC.cd <- AUC.cd(on.surv.rsp, on.surv.rsp.new, on.lp, on.lpnew, on.times)
on.AUC.cd
names(on.AUC.cd)
on.AUC.cd$iauc
plot(on.AUC.cd)

### AUC calculations for snow off ###
off.train.fit <- survival::coxph(survival::Surv(exit,event) ~ entropy, x=TRUE, data=off.train)
summary(off.train.fit)

off.lpnew <- predict(off.train.fit, newdata = off.test)
off.lp <- predict(off.train.fit)

off.surv.rsp <- survival::Surv(off.train$exit, off.train$event)
off.surv.rsp.new <- survival::Surv(off.test$exit, off.test$event)

off.times <- seq(from = 1, to = max(data.off$exit), by = 3)

# Concordance index
off.Cstat <- UnoC(off.surv.rsp, off.surv.rsp.new, off.lpnew)
off.Cstat

# AUC estimators
off.AUC.Uno <- AUC.uno(off.surv.rsp, off.surv.rsp.new, off.lpnew, off.times)
off.AUC.Uno
names(off.AUC.Uno)
off.AUC.Uno$iauc
plot(off.AUC.Uno)

off.AUC.sh <- AUC.sh(off.surv.rsp, off.surv.rsp.new, off.lp, off.lpnew, off.times)
off.AUC.sh
names(off.AUC.sh)
off.AUC.sh$iauc
plot(off.AUC.sh)

off.AUC.cd <- AUC.cd(off.surv.rsp, off.surv.rsp.new, off.lp, off.lpnew, off.times)
off.AUC.cd
names(off.AUC.cd)
off.AUC.cd$iauc
plot(off.AUC.cd)

## -----------------------------------------------------------
##
## Begin Section: Fitting models with the new spatCV folds
##
## -----------------------------------------------------------

## Cox survival models for Snow off period
### AUC calculations for snow off ###
sp.cv.off[[1]]
sp.fit.list <- vector(mode = "list", length = folds)
for(dat in 1:length(sp.cv.off)){
  fit.dat <- as.data.frame(sp.cv.off[[dat]][1])
  new.dat <- as.data.frame(sp.cv.off[[dat]][2])
  cox.fit <- survival::coxph(survival::Surv(exit,event) ~ entropy, x=TRUE, data=fit.dat)

  ## Predict to the training data
  off.lp <- predict(cox.fit)
  off.lpnew <- predict(cox.fit, newdata = new.dat)

  ##
  off.surv.rsp <- survival::Surv(fit.dat$exit, fit.dat$event)
  off.surv.rsp.new <- survival::Surv(new.dat$exit, new.dat$event)

  # Concordance index
  off.Cstat <- UnoC(off.surv.rsp, off.surv.rsp.new, off.lpnew)
  off.Cstat

  # AUC estimators
  off.times <- seq(from = 1, to = max(fit.dat$exit), by = 3)
  off.AUC.Uno <- AUC.uno(off.surv.rsp, off.surv.rsp.new, off.lpnew, off.times)
  off.AUC.sh <- AUC.sh(off.surv.rsp, off.surv.rsp.new, off.lp, off.lpnew, off.times)
  off.AUC.cd <- AUC.cd(off.surv.rsp, off.surv.rsp.new, off.lp, off.lpnew, off.times)
  p <- plot(off.AUC.Uno)

  ## Return all the objects
  l.tmp <- list(Fold = dat,
                Fit = cox.fit,
                PredCox = off.lpnew,
                Cstat = off.Cstat,
                AUC.uno = off.AUC.Uno,
                AUC.sh = off.AUC.sh,
                AUC.cd = off.AUC.cd
                )

  sp.fit.list[[dat]] <- l.tmp

}

str(sp.fit.list)

## View the first fold
mean(c(sp.fit.list[[1]]$AUC.uno$iauc,
sp.fit.list[[2]]$AUC.uno$iauc,
sp.fit.list[[3]]$AUC.uno$iauc,
sp.fit.list[[4]]$AUC.uno$iauc,
sp.fit.list[[5]]$AUC.uno$iauc
))

mean(c(sp.fit.list[[1]]$AUC.sh$iauc,
       sp.fit.list[[2]]$AUC.sh$iauc,
       sp.fit.list[[3]]$AUC.sh$iauc,
       sp.fit.list[[4]]$AUC.sh$iauc,
       sp.fit.list[[5]]$AUC.sh$iauc
))

mean(c(sp.fit.list[[1]]$Cstat,
     sp.fit.list[[2]]$Cstat,
     sp.fit.list[[3]]$Cstat,
     sp.fit.list[[4]]$Cstat,
     sp.fit.list[[5]]$Cstat
))

## -----------------------------------------------------------
##
## Begin Section: Snow on
##
## -----------------------------------------------------------
sp.cv.on[[1]]
sp.fiton.list <- vector(mode = "list", length = folds)
for(dat in 1:length(sp.cv.off)){
  fit.dat <- as.data.frame(sp.cv.on[[dat]][1])
  new.dat <- as.data.frame(sp.cv.on[[dat]][2])
  cox.fit <- survival::coxph(survival::Surv(exit,event) ~ snow, x=TRUE, data=fit.dat)

  ## Predict to the training data
  on.lp <- predict(cox.fit)
  on.lpnew <- predict(cox.fit, newdata = new.dat)

  ##
  on.surv.rsp <- survival::Surv(fit.dat$exit, fit.dat$event)
  on.surv.rsp.new <- survival::Surv(new.dat$exit, new.dat$event)

  # Concordance index
  on.Cstat <- UnoC(on.surv.rsp, on.surv.rsp.new, on.lpnew)
  on.Cstat

  # AUC estimators
  on.times <- seq(from = 1, to = max(fit.dat$exit), by = 3)
  on.AUC.Uno <- AUC.uno(on.surv.rsp, on.surv.rsp.new, on.lpnew, on.times)
  on.AUC.sh <- AUC.sh(on.surv.rsp, on.surv.rsp.new, on.lp, on.lpnew, on.times)
  on.AUC.cd <- AUC.cd(on.surv.rsp, on.surv.rsp.new, on.lp, on.lpnew, on.times)
  p <- plot(on.AUC.Uno)

  ## Return all the objects
  l.tmp <- list(Fold = dat,
                Fit = cox.fit,
                PredCox = on.lpnew,
                Cstat = on.Cstat,
                AUC.uno = on.AUC.Uno,
                AUC.sh = on.AUC.sh,
                AUC.cd = on.AUC.cd
  )

  sp.fiton.list[[dat]] <- l.tmp

}

str(sp.fiton.list)

## View the first fold
mean(c(sp.fiton.list[[1]]$AUC.uno$iauc,
       sp.fiton.list[[2]]$AUC.uno$iauc,
       sp.fiton.list[[3]]$AUC.uno$iauc,
       sp.fiton.list[[4]]$AUC.uno$iauc,
       sp.fiton.list[[5]]$AUC.uno$iauc
))

mean(c(sp.fiton.list[[1]]$AUC.sh$iauc,
       sp.fiton.list[[2]]$AUC.sh$iauc,
       sp.fiton.list[[3]]$AUC.sh$iauc,
       sp.fiton.list[[4]]$AUC.sh$iauc,
       sp.fiton.list[[5]]$AUC.sh$iauc
))

mean(c(sp.fiton.list[[1]]$Cstat,
       sp.fiton.list[[2]]$Cstat,
       sp.fiton.list[[3]]$Cstat,
       sp.fiton.list[[4]]$Cstat,
       sp.fiton.list[[5]]$Cstat
))



summary(off.train.fit)

off.lpnew <- predict(off.train.fit, newdata = off.test)
off.lp <- predict(off.train.fit)

off.surv.rsp <- survival::Surv(off.train$exit, off.train$event)
off.surv.rsp.new <- survival::Surv(off.test$exit, off.test$event)

off.times <- seq(from = 1, to = max(data.off$exit), by = 3)

# Concordance index
off.Cstat <- UnoC(off.surv.rsp, off.surv.rsp.new, off.lpnew)
off.Cstat

# AUC estimators
off.AUC.Uno <- AUC.uno(off.surv.rsp, off.surv.rsp.new, off.lpnew, off.times)
off.AUC.Uno
names(off.AUC.Uno)
off.AUC.Uno$iauc
plot(off.AUC.Uno)

off.AUC.sh <- AUC.sh(off.surv.rsp, off.surv.rsp.new, off.lp, off.lpnew, off.times)
off.AUC.sh
names(off.AUC.sh)
off.AUC.sh$iauc
plot(off.AUC.sh)

off.AUC.cd <- AUC.cd(off.surv.rsp, off.surv.rsp.new, off.lp, off.lpnew, off.times)
off.AUC.cd
names(off.AUC.cd)
off.AUC.cd$iauc
plot(off.AUC.cd)


