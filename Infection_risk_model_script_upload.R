#Infection risk model for STB in UK wheat
#Calculates hourly germination and growth risk of Z.tritici ascospores and pycnidiospores at spatial resolution of 0.1° x 0.1°
#Adaptation of model from Chaloner et al. (2019) doi: 10.1098/rstb.2018.0266
#weather vector inputs: temp_by_pixel, wet_by_pixel, precip_by_pixel
#UK growing season runs from Oct-July
#run for all pixels across one growing season


########################### Load weather data -----------------------------------------------
# Load libraries
library(ggplot2); library(R.utils); library(accelerometry); library(abind); library(terra); library(conflicted);library(future.apply)

setwd("~")

# Format raster data ------------------------------------------------------
weather_raster <- rast(".grib") #change this only - combined weather raster for chosen time period
#Layer 1 = temp
#Layer 2 = wetness
#Layer 3 = precip

#temperature
temp_layer <- weather_raster[[seq(1, nlyr(weather_raster), by = 3)]]
writeRaster(temp_layer, "temp_layer.tif", overwrite=TRUE)

#wetness
wet_layer <- weather_raster[[seq(2, nlyr(weather_raster), by = 3)]]
writeRaster(wet_layer, "wet_layer.tif", overwrite=TRUE)

#precip
precip_layer <- weather_raster[[seq(3, nlyr(weather_raster), by = 3)]]
writeRaster(precip_layer, "precip_layer.tif", overwrite=TRUE)

# Load in data ------------------------------------------------------------
#load for all pixels (same geographical area for each)
temp_gs_rast <- rast("temp_layer.tif")  
wet_gs_rast <- rast("wet_layer.tif")
precip_gs_rast <-rast("precip_layer.tif")

#convert to vectors 
temp_gs <- as.vector(temp_gs_rast) #takes values from top-left, across row, down to next row etc until at bottom right
wet_gs <- as.vector(wet_gs_rast)
precip_gs <-as.vector(precip_gs_rast)

temp_gs <- temp_gs - 273.15 #temp is in Kelvin to need to convert to degrees C (K-273.15 = C)

#split into vectors for each pixel
chunk_size <- ncell(temp_gs_rast) #define size of chunks to split vectors into
temp_cells_matrix <- matrix(temp_gs, ncol = chunk_size, byrow = TRUE) #number of cols = ncell(), fill row-wise
temp_by_pixel <- split(temp_cells_matrix, col(temp_cells_matrix))  # Split into list of vectors, each is one pixel across entire growing season

wet_cells_matrix <- matrix(wet_gs, ncol = chunk_size, byrow = TRUE) 
wet_by_pixel <- split(wet_cells_matrix, col(wet_cells_matrix)) 

precip_cells_matrix <- matrix(precip_gs, ncol = chunk_size, byrow = TRUE)
precip_by_pixel <- split(precip_cells_matrix, col(precip_cells_matrix)) 

#specify tmin, topt, tmax, α (a) and γ (g), for germination OR growth depending on which iteration of the germ function you are running 

#tmin <-  9.92; topt <-  19.1 ; tmax <- 32.2; a <- 58.5; g <- 1.3 #germination
tmin <-  -10.35; topt <-  17.17; tmax <- 19.77; a <- 189; g <- 2.2 #growth

############################### Thermal death survival function - the act of drying (hour 1 of dry period) ---------------------------
survival_function_dry <- function(x, c = 1.47885, m = 0.16889){ # c and m calculated from regression analysis of percentage Z tritici thermal death, x is temperature
  prop_dead_1h <- ((m*x + c) / 100) # turn death into a proportion (/ 100)
  prop_dead_1h[prop_dead_1h < 0] <- 0 # if regression predicts negative death, change to 0 (cannot have negative death - this is at very cold (< -8 oC) temperatures)
  prop_dead_1h <- prop_dead_1h + 0.3612 # 0.3612 calculated as the average death due to the act of drying (as opposed to the temperature-dependent death during the dry period)
  surv <- 1 - prop_dead_1h # proportion survival is 1 - death
  surv
}
print(head((survival_function_dry(temp_by_pixel[[1]], c = 1.47885, m = 0.16889))))

############################### Thermal death survival function - all hours of dry after hour 1 ----------------------------------
# Same as above, excluding the 0.3612 additional death
survival_function_dry_continous <- function(x, c = 1.47885, m = 0.16889){
  prop_dead_1h <- ((m*x + c) / 100)
  prop_dead_1h[prop_dead_1h < 0] <- 0 
  surv <- 1 - prop_dead_1h
  surv
}
print(head((survival_function_dry_continous(temp_by_pixel[[1]], c = 1.47885, m = 0.16889))))
# No death under wet conditions

############################### Beta function -----------------------------------------------------------
# Beta function to calculate temperature dependent rate (0, 1) Magarey et al., (2005) (https://apsjournals.apsnet.org/doi/abs/10.1094/PHYTO-95-0092)
rt <- function(x, tmin, topt, tmax){
  r <- ((tmax - x)/(tmax - topt))*((x - tmin)/(topt - tmin))^((topt-tmin)/(tmax - topt))
  r[r < 0] <- 0; r[is.na(r)] <- 0 # Where rates are calculated as <0, set to 0 (this occurs <tmin and >tmax).
  r
}

############################### Timestamps --------------------------------------------------------------
#Retrieve timestamps from SpatRaster
timestamps <- time(temp_gs_rast)#this raster has the times stored in it

# Convert timestamps to POSIXlt format
times <- as.POSIXlt(timestamps)

# Extract the year and month from the hourly timestamps
years <- times$year + 1900  # Adding 1900 to get the correct year as stores as years since 1900
months <- times$mon + 1     # Months are 0-indexed, so add 1

years_in_season <- unique(years) #should be two years per season (Oct-July)

# Define the years `j` and `j+1` 
year_j <- years_in_season[1]
year_j_plus_1 <- years_in_season[2]

#store hours in the season as a vector
temp_first_pixel <- as.numeric(temp_gs_rast[1,1,]) - 273.15 #K-273.15 = C
hour_season <- 1:length(temp_first_pixel)

################################################## Functions for spore germ and growth used in germ ------------------------

############### Ascospore germ ----------------------------------------------------------
calculate_ascospore_burden <- function(n, t, h, Fm, gs_length) {
  ascospore_burden <- (n*t^2*exp(-h*t)) / max(n*t^2*exp(-h*t)) # seasonal ascospore burden determined from equation provided by Kitchen et al., (2016), with modified parameters (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0161887)
  ascospore_burden_matrix <- matrix(ascospore_burden[1:gs_length], nc = gs_length, nr = gs_length, byrow = TRUE)
  Fm <- Fm * ascospore_burden_matrix # multiply germination and ascospore cohort size matrices together. I.e. at peak winter cohort size = 1, so germination is * 1.
  Fm
}

# Pycnidiospore germ ------------------------------------------------------
#Subset hours during october and november to restrict pycnidiospores so that they cannot land during this period.
#Still include these periods in model so that dimensions of ascospores and pycnidiospores are the same
#Make an index for year j, oct and nov to use in calculate_pyc_burden function
index_oct_nov <- (years == year_j) & (months %in% c(10, 11))
cts_october_nov_interest <- temp_gs_rast[index_oct_nov] # cts_october_nov will be used to restrict pycnidiospore germination to 1st December onwards

#relies on precipitation (spores are released from pycnidia and land on available leaves only when it rains)
calculate_pyc_burden <- function(precip_by_pixel, gs_length, Fm, cts_october_nov_interest, tmp) {
  precip_by_pixel[precip_by_pixel > 0] <- 1 # convert precip to binary
  precip_by_pixel <- precip_by_pixel[1:(length(precip_by_pixel) - 1)] # removes last element from vector
  precip_status_cohort <- matrix(precip_by_pixel[1:length(gs_length)], nc = gs_length, nr = gs_length, byrow = FALSE) # create precip matrix. byrow = F because rows represent hours so fill column-wise
  Fm <- Fm * precip_status_cohort # Germination only occurs during periods of precip
  
  # Remove pycnidia landing on leaf in oct or nov (we assume that there is a 2 month delay, the approximate latent period for infection from Bernard et al., (2013) for average temperature experienced during october).
  hours_to_remove_pyc <- length(cts_october_nov_interest) # length of period
  PYC_REMOVE <- rep(0, hours_to_remove_pyc) # cohort size of 0
  PYC_INCLUDE <- rep(1, length(tmp) - hours_to_remove_pyc - 1) # relative cohort size of 1 for 1st December onwards
  pyc_combine <- c(PYC_REMOVE, PYC_INCLUDE)
  plot(pyc_combine)
  
  # create a matrix by row that removes spore cohorts created during october or november
  delete_octnov_pyc_cohorts <- matrix(pyc_combine[1:(length(tmp) - 1)], nc = (length(tmp) - 1), nr = (length(tmp) - 1))
  Fm <- Fm * delete_octnov_pyc_cohorts # Periods before the 1st December are recalculated as 0 (*0), periods on or after 1st December are unchanged (*1)
  Fm
}

# Ascospore OR pycnidiospore growth  --------------------------------------
# Multiply by cohorts that have germinated - gives fraction of spores available to go on to next stage (growth)
calculate_growth <- function(array_germ_risk_vector, Fm) {
  array_germ_risk_vector <- array_germ_risk_vector[1:(length(array_germ_risk_vector) - 1)] # create a vector of germination for hours 1 to hour N, removes last value from initial vector
  array_germ_risk_vector[1] <- 0 # first vector value is "NA", set to 0. this shifts it all by 1 hour, so what germinated between hour 1 and 2 can now grow between hours 2 and 3
  Fd_germination <- matrix(array_germ_risk_vector[1:length(array_germ_risk_vector)], nc = length(array_germ_risk_vector), nr = length(array_germ_risk_vector), byrow = TRUE)# array_germ_risk_vector is from germination
  array_germ_risk_vector <- NULL
  Fm <- Fm * Fd_germination #corrects growth cohort size to be dependent on germination having to have occurred 
  Fm
}


################################################# Germ function -----------------------------------------------------------

germ <- function(canopy.wetness, tmp, hour, tmin, topt, tmax, a, g, array_germ_risk_vector, precip_by_pixel, cts_october_nov_interest){
  # canopy.wetness = wetness data
  # tmp = temperature data
  # hour = list of all hours
  # tmin, topt, tmax = beta parameters for temperature dependent rate 
  # a, g (α and γ) = scale and shape parameter for Weibull distribution 
  # array_germ_risk_vector = germination values - same str as canopy.wetness and tmp, vector of rt function (germ_rate)
  # precip_by_pixel = precip data
  # cts_october_nov_interest - subset of hours for october and november
  
  if(all(is.na(canopy.wetness))) return(NA) else # exit the function if all values of canopy.wetness are NA, i.e. because pixel is in the sea
    canopy.wetness[canopy.wetness > 0] <- 1 # convert wetness to binary (is it wet or not)
  
  h <-1:length(canopy.wetness)#a list of all hour values in a growing season
  print(head(h)); print(tail(h))
  ct <- tmp # rename - temperature
  wt <- canopy.wetness # rename - wetness
  
  # Temperature matrix filled with vectors
  gs_length <<- length(temp_first_pixel)-1 # gs_length = length of growth season - 1 because diff(ct) works pairwise (so you lose one)
  ti_temp <- matrix(ct[1:gs_length]+diff(ct)/2, nc = gs_length, nr = gs_length, byrow = F) # matrix of temperature mid-points between readings
  rate <- rt(ti_temp, tmin, topt, tmax) # calculates the relative rates (0,1) along the beta function of each temperature
  im_temp <- matrix(1:gs_length, nc = gs_length, nr =  gs_length) - rep((1:gs_length)-1,each=gs_length) #shifts each row to the matrix 1 along
  im0 <- im_temp-1 # takes 1 off each value of the matrix
  print(head(ti_temp)) 
  print(head(im0)) 
  
  #Wetness matrix filled with vectors
  ti_wet <- matrix(wt[1:gs_length], nc = gs_length, nr = gs_length, byrow = F) # create a matrix for wetness, where last value is removed (to match ti_temp) 
  # midpoints not used as with temperature, values are binary (dry or wet, 0 or 1)
  print(head(ti_wet))
  
  # Hazard function
  Hm_temp <- rate*((im_temp/a)^g - (im0/a)^g) # Probability of germination/growth multiplied by temperature dependent rate
  dimnames(Hm_temp) <- list("i" = 1:gs_length, "j" = 1:gs_length) #i is referring to the rows of the matrix, j is the columns
  Hm_temp[ti_wet == 0] <- 0 # if dry then no probability of germination/growth
  Hm_temp[im0 < 0] <- 0 # make upper triangle 0, we are interested in the lower part of the triangle (where i > j)
  #zeroing out the upper triangle of the matrix, focusing on the relevant lower part (which corresponds to past or present conditions where spore cohorts have actually landed)
  Hm <- Hm_temp # instantaneous hazard at a given hour, given all hours that have passed
  print(head(Hm))
  
  prop_survive_dry <- survival_function_dry(ti_temp) # temperature dependent survival for dry hour, when wet previous hour
  prop_dead_dry <- 1 - prop_survive_dry 
  print(head(prop_dead_dry))
  
  perc_survive_dry_continous <- survival_function_dry_continous(ti_temp) # temperature dependent survival for dry hour, when it was dry previous hour (no additional death from drying event)
  print(head(perc_survive_dry_continous))
  
  #Turn all the top right hand corner (upper triangle is future cohorts, haven't landed yet)
  #to 1s (survival) and 0s (death) - ensures cumprod equation is not effected by cohorts that have not yet landed on a leaf surface
  prop_survive_dry[im0 < 0] <- 1 # survival = 1, spores have not yet landed
  perc_survive_dry_continous[im0 < 0] <- 1 # survival = 1, spores have not yet landed
  prop_dead_dry[im0 < 0] <- 0 # death = 0, spores have not yet landed
  
  prop_survive_wet <- prop_survive_dry; prop_survive_wet[,] <- 1 # survival under wet conditions = 1 (no death under wet conditions)
  prop_dead_wet <- prop_dead_dry; prop_dead_wet[,] <- 0 # death under wet conditions = 0 (no death under wet conditions)
  
  # Determine starting matrix depending on whether it is wet or dry in hour 1
  if(ti_wet[1,1] == 1){
    av <- prop_survive_wet # number available to transition in hour 1
    perc_dead <- prop_dead_wet #no death under wet conditions
  } else {
    av <- prop_survive_dry # number available to transition in hour 1
    perc_dead <- prop_dead_dry
  }
  
  # Determine hazard and thermal death in subsequent hours. Dependent on previous hour
  Fm <- ti_temp * 0 # create a new empty hazard matrix, with the same dimensions as ti_temp/ti_wet
  Fm[1,] <- av[1,]*(1-exp(-Hm[1,])) # Fm[1,] = number that transition to next stage of life cycle, given that x died
  
  i2 <- vector()
  N <- nrow(ti_temp)
  
  # Each hour will have a hazard and thermal death, dependent on whether it is wet or dry + temperature
  # Loop to process each hour in the growth season
  for(i in 2:N){ # from 2nd hour in growth season to Nth hour
    i2 <- i - 1 # id for previous hour
    if(ti_wet[i,1] == 1){ # if  wet in current hour
      av[i,] <- (av[i2,] - Fm[i2,]) # the number available to germinate in hour i = (available to germ in i2 - those that did germ in i2) (* 1 (because no death under wet conditions))
    } else { # current hour is dry
      if(ti_wet[i2,1] == 0){ # if dry in previous hour
        av[i,] <- (av[i2,] - Fm[i2,]) * perc_survive_dry_continous[i,] # the number available to germinate in hour i = (available to germ in i2 - those that did germ in i2) * the fraction that survived dryness in that hour (calculated from function above) i.e. 0.6 i.e. 60%
      } else { # previous hour was wet (the act of drying has occurred)
        av[i,] <- (av[i2,] - Fm[i2,]) * prop_survive_dry[i,] # same as above but uses function that includes act of drying in death
      }
    }
    perc_dead[i,] <- (av[i2,] - Fm[i2,]) - av[i,] # number that didn't germinate in previous hour - those now available to germinate
    Fm[i,] <- av[i,]*(1-exp(-Hm[i,])) # proportion that transition (germinate or grow and penetrate stomata). Calculated from the cohort size available * the temperature dependent instantaneous hazard
  }
  

############ Spore and process selection ---------------------------------------------
  #for following, only run asco germ, pyc germ OR growth at once

  ###Ascospore burden
  n <- 1
  t <- seq(0, gs_length, by = 1)
  h <- (0.0035 / 2.2) # parameters chosen to provide an approximate fit of peak during winter and tail to summer

# !!! Ascospore germination -----------------------------------------------
  #Fm <- calculate_ascospore_burden(n, t, h, Fm, gs_length)
  

  
# !!! Pycnidiospore germination -------------------------------------------
  # relies on precipitation (spores are released from pycnidia and land on available leaves only when it rains)
  #Fm <- calculate_pyc_burden(precip_by_pixel, gs_length, Fm, cts_october_nov_interest, tmp)
 
  
# !!! Ascospore OR pycnidiospore growth -----------------------------------
  Fm <- calculate_growth(array_germ_risk_vector, Fm)
  
  
# Cumulative
  Fm_cum <- apply(Fm, 2, cumsum)# cumulative germination/growth across intervals
  print(head(Fm_cum))
  Fd <- diff(c(0, rowSums(Fm_cum))) 
  print(head(Fd))# total germination / growth in a specific hour
  Fd <- c(NA, Fd) # NA is added for hour 1. This is because the number of hours in the model is N - 1 growing season, because temperature is interpolated. Adding "NA" to hour 1 corrects for this, so the model output has the same dimensions as the climate data.
  Fd
  #Fm_cum is the cumulative risk for each cohort (as a matrix), Fd is the cumulative total growth (combined cohorts, as a vector)
}
##END OF GERM

############### Run Model  -------------------------------------------------------------------------
#iterates over list elements in parallel
# Use fewer cores to prevent crashes (max 4 workers)
plan(multisession, workers = min(4, availableCores() - 1)) 

# Run only on the first 100 elements to check for memory issues before running all
#test_indices <- seq_along(wet_by_pixel)[1]
#change wet_by_puxel (below) to test_indices if testing

#Run Model
All_summary_growth_risk <- future_lapply(seq_along(wet_by_pixel), function(i) {
  
  
  germ_result <- rt(temp_by_pixel[[i]], tmin, topt, tmax)
  
  germ(
    canopy.wetness = wet_by_pixel[[i]], 
    tmp = temp_by_pixel[[i]], 
    hour = hour_season[[i]], 
    tmin = tmin, 
    topt = topt, 
    tmax = tmax, 
    a = a, 
    g = g, 
    array_germ_risk_vector = germ_result,  
    precip_by_pixel = precip_by_pixel[[i]],  
    cts_october_nov_interest = cts_october_nov_interest[[i]]
  )
})


#remove NA value at start of each element (added earlier in code to make lengths correct)
cleaned_summary_growth_risk <- lapply(All_summary_growth_risk, na.omit)
 
plot(cleaned_summary_growth_risk[[1]], type = "l", xlab = "Time", ylab = "Pyc Germ Risk")

# Save model output as risk across the growth season, separated by pixel 
write.csv(cleaned_summary_growth_risk, file = paste0("Growth_17711", year_j, ".csv"))
