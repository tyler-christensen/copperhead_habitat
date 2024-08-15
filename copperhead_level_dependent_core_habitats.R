rm(list = ls())

##### Import and clean data #####

library(readxl)

##import data
data <- read_excel("2016_2019_copperhead_tracking_data_shareable.xlsx")

##select variables
names(data)

library(dplyr)

## count missing habitat structure values

cols <- c("CWD", "GND", "WDY", "RCK", "CAN")
colSums(is.na(data[, cols]))

##impute missing values using column means; scale and center variables

for (c in cols){
  col_mean <- mean(data[[c]], na.rm = TRUE)
  data[[c]][is.na(data[[c]])] <- col_mean
  data[[c]] <- as.vector(data[[c]])
}

## create unscaled data columns
data$CWD_unsc <- data$CWD
data$GND_unsc <- data$GND
data$WDY_unsc <- data$WDY
data$RCK_unsc <- data$RCK
data$CAN_unsc <- data$CAN

for (c in cols){
  data[[c]] <- scale(data[[c]])
  data[[c]] <- as.vector(data[[c]])
}

rm(c, col_mean, cols)

## Create variable "snake.season"

library(stringr)

data$snake.season <- str_c(data$year, "_", data$ind_id) ##concatenate year and individual so i = snake*season
names(data)
data$snake.season <- as.factor(data$snake.season) ##convert to factor for list (later)

## Remove snake.seasons with fewer than 15 observations

for(l in (levels(data$snake.season))){ ##k is now each level of snake.season
  test_ID <- l ##"test_ID" takes on the value of l
  if(nrow(dplyr::filter(data, snake.season == test_ID)) < 20) { ##condition: if data$snake.season < 15, then...
    data <- subset(data, data$snake.season != test_ID) ##retain only rows of data where snake.season does not equal test_ID
  }
}

rm(test_ID, l)

length(levels(data$snake.season)) #shows that R retained all original factor levels, messing up later for loops
data$snake.season <- droplevels(data$snake.season) ##drops factor levels without observations from snake.season
length(levels(data$snake.season)) ##check that factor levels were dropped

##### Snake-season statistics #####

ssndf <- data.frame("ind" = "", "snake.ssn" = "", "n" = "", "status" = "")

for(l in (levels(data$snake.season))){
  #l <- levels(data$snake.season)[[2]]
  index <- which(levels(data$snake.season) == l)
  test_df <- data %>% filter(snake.season == l)
  ind <- test_df$ind_id[1]
  obs <- nrow(test_df)
  stat <- test_df$year_repr[[1]]
  new_entry <- c(ind, l, obs, stat)
  ssndf[index,] <- new_entry
}

rm(stat, obs, new_entry, index, ind, test_df)

length(unique(ssndf$ind)) ## number of individuals
range(ssndf$n) ## range in number of observations per individual
table(ssndf$status) ## number of gravid females, males, and nongravid females
mean(as.numeric(ssndf$n)) ## average number of observations per snake-season
sd(as.numeric(ssndf$n)) ## sd in number of observations per snake-season

##### Define coordinate system, convert to UTM #####

library(sp)
library(sf)

sp::coordinates(data) <- ~long+lat
st_crs(data) #shows there is no coordinate system

sp::proj4string(data) <- CRS("+init=epsg:4326") #set current coordinate system (WGS84)
data <- st_as_sf(data)
data <- st_transform(data, CRS("+proj=utm +zone=18 +datum=NAD83 +units=m"))

coords <- st_coordinates(data) #retrieve coordinates
colnames(coords) <- c("long", "lat") #rename coordinate columns
data <- cbind(data, coords) #merge

rm(coords)

##### Create blank raster (grid) to define extents of kernel density UDs #####

library(raster)

##create blank raster
xmin <- round(min(data$long) - 200) #set minimum x extent
xmax <- round(max(data$long) + 200) #max x extent
ymin <- round(min(data$lat) - 200) #min y extent
ymax <- round(max(data$lat) + 200) #max y extent
ext <- extent(xmin, xmax, ymin, ymax) #extent of new raster according to point extents

cols <- xmax - xmin
rows <- ymax - ymin

r <- raster(ext, nrow = rows, ncol = cols) #create blank raster at new extent
r.sp <- as(r, "SpatialPixels")
rm(xmin, xmax, ymin, ymax, rows, l, cols)

##### Individual KDE UDs #####

## Prep for loop

## first, prep new use level variables:
data_UD <- data #new data frame for use level assignments to observations

## Individually-optimized kde UDs
data_UD$pd_ind <- 0 ## individual prob density raster vals. Raster values sum to 1.
#data_UD$vol_ind <- 0 ## individual vol contour (% height) vals. Raster values from 0 to 100. Low values = intense use.
data_UD$ext_ind <- NA

## Population optimized kde UD
data_UD$pd_pop <- 0
#data_UD$vol_pop <- 0
data_UD$ext_pop <- NA

h <- seq(6, 100, by = 1) ##range of possible values of h
bw_df <- data.frame(matrix(ncol= 2, nrow = 0)) ##create an empty data frame to store bandwidths
n <- c("snake.season", "h") #column names for bw_df
colnames(bw_df) <- n #column names for bw_df
idx <- 1 #index

KDE.stats <- data.frame("snake.season" = "", "repr.stat" = "", "home.range" = "", "core.area" = "", "core.prop" = "")

levels(data_UD$snake.season)

library(raster)
library(adehabitatHR)

for(i in levels(data_UD$snake.season)) {
  #i <- "2016_AGCO_03" #i is now the first snake.season (temporarily)
  temp_df <- base::subset(data_UD, snake.season == i) #data frame of observations of ith snake.season
  
  temp_sp <- temp_df %>% dplyr::select(geometry)
  temp_sp <- as(temp_df, "Spatial") #convert to Spatial Points Data Frame for mcp and kde
  bw <- vector() #vector for storing snake bandwidths
  
  ##next, for each (ith) snake.season, optimize bandwidth (j)
  for(j in 1:length(h)) { ##this loop optimizes the bandwidth for the ith individual
    m <- suppressWarnings(mcp(temp_sp, percent = 95, unin = "m", unout = "ha")) # MCP
    m.area <- as.numeric(m$area) ##turns the MCP area output into a single-value numeric vector
    k <- suppressWarnings(kernelUD(temp_sp, h = j, grid = 100, kern = "bivnorm"))
    k.area <- suppressWarnings(as.numeric(kernel.area(k, 95))) ##like mcp area, turns kde area output into a single-value numeric vector
    if (k.area < m.area) { ##if the kernel area is still less than the mcp area, then...
      new_element <- j ## ...create a new element with value j (the last test bandwidth)
      bw[[length(bw) + 1]] <- new_element ##add a new spot to bw (list of test bandwidths) and add new_element to it
    } 
  } 
  
  h_final <- tail(bw, n=1) ##save the last element of bw (the list of test bandwidths) to h_final

  ## create a temporary dataframe (new_obs) with the snake.season and its bandwidth, and add this to bw_df
  new_obs <- data.frame(i, h_final)
  bw_df <- rbind(bw_df, setNames(new_obs, names(bw_df))) ##binds bw_df with new_obs (and renames new_obs columns on the basis of the bw_df columns)
  
  ##### Now that there is a data frame for the ith individual (temp_coordsUTM), and it's h_final, 
  ##### assign each observation for ind_i its UD levels:
  
  ##The code below (1) generates a kde using the optimized bandwidth; (2) runs "getvolumeUD" (getvolumeUD creates a modified UD so that 
  ##the value of each pixel equals the percent volume of the UD at the contour where the pixel is located;
  ##and then(3) converts the kde UD surface to a raster.
  
  ## Generate new kde UD using the optimized bandwidth and create UD object for rasterization
  kde <- kernelUD(temp_sp, h = h_final, grid = r.sp, kern = "bivnorm") ##kernel UD using optimized bandwidth (equivalent to 'final' k of for loop)

  ## rasterize
  kde.rast <- raster(kde)
  sum(values(kde.rast)) ## check that raster correctly sums to 1
  
  ## extract kde values to points
  pd.i <- raster::extract(kde.rast, temp_sp)
  temp_df$pd_ind <- pd.i
  
  ## evaluate which points fall within the upper and lower third of the UD volume
  values <- values(kde.rast) # list of values in kde raster
  values <- sort(values, decreasing = TRUE) # sort descending (highest probability, lowest volume first)
  cumsum <- cumsum(values) # create vector cumulative sum of raster values
  
  upper <- max(which(cumsum < 0.3)) #which prob dens value is at the ith percentile?
  lower <- min(which(cumsum > 0.3)) #which prob dens value is at the ith percentile?
  
  #upper <- max(which(cumsum < 0.5))
  #lower <- min(which(cumsum > 0.5))
  
  temp_df$ext_ind[temp_df$pd_ind > values[upper]] <- 1 #if prob dens is > than this value, 1
  temp_df$ext_ind[temp_df$pd_ind < values[lower]] <- 0 #if prob dens is < than this value, 0
  
  name <- paste(i, ".data", sep = "") #prepare name for random-season data frame
  assign(name, temp_df) #rename i.df
  
  ## populate new entry for KDE.stats
  t.df <- base::subset(data_UD, snake.season == i) ##another temp data frame
  repr <- t.df$year_repr[1]
  ind <- i
  hr <- suppressWarnings(as.numeric(kernel.area(kde, percent = 95)))
  core <- suppressWarnings(as.numeric(kernel.area(kde, percent = 33)))
  core.prop <- core/hr
  new_entry <- c(i, repr, hr, core, core.prop)
  KDE.stats[idx,] <- new_entry
  idx <- idx + 1
  
  print(i)
  
} ##closes main for loop (over m snake.seasons)

KDE.stats
bw_df

##### Core area stats #####

paste("mean core home range =", mean(as.numeric(KDE.stats$core.area)))
paste("SD core home range =", sd(as.numeric(KDE.stats$core.area)))
paste("range of core home range [low, high] =", range(as.numeric(KDE.stats$core.area)))

##### Home range stats #####

G.hr <- KDE.stats %>% filter(repr.stat == "GVD")
G.n <- nrow(G.hr)
G.mean <- mean(as.numeric(G.hr$home.range))
G.sd <- sd(as.numeric(G.hr$home.range))
G.se <- G.sd/sqrt(nrow(G.hr))
G.CI <- G.se*1.96
G.range <- range(as.numeric(G.hr$home.range))

N.hr <- KDE.stats %>% filter(repr.stat == "NGR")
N.n <- nrow(N.hr)
N.mean <- mean(as.numeric(N.hr$home.range))
N.sd <- sd(as.numeric(N.hr$home.range))
N.se <- N.sd/sqrt(nrow(N.hr))
N.CI <- N.se*1.96
N.range <- range(as.numeric(N.hr$home.range))

M.hr <- KDE.stats %>% filter(repr.stat == "MAL")
M.n <- nrow(M.hr)
M.mean <- mean(as.numeric(M.hr$home.range))
M.sd <- sd(as.numeric(M.hr$home.range))
M.se <- M.sd/sqrt(nrow(M.hr))
M.CI <- M.se*1.96
M.range <- range(as.numeric(M.hr$home.range))

KDE.repr.stats <- data.frame("repr" = as.character(), "n" = as.numeric(), "mean" = as.numeric(), "SD" = as.numeric(), "range.min" = as.numeric(), "range.max" = as.numeric())

KDE.repr.stats[1,] <- c("GVD", G.n, G.mean, G.sd, G.range[1], G.range[2])
KDE.repr.stats[2,] <- c("NGR", N.n, N.mean, N.sd, N.range[1], N.range[2])
KDE.repr.stats[3,] <- c("MAL", M.n, M.mean, M.sd, M.range[1], M.range[2])

KDE.repr.stats

rm(G.hr, k, kde, M.hr, N.hr, t.df, vol.i.rast, core, core.prop, G.CI, G.mean, G.n, G.range, G.sd,
   G.se, hr, idx, ind, M.CI, M.mean, M.n, M.range, m, new_obs, temp_df, temp_sp, bw, cond2, h, kde.rast, vol,
   M.sd, M.se, N.CI, h_final, i, j, k.area, m.area, n, new_element, name, ext, pd.i, 
   N.mean, N.n, N.range, N.sd, N.se, new_entry, r_max, repr, vol.i)

##### Prepare for loop for population level UD model #####

## for loop to obtain the population-level bandwidth. The 95% population-level MCP area changed only slightly
## across sets when each set consisted of one randomly selected season per snake. Therefore, to save computing time,
## obtain a single bandwidth to use on all iterations of the population-level for loop.

hseq <- seq(10, 100, by = 1) ##range of possible values of h
pop_bw_df <- data.frame(matrix(ncol= 1, nrow = 0)) ##create an empty data frame to store bandwidths
cols <- c("h") #column names for bw_df
colnames(pop_bw_df) <- cols #column names for bw_df

## Calculate MCP area for entire population
data_UDgeo <- data_UD %>% dplyr::select(geometry)
data_UDsp <- as(data_UDgeo, "Spatial") #convert to Spatial Points Data Frame for mcp and kde

mpop <- suppressWarnings(mcp(data_UDsp, percent = 95, unin = "m", unout = "ha")) # MCP
mpop.area <- as.numeric(mpop$area) ##turns the MCP area output into a single-value numeric vector


for (h in hseq) {
  #h <- hseq[1]
  k <- suppressWarnings(kernelUD(data_UDsp, h = h, grid = 100, kern = "bivnorm"))
  kpop.area <- suppressWarnings(as.numeric(kernel.area(k, 95))) ##like mcp area, turns kde area output into a single-value numeric vector
  if (kpop.area < mpop.area) { ##if the kernel area is still less than the mcp area, then...
    new_element <- h ## ...create a new element with value h (the last test bandwidth)
    pop_bw_df[nrow(pop_bw_df) + 1,1] <- new_element ##add a new spot to bw (list of test bandwidths) and add new_element to it
  } else {
    break
  }
}

pop_bw <- pop_bw_df[nrow(pop_bw_df),1]

rm(new_element, mpop.area, kpop.area, hseq, h, cols, mpop, k, data_UDgeo)

## Create list where each element is a list of 1 - 3 snake-seasons of data

m03 <- list(`2016_AGCO_03.data`, `2017_AGCO_03.data`, `2018_AGCO_03.data`)
m06 <- list(`2017_AGCO_06.data`, `2018_AGCO_06.data`, `2019_AGCO_06.data`) 
m07 <- list(`2017_AGCO_07.data`, `2018_AGCO_07.data`) 
m09 <- list(`2017_AGCO_09.data`, `2018_AGCO_09.data`, `2019_AGCO_09.data`) 
m11 <- list(`2018_AGCO_11.data`, `2019_AGCO_11.data`) 
m17 <- list(`2018_AGCO_17.data`, `2019_AGCO_17.data`) 
m18 <- list(`2018_AGCO_18.data`, `2019_AGCO_18.data`) 
m21 <- list(`2018_AGCO_21.data`, `2019_AGCO_21.data`) 

s.data_UD <- rbind(`2017_AGCO_05.data`, `2017_AGCO_08.data`, `2018_AGCO_10.data`, `2018_AGCO_20.data`,
                   `2018_AGCO_22.data`, `2019_AGCO_27.data`, `2019_AGCO_31.data`, `2019_AGCO_33.data`)

mlists <- list(m03, m06, m07, m09, m11, m17, m18, m21)

#list of individual snake seasons. When a snake had multiple seasons, the element for that individual is a list.

data_UD <- rbind(`2016_AGCO_03.data`, `2017_AGCO_03.data`, `2018_AGCO_03.data`, `2017_AGCO_06.data`, 
   `2018_AGCO_06.data`, `2019_AGCO_06.data`, `2017_AGCO_07.data`, `2018_AGCO_07.data`,
   `2017_AGCO_09.data`, `2018_AGCO_09.data`, `2019_AGCO_09.data`, `2018_AGCO_11.data`,
   `2019_AGCO_11.data`, `2018_AGCO_17.data`, `2019_AGCO_17.data`, `2018_AGCO_18.data`, 
   `2019_AGCO_18.data`, `2018_AGCO_21.data`, `2019_AGCO_21.data`, `2019_AGCO_27.data`, 
   `2017_AGCO_05.data`, `2017_AGCO_08.data`, `2018_AGCO_10.data`, `2018_AGCO_20.data`, 
   `2018_AGCO_22.data`, `2019_AGCO_31.data`, `2019_AGCO_33.data`)

##### Generate average population raster from MCMC-sampled individual-level (snake-season) UDs #####

MCtrials <- 1000 #number of MCMC trials

blank_data_UD <- data.frame(matrix(nrow = 0, ncol = ncol(data_UD)))
colnames(blank_data_UD) <- colnames(data_UD)

## Create blank raster with extent / resolution matching kernel UDs; gets added to upon each iteration of loop
slate <- r
slate[] <- 0

for (t in 1:MCtrials){
  #t <- 5
  
  m.data_UD <- blank_data_UD #clear values of m.data_UD
  
  for (m in 1:length(mlists)){ ## randomly select one snake-season per multi-season snake and merge data
    #m <- 2 #m is an index for the a snake
    m.l <- mlists[m][[1]] #elements of m.l are the data for each season for that snake 
    m.r <- sample(m.l, 1)[[1]] #selects one random snake-season from snake m
    m.data_UD <- rbind(m.data_UD, m.r) #add each random snake season to m.data_UD
  }
  
  ## merge all single-season data with current set of randomly-selected multiple-season data
  t_data <- rbind(m.data_UD, s.data_UD)
  t_data_sp <- as(t_data, "Spatial")
  
  ## create kde UD of t_data at the pop-optimized bandwidth (pop_bw)
  k <- suppressWarnings(kernelUD(t_data_sp, h = pop_bw, grid = r.sp, kern = "bivnorm"))
  krast <- raster(k)
  
  slate <- slate+krast
  
  ## result is one population-level UD, the sum of t rasters
}

## scale raster so the values sum to 1
r_sum <- sum(values(slate))
sc.pop.rast <- slate / r_sum
plot(sc.pop.rast)
sum(values(sc.pop.rast)) ## correctly sums to 1

## extract pd.p1 values to points
pd.pop <- raster::extract(sc.pop.rast, data_UD)
data_UD$pd_pop <- pd.pop

plot(data_UD$pd_ind, data_UD$pd_pop)

## evaluate which points fall within the upper and lower third of the UD volume
values <- values(sc.pop.rast) # list of values in kde raster
values <- sort(values, decreasing = TRUE) # sort descending (highest probability, lowest volume first)
cumsum <- cumsum(values) # create vector cumulative sum of raster values

upper_index <- max(which(cumsum < 0.3))
lower_index <- min(which(cumsum > 0.3))

#upper_index <- max(which(cumsum < 0.5))
#lower_index <- min(which(cumsum > 0.5))

pop_upper <- values[upper_index]
pop_lower <- values[lower_index]

contour_lines <- rasterToContour(sc.pop.rast, levels = c(pop_upper, pop_lower))

plot(sc.pop.rast)
plot(contour_lines, add = TRUE)

data_UD$ext_pop[data_UD$pd_pop > pop_upper] <- 1
data_UD$ext_pop[data_UD$pd_pop < pop_lower] <- 0

rm(`2016_AGCO_03.data`, `2017_AGCO_03.data`, `2018_AGCO_03.data`, `2017_AGCO_06.data`, 
   `2018_AGCO_06.data`, `2019_AGCO_06.data`, `2017_AGCO_07.data`, `2018_AGCO_07.data`,
   `2017_AGCO_09.data`, `2018_AGCO_09.data`, `2019_AGCO_09.data`, `2018_AGCO_11.data`,
   `2019_AGCO_11.data`, `2018_AGCO_17.data`, `2019_AGCO_17.data`, `2018_AGCO_18.data`, 
   `2019_AGCO_18.data`, `2018_AGCO_21.data`, `2019_AGCO_21.data`, `2019_AGCO_27.data`, 
   `2017_AGCO_05.data`, `2017_AGCO_08.data`, `2018_AGCO_10.data`, `2018_AGCO_20.data`, 
   `2018_AGCO_22.data`, `2019_AGCO_31.data`, `2019_AGCO_33.data`)

rm(contour_lines, blank_data_UD, k, krast, m.data_UD, m.l, m.r, m03, m06, m07, m09, m11, m17, m18, m21, mlists,
   s.data_UD, slate, t_data, t_data_sp, cumsum, lower, m, MCtrials, pd.pop, r_sum, t, upper, values)

require(ggplot2)

ggplot(data = data_UD, aes(x = RCK_unsc, y = ext_ind)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  geom_smooth(aes(x = RCK_unsc, y = ext_pop), method = "glm", method.args = list(family = "binomial"), colour = "red") +
  ylim(0,1)

##### Model fitting #####

attach(data_UD)
par(mfrow = c(1, 2))

boxplot(CWD ~ ext_pop, main = "pop CWD")
boxplot(CWD ~ ext_ind, main = "ind CWD")

boxplot(GND ~ ext_pop, main = "pop GND")
boxplot(GND ~ ext_ind, main = "ind GND")

boxplot(WDY ~ ext_pop, main = "pop WDY")
boxplot(WDY ~ ext_ind, main = "ind WDY")

boxplot(RCK ~ ext_pop, main = "pop RCK")
boxplot(RCK ~ ext_ind, main = "ind RCK")

boxplot(CAN ~ ext_pop, main = "pop CAN")
boxplot(CAN ~ ext_ind, main = "ind CAN")

dev.off()

detach(data_UD)

##### Linear Models #####

library(MuMIn)
library(lme4)

dat <- data_UD
dat$geometry <- NULL
dat <- data.frame(dat)
ind_dat <- subset(dat, !is.na(ext_ind))
pop_dat <- subset(dat, !is.na(ext_pop))

ind_glob <- glm(ext_ind ~ CWD + GND + WDY + RCK + CAN, family = binomial(link = "logit"), data = ind_dat, na.action = "na.fail")

pop_glob <- glmer(ext_pop ~ CWD + GND + WDY + RCK + CAN + (1 | snake.season), 
                  family = "binomial", data = pop_dat, na.action = "na.fail",
                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

ind_dredge <- dredge(ind_glob) # Dredge individual global model
pop_dredge <- dredge(pop_glob) # Dredge population global model

top_ind <- subset(ind_dredge, delta < 2)
ind_avg <- model.avg(top_ind)
ind_coefs <- stats::coef(ind_avg) # Extract coefficients from average individual model
ind_confint <- stats::confint(ind_avg) # Extract CIs from average individual model
ind_summ <- data.frame(ind_confint)
ind_summ <- ind_summ[order(rownames(ind_summ)), ]
ind_summ$coefs <- ind_coefs
colnames(ind_summ) <- c("lower", "upper", "coef")

top_pop <- subset(pop_dredge, delta < 2)
pop_avg <- model.avg(top_pop)
pop_coefs <- stats::coef(pop_avg) # Extract coefficients from average population model
pop_confint <- confint(pop_avg) # Extract CIs from average population model
pop_summ <- data.frame(pop_confint)
pop_summ <- pop_summ[order(rownames(pop_summ)), ]
pop_summ$coefs <- as.numeric(pop_coefs)
colnames(pop_summ) <- c("lower", "upper", "coef")

ind_output <- summary(ind_avg)
ind_output2 <- data.frame(Estimate = ind_output$coefmat.full[,"Estimate"], SE = ind_output$coefmat.full[,"Adjusted SE"])
ind_output2$lower_CI <- ind_output2$Estimate - (1.96*ind_output2$SE)
ind_output2$upper_CI <- ind_output2$Estimate + (1.96*ind_output2$SE)
ind_output2 <- ind_output2[order(rownames(ind_output2)), , drop = FALSE]

ind_output2

ind_output
#Model-averaged coefficients:  
#  (full average) 
#              Estimate Std. Error Adjusted SE z value Pr(>|z|)    
#  (Intercept)  0.32528    0.07171     0.07181   4.530  5.9e-06 ***
#  CAN         -0.43715    0.08471     0.08483   5.153  3.0e-07 ***
#  CWD          0.30754    0.08132     0.08142   3.777 0.000159 ***
#  RCK          0.43273    0.08980     0.08992   4.812  1.5e-06 ***
#  WDY         -0.16140    0.07683     0.07694   2.098 0.035917 *  
#  GND         -0.03933    0.08002     0.08008   0.491 0.623378    

pop_output <- summary(pop_avg)
pop_output2 <- data.frame(Estimate = pop_output$coefmat.full[,"Estimate"], SE = pop_output$coefmat.full[,"Adjusted SE"])
pop_output2$lower_CI <- pop_output2$Estimate - (1.96*pop_output2$SE)
pop_output2$upper_CI <- pop_output2$Estimate + (1.96*pop_output2$SE)
pop_output2 <- pop_output2[order(rownames(pop_output2)), , drop = FALSE]

pop_output2

pop_output
#Model-averaged coefficients:  
#  (full average) 
#             Estimate Std. Error Adjusted SE z value Pr(>|z|)    
#  (Intercept)  -0.1561     0.2502      0.2505   0.623 0.533281    
#  CAN           0.4436     0.1181      0.1182   3.752 0.000175 ***
#  CWD           0.5088     0.1196      0.1198   4.247 2.17e-05 ***
#  GND          -0.7928     0.1349      0.1351   5.867  < 2e-16 ***
#  RCK           0.6475     0.1224      0.1225   5.285 1.00e-07 ***
#  WDY          -0.1055     0.1193      0.1194   0.884 0.376829   

##### Unscaled models #####

ind_glob_unsc <- glm(ext_ind ~ CWD_unsc + GND_unsc + WDY_unsc + RCK_unsc + CAN_unsc, family = binomial(link = "logit"), data = ind_dat, na.action = "na.fail")

pop_glob_unsc <- glmer(ext_pop ~ CWD_unsc + GND_unsc + WDY_unsc + RCK_unsc + CAN_unsc + (1 | snake.season), 
                  family = "binomial", data = pop_dat, na.action = "na.fail",
                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

ind_dredge_unsc <- dredge(ind_glob_unsc) # Dredge individual global model
pop_dredge_unsc <- dredge(pop_glob_unsc) # Dredge population global model

top_ind_unsc <- subset(ind_dredge_unsc, delta < 2)
ind_avg_unsc <- model.avg(top_ind_unsc)
ind_coefs_unsc <- data.frame(stats::coef(ind_avg_unsc)) # Extract coefficients from average individual model
ind_coefs_unsc <- ind_coefs_unsc[order(rownames(ind_coefs_unsc)), , drop = FALSE]
ind_confint_unsc <- confint(ind_avg_unsc) # Extract CIs from average individual model
ind_confint_unsc <- ind_confint_unsc[order(rownames(ind_confint_unsc)), ]
ind_summ_unsc <- data.frame(ind_confint_unsc)
ind_summ_unsc$coefs <- unlist(ind_coefs_unsc)
colnames(ind_summ_unsc) <- c("lower", "upper", "coef")

top_pop_unsc <- subset(pop_dredge_unsc, delta < 2)
pop_avg_unsc <- model.avg(top_pop_unsc)
pop_coefs_unsc <- data.frame(stats::coef(pop_avg_unsc)) # Extract coefficients from average population model
pop_coefs_unsc <- pop_coefs_unsc[order(rownames(pop_coefs_unsc)), , drop = FALSE]
pop_confint_unsc <- confint(pop_avg_unsc) # Extract CIs from average population model
pop_confint_unsc <- pop_confint_unsc[order(rownames(pop_confint_unsc)), ]
pop_summ_unsc <- data.frame(pop_confint_unsc)
pop_summ_unsc$coefs_unsc <- unlist(pop_coefs_unsc)
colnames(pop_summ_unsc) <- c("lower", "upper", "coef")

summary(pop_avg_unsc)

##### Probability plots #####

# Load coefficients and 95% confidence intervals
ind.coefs <- ind_summ_unsc
pop.coefs <- pop_summ_unsc

# Create data frame of variable names for figures
vars <- c("CAN", "CWD", "GND", "RCK", "WDY")
var_names <- data.frame("Canopy (%)", "Coarse Woody Debris (%)", "Ground Layer Veg. (%)", "Rock (%)", "Shrubs (%)")
colnames(var_names) <- vars

# Create list of variable names matching coefs data frames
varsE <- c("(Intercept)", "CAN", "CWD", "GND", "RCK", "WDY")

# Create blank data frame to populate with predicted values; splits into one df per variable below
df <- data.frame(matrix(nrow = 100, ncol = 4))
colnames(df) <- c("x", "mean", "lower", "upper")
df$x <- seq(1, 100, by = 1)


#plot.df <- data.frame(matrix(nrow = 6, ncol = 4))
#colnames(plot.df) <- c("var", "coef", "lower", "upper")
#plot.df$var <- varsE

#pdf <- data.frame(plot.df, row.names = 1)
#idf <- data.frame(plot.df, row.names = 1)
#pdf$coef <- pop.coefs$coef
#idf$coef <- ind.coefs$coef

####starting with population: ###################################################

# Create empty data frame for functions to calculate mean, upper, and lower predictions using beta coefficients and their upper and lower CIs
pFUNdf <- data.frame(matrix(nrow = 5, ncol = 4))
colnames(pFUNdf) <- c("variable", "estimate", "lower", "upper")
pFUNdf$variable <- vars

##Next, populate FUNdf with functions for each variable's predicted means
pCANfun <- function(x) {1/(1+exp(-(pB0 + (pBcan * x) + pTcwd + pTgnd + pTrck + pTwdy)))}
pCWDfun <- function(x) {1/(1+exp(-(pB0 + pTcan + (pBcwd * x) + pTgnd + pTrck + pTwdy)))}
pGNDfun <- function(x) {1/(1+exp(-(pB0 + pTcan + pTcwd + (pBgnd * x) + pTrck + pTwdy)))}
pRCKfun <- function(x) {1/(1+exp(-(pB0 + pTcan + pTcwd + pTgnd + (pBrck * x) + pTwdy)))}
pWDYfun <- function(x) {1/(1+exp(-(pB0 + pTcan + pTcwd + pTgnd + pTrck + (pBwdy * x))))}

pFUNdf$estimate <- c(pCANfun, pCWDfun, pGNDfun, pRCKfun, pWDYfun)
rm(pCANfun, pCWDfun, pGNDfun, pRCKfun, pWDYfun)

##Next, populate FUNdf with functions for each variable's lower CI
pCANfun.l <- function(x) {1/(1+exp(-(pB0.l + (pBcan.l * x) + pTcwd + pTgnd + pTrck + pTwdy)))}
pCWDfun.l <- function(x) {1/(1+exp(-(pB0.l + pTcan + (pBcwd.l * x) + pTgnd + pTrck + pTwdy)))}
pGNDfun.l <- function(x) {1/(1+exp(-(pB0.l + pTcan + pTcwd + (pBgnd.l * x) + pTrck + pTwdy)))}
pRCKfun.l <- function(x) {1/(1+exp(-(pB0.l + pTcan + pTcwd + pTgnd + (pBrck.l * x) + pTwdy)))}
pWDYfun.l <- function(x) {1/(1+exp(-(pB0.l + pTcan + pTcwd + pTgnd + pTrck + (pBwdy.l * x))))}

pFUNdf$lower <- c(pCANfun.l, pCWDfun.l, pGNDfun.l, pRCKfun.l, pWDYfun.l)
rm(pCANfun.l, pCWDfun.l, pGNDfun.l, pRCKfun.l, pWDYfun.l)

##Next, populate pFUNdf with functions for each variable's upper CI
pCANfun.u <- function(x) {1/(1+exp(-(pB0.u + (pBcan.u * x) + pTcwd + pTgnd + pTrck + pTwdy)))}
pCWDfun.u <- function(x) {1/(1+exp(-(pB0.u + pTcan + (pBcwd.u * x) + pTgnd + pTrck + pTwdy)))}
pGNDfun.u <- function(x) {1/(1+exp(-(pB0.u + pTcan + pTcwd + (pBgnd.u * x) + pTrck + pTwdy)))}
pRCKfun.u <- function(x) {1/(1+exp(-(pB0.u + pTcan + pTcwd + pTgnd + (pBrck.u * x) + pTwdy)))}
pWDYfun.u <- function(x) {1/(1+exp(-(pB0.u + pTcan + pTcwd + pTgnd + pTrck + (pBwdy.u * x))))}

pFUNdf$upper <- c(pCANfun.u, pCWDfun.u, pGNDfun.u, pRCKfun.u, pWDYfun.u)
rm(pCANfun.u, pCWDfun.u, pGNDfun.u, pRCKfun.u, pWDYfun.u)

##Next, need to define values of coefficients (B's) and terms (B's * x's)
pB0 <- pop.coefs["(Intercept)", "coef"]
pB0.l <- pop.coefs["(Intercept)", "lower"]
pB0.u <- pop.coefs["(Intercept)", "upper"]

pBcan <- pop.coefs["CAN_unsc", "coef"]
pBcan.l <- pop.coefs["CAN_unsc", "lower"]
pBcan.u <- pop.coefs["CAN_unsc", "upper"]

pBcwd <- pop.coefs["CWD_unsc", "coef"]
pBcwd.l <- pop.coefs["CWD_unsc", "lower"]
pBcwd.u <- pop.coefs["CWD_unsc", "upper"]

pBgnd <- pop.coefs["GND_unsc", "coef"] 
pBgnd.l <- pop.coefs["GND_unsc", "lower"]
pBgnd.u <- pop.coefs["GND_unsc", "upper"]

pBrck <- pop.coefs["RCK_unsc", "coef"]
pBrck.l <- pop.coefs["RCK_unsc", "lower"]
pBrck.u <- pop.coefs["RCK_unsc", "upper"]

pBwdy <- pop.coefs["WDY_unsc", "coef"]
pBwdy.l <- pop.coefs["WDY_unsc", "lower"]
pBwdy.u <- pop.coefs["WDY_unsc", "upper"]

pTcan <- pop.coefs["CAN_unsc", "coef"] * mean(pop_dat$CAN_unsc)
pTcwd <- pop.coefs["CWD_unsc", "coef"] * mean(pop_dat$CWD_unsc, na.rm = TRUE)
pTgnd <- pop.coefs["GND_unsc", "coef"] * mean(pop_dat$GND_unsc, na.rm = TRUE)
pTrck <- pop.coefs["RCK_unsc", "coef"] * mean(pop_dat$RCK_unsc)
pTwdy <- pop.coefs["WDY_unsc", "coef"] * mean(pop_dat$WDY_unsc)

pFUNdf <- data.frame(pFUNdf, row.names = 1)

#####Now the same for individuals

iFUNdf <- data.frame(matrix(nrow = 5, ncol = 4)) #data frame for mean, upper, and lower functions
colnames(iFUNdf) <- c("variable", "estimate", "lower", "upper")
vars <- c("CAN", "CWD", "GND", "RCK", "WDY")
iFUNdf$variable <- vars

##Next, populate FUNdf with functions for each variable's predicted means
iCANfun <- function(x) {1/(1+exp(-(iB0 + (iBcan * x) + iTcwd + iTgnd + iTrck + iTwdy)))}
iCWDfun <- function(x) {1/(1+exp(-(iB0 + iTcan + (iBcwd * x) + iTgnd + iTrck + iTwdy)))}
iGNDfun <- function(x) {1/(1+exp(-(iB0 + iTcan + iTcwd + (iBgnd * x) + iTrck + iTwdy)))}
iRCKfun <- function(x) {1/(1+exp(-(iB0 + iTcan + iTcwd + iTgnd + (iBrck * x) + iTwdy)))}
iWDYfun <- function(x) {1/(1+exp(-(iB0 + iTcan + iTcwd + iTgnd + iTrck + (iBwdy * x))))}

iFUNdf$estimate <- c(iCANfun, iCWDfun, iGNDfun, iRCKfun, iWDYfun)
rm(iCANfun, iCWDfun, iGNDfun, iRCKfun, iWDYfun)

##Next, populate FUNdf with functions for each variable's lower CI
iCANfun.l <- function(x) {1/(1+exp(-(iB0.l + (iBcan.l * x) + iTcwd + iTgnd + iTrck + iTwdy)))}
iCWDfun.l <- function(x) {1/(1+exp(-(iB0.l + iTcan + (iBcwd.l * x) + iTgnd + iTrck + iTwdy)))}
iGNDfun.l <- function(x) {1/(1+exp(-(iB0.l + iTcan + iTcwd + (iBgnd.l * x) + iTrck + iTwdy)))}
iRCKfun.l <- function(x) {1/(1+exp(-(iB0.l + iTcan + iTcwd + iTgnd + (iBrck.l * x) + iTwdy)))}
iWDYfun.l <- function(x) {1/(1+exp(-(iB0.l + iTcan + iTcwd + iTgnd + iTrck + (iBwdy.l * x))))}

iFUNdf$lower <- c(iCANfun.l, iCWDfun.l, iGNDfun.l, iRCKfun.l, iWDYfun.l)
rm(iCANfun.l, iCWDfun.l, iGNDfun.l, iRCKfun.l, iWDYfun.l)

##Next, populate iFUNdf with functions for each variable's upper CI
iCANfun.u <- function(x) {1/(1+exp(-(iB0.u + (iBcan.u * x) + iTcwd + iTgnd + iTrck + iTwdy)))}
iCWDfun.u <- function(x) {1/(1+exp(-(iB0.u + iTcan + (iBcwd.u * x) + iTgnd + iTrck + iTwdy)))}
iGNDfun.u <- function(x) {1/(1+exp(-(iB0.u + iTcan + iTcwd + (iBgnd.u * x) + iTrck + iTwdy)))}
iRCKfun.u <- function(x) {1/(1+exp(-(iB0.u + iTcan + iTcwd + iTgnd + (iBrck.u * x) + iTwdy)))}
iWDYfun.u <- function(x) {1/(1+exp(-(iB0.u + iTcan + iTcwd + iTgnd + iTrck + (iBwdy.u * x))))}

iFUNdf$upper <- c(iCANfun.u, iCWDfun.u, iGNDfun.u, iRCKfun.u, iWDYfun.u)
rm(iCANfun.u, iCWDfun.u, iGNDfun.u, iRCKfun.u, iWDYfun.u)

##Next, need to define values of coefficients (B's) and terms (B's * x's)
iB0 <- ind.coefs["(Intercept)", "coef"]
iB0.l <- ind.coefs["(Intercept)", "lower"]
iB0.u <- ind.coefs["(Intercept)", "upper"]

iBcan <- ind.coefs["CAN_unsc", "coef"]
iBcan.l <- ind.coefs["CAN_unsc", "lower"]
iBcan.u <- ind.coefs["CAN_unsc", "upper"]

iBcwd <- ind.coefs["CWD_unsc", "coef"]
iBcwd.l <- ind.coefs["CWD_unsc", "lower"]
iBcwd.u <- ind.coefs["CWD_unsc", "upper"]

iBgnd <- ind.coefs["GND_unsc", "coef"] 
iBgnd.l <- ind.coefs["GND_unsc", "lower"]
iBgnd.u <- ind.coefs["GND_unsc", "upper"]

iBrck <- ind.coefs["RCK_unsc", "coef"]
iBrck.l <- ind.coefs["RCK_unsc", "lower"]
iBrck.u <- ind.coefs["RCK_unsc", "upper"]

iBwdy <- ind.coefs["WDY_unsc", "coef"]
iBwdy.l <- ind.coefs["WDY_unsc", "lower"]
iBwdy.u <- ind.coefs["WDY_unsc", "upper"]

iTcan <- ind.coefs["CAN_unsc", "coef"] * mean(ind_dat$CAN_unsc, na.rm = TRUE)
iTcwd <- ind.coefs["CWD_unsc", "coef"] * mean(ind_dat$CWD_unsc, na.rm = TRUE)
iTgnd <- ind.coefs["GND_unsc", "coef"] * mean(ind_dat$GND_unsc, na.rm = TRUE)
iTrck <- ind.coefs["RCK_unsc", "coef"] * mean(ind_dat$RCK_unsc, na.rm = TRUE)
iTwdy <- ind.coefs["WDY_unsc", "coef"] * mean(ind_dat$WDY_unsc)

iFUNdf <- data.frame(iFUNdf, row.names = 1)

#####

df <- data.frame(matrix(nrow = 100, ncol = 4))
colnames(df) <- c("x", "mean", "upper", "lower")
df$x <- seq(1, 100, by = 1)

vars
#var_df

for (v in vars) {
  #v <- vars[1] #vars is now, e.g., "CAN"
  p.df <- df
  i.df <- df
  
  p.m.func <- eval(parse(text=pFUNdf[v,]$estimate))
  p.l.func <- eval(parse(text=pFUNdf[v,]$lower))
  p.u.func <- eval(parse(text=pFUNdf[v,]$upper))
  p.df$mean <- sapply(p.df$x, p.m.func)
  p.df$lower <- sapply(p.df$x, p.l.func)
  p.df$upper <- sapply(p.df$x, p.u.func)
  
  i.m.func <- eval(parse(text=iFUNdf[v,]$estimate))
  i.l.func <- eval(parse(text=iFUNdf[v,]$lower))
  i.u.func <- eval(parse(text=iFUNdf[v,]$upper))
  i.df$mean <- sapply(i.df$x, i.m.func)
  i.df$lower <- sapply(i.df$x, i.l.func)
  i.df$upper <- sapply(i.df$x, i.u.func)
  
  assign(paste(v, "_plot", sep = ""),  ggplot(data = p.df, aes(x = x)) +
           xlim(0, 100) +
           ylim(0,1) +
           ylab("Probability") +
           xlab(var_names[v]) +
           theme_bw() +
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #grid pattern
           theme(axis.text = element_text(size = 30)) + #axis text size
           theme(axis.title = element_text(size = 40)) + #axis title size
           theme(axis.text.x= element_text(colour = "#212121")) + #axis number color
           theme(axis.text.y= element_text(colour = "#212121")) + #axis number color
           theme(panel.border = element_rect(colour = "black", fill=NA, size=2)) +
           theme(plot.title = element_text(hjust = 0.5)) + #plot title
           theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) +
           ##population probability curve
           geom_ribbon(data = p.df, aes(x = x, ymin = lower, ymax = upper), fill = "dodgerblue4", alpha = 0.5) + ## ribbon
           geom_line(data = p.df, aes(x = x, y = mean), size = 1.2, colour = "blue") + ##prediction line
           ##individual probability curve
           geom_ribbon(data = i.df, aes(x = x, ymin = lower, ymax = upper), fill = "firebrick", alpha = 0.5) + ## ribbon
           geom_line(data = i.df, aes(x = x, y = mean), size = 1.2, colour = "red3") ##prediction line
  )
  
}

grid.arrange(CAN_plot, CWD_plot, GND_plot, RCK_plot, ncol=2)
dev.print(png, width = 1500, height = 1000, "prob_plots.PNG")

rm(CAN_plot, CWD_plot, dat, data_UDsp, df, GND_plot, i.df, iFUNdf, ind_dat, ind_avg_unsc,
   ind_coefs_unsc, ind_confint, ind_confint_unsc, ind_dredge, ind_dredge_unsc, ind_glob, ind_glob_unsc,
   ind_summ, ind_summ_unsc, ind_coefs, p.df, pFUNdf, pop_avg_unsc, pop_coefs_unsc,
   pop_confint, pop_confint_unsc, pop_dat, pop_dredge, pop_dredge_unsc, pop_glob, pop_glob_unsc, 
   pop.coefs, r, r.sp, RCK_plot, sc.pop.rast, ssndf, top_ind, top_ind_unsc, top_pop, top_pop_unsc,
   var_names, WDY_plot, pop_summ, pop_summ_unsc)

rm(iB0, iB0.l, iB0.u, iBcan, iBcan.l, iBcan.u, iBcwd, iBcwd.l, iBcwd.u, iBgnd, iBgnd.l, iBgnd.u,
   iBrck, iBrck.l, iBrck.u, iBwdy, iBwdy.l, iBwdy.u, ind.coefs, iTcan, iTcwd, iTgnd, iTrck, iTwdy,
   lower_index, pB0, pB0.l, pB0.u, pBcan, pBcan.l, pBcan.u, pBcwd, pBcwd.l, pBcwd.u, pBgnd, pBgnd.l,
   pBgnd.u, pBrck, pBrck.l, pBrck.u, pBwdy, pBwdy.l, pBwdy.u, pop_coefs, pop_lower, pop_upper,
   pTcan, pTcwd, pTgnd, pTrck, pTwdy, upper_index, v, vars, varsE, i.l.func, i.m.func, i.u.func,
   p.l.func, p.m.func, p.u.func)

##### Coefficient equality tests #####

ind_summ <- summary(ind_avg)
ind_coefs <- data.frame(ind_summ$coefmat.full[,1:2])
ind_coefs <- ind_coefs[order(rownames(ind_coefs)), ]

pop_summ <- summary(pop_avg)
pop_coefs <- data.frame(pop_summ$coefmat.full[,1:2])
pop_coefs <- pop_coefs[order(rownames(pop_coefs)), ]

vars <- rownames(pop_coefs)

for (v in rownames(ind_coefs)) {
  #v <- rownames(ind_coefs)[2]
  beta1 <- ind_coefs[v, "Estimate"]
  beta2 <- pop_coefs[v, "Estimate"]
  
  SEbeta1 <- ind_coefs[v, "Std..Error"]
  SEbeta2 <- pop_coefs[v, "Std..Error"]
  
  Wald.stat <- ((beta1 - beta2)^2) / (SEbeta1^2 + SEbeta2^2)
  
  print(paste(v, Wald.stat))
  
  p_value <- 1 - pchisq(Wald.stat, df = 1)
  
  name <- paste(v, "_pval", sep = "")
  
  assign(name, p_value)
  
}

Wald.sig <- t(data.frame(CWD_pval, GND_pval, WDY_pval, RCK_pval, CAN_pval))
Wald.sig <- data.frame(Wald.sig[order(rownames(Wald.sig)), ])

Wald.sig

rm(`(Intercept)_pval`, beta1, beta2, CAN_pval, CWD_pval, GND_pval, name, p_value, RCK_pval,
   SEbeta1, SEbeta2, v, vars, Wald.stat, WDY_pval)

##### Bootstrapped covariate means at high and low use areas #####

ind.high <- data_UD %>% dplyr::filter(ext_ind == 1)
ind.low <- data_UD %>% dplyr::filter(ext_ind == 0)

pop.high <- data_UD %>% dplyr::filter(ext_pop == 1)
pop.low <- data_UD %>% dplyr::filter(ext_pop == 0)

##### HIGH USE AREAS bootstrapped covariate means #####

library(boot)

ind.high <- data_UD %>% dplyr::filter(ext_ind == 1)
ind.low <- data_UD %>% dplyr::filter(ext_ind == 0)

pop.high <- data_UD %>% dplyr::filter(ext_pop == 1)
pop.low <- data_UD %>% dplyr::filter(ext_pop == 0)

vars <- c("CWD_unsc", "GND_unsc", "WDY_unsc", "RCK_unsc", "CAN_unsc")

for(v in vars) {
  #v <- vars[1]
  cov <- ind.high[[v]]
  hist(cov, breaks = 10, main = v)
  avg <- function(cov,i) mean(cov[i], na.rm = TRUE)
  boot <- boot(cov, avg, R = 1000)
  CI <- 1.96*sd(boot$t)
  print(paste("Ind High Use: ", v, ", mean = ", boot$t0, ", 95% CI = ", CI, sep = ""))
}

for(v in vars) {
  #v <- vars[1]
  cov <- ind.low[[v]]
  hist(cov, breaks = 10, main = v)
  avg <- function(cov,i) mean(cov[i], na.rm = TRUE)
  boot <- boot(cov, avg, R = 1000)
  CI <- 1.96*sd(boot$t)
  print(paste("Ind Low Use: ", v, ", mean = ", boot$t0, ", 95% CI = ", CI, sep = ""))
}

for(v in vars) {
  #v <- vars[1]
  cov <- pop.high[[v]]
  hist(cov, breaks = 10, main = v)
  avg <- function(cov,i) mean(cov[i], na.rm = TRUE)
  boot <- boot(cov, avg, R = 1000)
  CI <- 1.96*sd(boot$t)
  print(paste("Pop High Use: ", v, ", mean = ", boot$t0, ", 95% CI = ", CI, sep = ""))
}

for(v in vars) {
  #v <- vars[1]
  cov <- pop.low[[v]]
  hist(cov, breaks = 10, main = v)
  avg <- function(cov,i) mean(cov[i], na.rm = TRUE)
  boot <- boot(cov, avg, R = 1000)
  CI <- 1.96*sd(boot$t)
  print(paste("Pop Low Use: ", v, ", mean = ", boot$t0, ", 95% CI = ", CI, sep = ""))
}

##### PCA to compare high use areas for inds, pop, and repr classes #####

require(dplyr)

## Need to create three new columns containing "group" if observation was high-use in one or more groups

data <- data.frame(data_UD)

data$PC1 <- NA
data$PC2 <- NA
data$PC3 <- NA
PCA.df <- data %>% dplyr::select(year_repr, CWD, GND, WDY, RCK, CAN, ext_ind, ext_pop, PC1, PC2, PC3)

head(PCA.df)

## create columns
PCA.df$rep_group <- "" ##contains char c("GVD", "NGR", "MAL") if observation is of ind high use
PCA.df$ind_group <- "" ##contains char "ind" if observation is of ind high use
PCA.df$pop_group <- "" ## contains char "pop" if observation is of high pop use

for(r in 1:nrow(PCA.df)){
  #r <- 8
  
  ## Replaces rep_group with reproductive class IF obs is from ind high-use area
  cond1 <- PCA.df[r,"ext_ind"] == 1
  repl1 <- PCA.df[r,"year_repr"]
  PCA.df[r,"rep_group"] <- ifelse(cond1, repl1, NA)
  
  ## Replaces ind_group with "IND" if observation is from ind high-use area
  repl2 <- "IND"
  PCA.df[r, "ind_group"] <- ifelse(cond1, repl2, NA)
  
  ## Replaces pop_group with "POP" if observation is from pop high-use area
  cond3 <- PCA.df[r,"ext_pop"] == 1
  repl3 <- "POP"
  PCA.df[r, "pop_group"] <- ifelse(cond3, repl3, NA)
}

## Remove observations that did not appear in any high-use areas

PCA.df2 <- PCA.df %>% filter(!is.na(rep_group) | !is.na(ind_group) | !is.na(pop_group))

PCA <- prcomp(PCA.df2[,c("CAN", "RCK", "CWD", "GND", "WDY")], scale. = TRUE)
summary(PCA)
PCA

PCA.df2$PC1 <- PCA$x[,1]
PCA.df2$PC2 <- PCA$x[,2]
PCA.df2$PC3 <- PCA$x[,3]

hist(PCA.df2$PC1)
hist(PCA.df2$PC2)

## Need to create a new data frame with a single "group" column; observations that
## appear in multiple groups will be duplicated

rep <- PCA.df2 %>% filter(!is.na(rep_group))
ind <- PCA.df2 %>% filter(!is.na(ind_group))
pop <- PCA.df2 %>% filter(!is.na(pop_group))

rep$group <- rep$rep_group
ind$group <- "IND"
pop$group <- "POP"

PCA.df3 <- rbind(rep, ind, pop)

aov1 <- aov(PC1 ~ group, data = PCA.df3)
tukey <- TukeyHSD(aov1)

p1 <- tukey$`group`[, "p adj"]
p1.names <- names(p1)

adjusted_p_values1 <- p.adjust(p1, method = "bonferroni")

## summary of PC1 results
pc1df <- data.frame(PC = rep("PC1", length(p1)), comp = p1.names, b.p.val = adjusted_p_values1)

aov2 <- aov(PC2 ~ group, data = PCA.df3)
tukey2 <- TukeyHSD(aov2)

p2 <- tukey2$`group`[, "p adj"]
p2.names <- names(p2)

adjusted_p_values2 <- p.adjust(p2, method = "bonferroni")

## summary of PC2 results
pc2df <- data.frame(PC = rep("PC2", length(p1)), comp = p2.names, b.p.val = adjusted_p_values2)

PC1andPC2 <- rbind(pc1df, pc2df)
PC1andPC2

write.table(PC1andPC2, "2024-08-13_summary_of_PCA_ANOVA_Tukey_Tests.xlsx")

#### PC1 and PC2 stats for plotting #####

PC.list <- c("PC1", "PC2", "PC3")

PCA.df3$group <- factor(PCA.df3$group)
boxplot(PCA.df3$CAN ~ PCA.df3$group)

for(g in levels(PCA.df3$group)) {
  #g <- levels(PCA.df3$group)[5]
  gdf <- PCA.df3 %>% filter(group == g)
  
  for (p in PC.list) {
    #p <- PC.list[1]
    gXp.mean <-  mean(gdf[[p]])
    gXp.name1 <- paste(g, ".", p, ".mean", sep = "")
    assign(gXp.name1, gXp.mean)
    
    gXp.se <- sd(gdf[[p]])/sqrt(nrow(gdf))
    gXp.name2 <- paste(g, ".", p, ".se", sep = "")
    assign(gXp.name2, gXp.se)
    
    gXp.CI <- 1.96*gXp.se
    gXp.name3 <- paste(g, ".", p, ".CI", sep = "")
    assign(gXp.name3, gXp.CI)
  }
  
}

PC1.full <- data.frame(group = c("POP", "IND", "GVD", "NGR", "MAL"), 
                       mean = c(POP.PC1.mean, IND.PC1.mean, GVD.PC1.mean, NGR.PC1.mean, MAL.PC1.mean),
                       se = c(POP.PC1.se, IND.PC1.se, GVD.PC1.se, NGR.PC1.se, MAL.PC1.se),
                       CI = c(POP.PC1.CI, IND.PC1.CI, GVD.PC1.CI, NGR.PC1.CI, MAL.PC1.CI))

PC2.full <- data.frame(group = c("POP", "IND", "GVD", "NGR", "MAL"), 
                       mean = c(POP.PC2.mean, IND.PC2.mean, GVD.PC2.mean, NGR.PC2.mean, MAL.PC2.mean),
                       se = c(POP.PC2.se, IND.PC2.se, GVD.PC2.se, NGR.PC2.se, MAL.PC2.se),
                       CI = c(POP.PC2.CI, IND.PC2.CI, GVD.PC2.CI, NGR.PC2.CI, MAL.PC2.CI))

PC3.full <- data.frame(group = c("POP", "IND", "GVD", "NGR", "MAL"), 
                       mean = c(POP.PC3.mean, IND.PC3.mean, GVD.PC3.mean, NGR.PC3.mean, MAL.PC3.mean),
                       se = c(POP.PC3.se, IND.PC3.se, GVD.PC3.se, NGR.PC3.se, MAL.PC3.se),
                       CI = c(POP.PC3.CI, IND.PC3.CI, GVD.PC3.CI, NGR.PC3.CI, MAL.PC3.CI))

PC1.full$group <- factor(PC1.full$group, levels = c("MAL", "NGR", "GVD", "IND", "POP"))
PC2.full$group <- factor(PC2.full$group, levels = c("MAL", "NGR", "GVD", "IND", "POP"))

plot1 <- ggplot(PC1.full, aes(x=group, y = mean)) +
  scale_colour_manual(values = c("violetred", "mediumaquamarine", "darkgoldenrod1", "red3","royalblue4")) +
  coord_flip() +
  geom_hline(yintercept = c(-1, 0, 1), color = "gray", size = 1) +
  geom_errorbar(aes(ymin = mean - CI, ymax = mean + CI, colour = factor(group)), width = 0, linewidth = 4, alpha = 0.5) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, colour = factor(group)), width = 0, linewidth = 6, alpha = 1) +
  geom_point(aes(colour = factor(group)), size = 16) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "none",
        axis.text = element_text(size = 30),
        axis.text.x= element_text(colour = "#212121"),
        axis.text.y= element_text(colour = "#212121"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) +
  ylab(NULL) +
  xlab(NULL) +
  scale_y_continuous(position = "right", limits = c(-1.1, 1.1)) +
  scale_x_discrete(position = "top",
                   labels = c("POP" = "Population", 
                              "IND" = "Individual", 
                              "GVD" = "Gravid",
                              "NGR" = "Nongravid",
                              "MAL" = "Male"))

plot2 <- ggplot(PC2.full, aes(x=factor(group, levels = c("MAL", "NGR", "GVD", "IND", "POP")), y = mean)) +
  scale_colour_manual(values = c("violetred", "mediumaquamarine", "darkgoldenrod1", "red3","royalblue4")) +
  coord_flip() +
  geom_hline(yintercept = c(-1, 0, 1), color = "gray", size = 1) +
  geom_errorbar(aes(ymin = mean - CI, ymax = mean + CI, colour = factor(group)), width = 0, linewidth = 4, alpha = 0.5) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, colour = factor(group)), width = 0, linewidth = 6, alpha = 1) +
  geom_point(aes(colour = factor(group)), size = 16) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "none",
        axis.text = element_text(size = 30),
        axis.text.x= element_text(colour = "#212121"),
        axis.text.y= element_text(colour = "#212121"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) +
  ylab(NULL) +
  xlab(NULL) +
  scale_y_continuous(position = "right", limits = c(-1.1, 1.1)) +
  scale_x_discrete(position = "top",
                   labels = c("POP" = "Population", 
                              "IND" = "Individual", 
                              "GVD" = "Gravid",
                              "NGR" = "Nongravid",
                              "MAL" = "Male"))

library(gridExtra)
grid.arrange(plot1, plot2, ncol = 1)
dev.print(png, width = 1000, height = 800, "PCA_plots.PNG")

rm(aov1, aov2, gdf, ind, ind_avg, ind_summ, ind_coefs, PC1.full, PC1andPC2, pc1df, PC2.full, pc2df, PC3.full,
   PCA, PCA.df, PCA.df2, PCA.df3, pop, pop_avg, pop_coefs, pop_summ, rep, tukey, tukey2, Wald.sig, 
   adjusted_p_values1, adjusted_p_values2, cond1, cond3, g, GVD.PC1.CI, GVD.PC1.mean, GVD.PC1.se,
   GVD.PC2.CI, GVD.PC2.mean, GVD.PC2.se, GVD.PC3.CI, GVD.PC3.mean, GVD.PC3.se, gXp.CI, gXp.mean, gXp.name1,
   gXp.name2, gXp.name3, gXp.se, IND.PC1.CI, IND.PC1.mean, IND.PC1.se, IND.PC2.CI, IND.PC2.mean,
   IND.PC2.se, IND.PC3.CI, IND.PC3.mean, IND.PC3.se, MAL.PC1.mean, MAL.PC1.CI, MAL.PC1.se, MAL.PC2.CI,
   MAL.PC2.mean, MAL.PC2.se, MAL.PC3.CI, MAL.PC3.mean, MAL.PC3.se, NGR.PC2.mean, NGR.PC2.se, NGR.PC3.CI, 
   NGR.PC3.mean, NGR.PC3.se, NGR.PC1.mean, NGR.PC1.CI, NGR.PC1.se, NGR.PC2.CI, p, p1, p1.names, p2, 
   p2.names, PC.list, POP.PC1.CI, POP.PC1.mean, POP.PC1.se, POP.PC2.CI, POP.PC2.mean, POP.PC2.se, 
   POP.PC3.CI, POP.PC3.mean, POP.PC3.se, r, repl1, repl2, repl3, plot1, plot2)

##### Supplemental Analyses #####

## spaMM models (2024-08-11)

library(spaMM)

datsp <- as.data.frame(data_UD)
coordinates(datsp) <- c("long", "lat")
datsp <- st_as_sf(datsp, CRS("+init=epsg:26918"))

coords <- st_coordinates(datsp) #retrieve coordinates
colnames(coords) <- c("long", "lat") #rename coordinate columns
datsp <- cbind(datsp, coords) #merge

pop.glob <- fitme(ext_pop ~ CWD + GND + WDY + RCK + CAN + Matern(1|long+lat), data=datsp, family = "binomial")
summary(pop.glob, details = list(p_value = TRUE))

##### semivariograms for supplement #####

library(geoR)

data2 <- data
data2$geometry <- NULL

varnames <- c("CAN", "CWD", "GND", "RCK", "WDY")

for (v in varnames) {
  #v <- varnames[1]
  vario <- variog(coords = matrix(c(data$long, data$lat), ncol = 2), data = data2[,v], max.dist = 300, option = "smooth")
  plot.name <- paste(v, "_vario", sep = "")
  assign(plot.name, vario)
}

par(mfrow = c(2, 3), 
    cex.main = 4,  # Main title font size
    cex.lab = 3,   # Axis labels font size
    cex.axis = 2,    # Axis tick labels font size
    mar = c(5, 5, 5, 5))

plot(CAN_vario, main = "Canopy Cover (%)")
plot(CWD_vario, main = "Woody Debris (%)")
plot(GND_vario, main = "Ground Layer (%)")
plot(RCK_vario, main = "Rock (%)") 
plot(WDY_vario, main = "Woody Shrubs (%)")

dev.print(png, width = 1000, height = 700, "semivariograms.PNG")


##### spaMM models for supplement #####

library(spaMM)
library(RSpectra)
library(ROI.plugin.glpk)

dat <- data_UD

datsp <- as.data.frame(dat)
coordinates(datsp) <- c("long", "lat")
datsp <- st_as_sf(datsp, CRS("+init=epsg:26918"))

coords <- st_coordinates(datsp) #retrieve coordinates
colnames(coords) <- c("long", "lat") #rename coordinate columns
datsp <- cbind(datsp, coords) #merge

ind.glob <- fitme(ext_ind ~ CWD + GND + WDY + RCK + CAN + Matern(1|long+lat), data=datsp, family = "binomial")
summary(ind.glob, details = list(p_value = TRUE))

#------------ Fixed effects (beta) ------------
#            Estimate Cond. SE t-value   p-value
#(Intercept)  -1.8551   0.3588 -5.1700 2.911e-07
#CWD           0.1458   0.1724  0.8457 3.980e-01
#GND           0.3108   0.1743  1.7828 7.497e-02
#WDY          -0.2066   0.1457 -1.4177 1.566e-01
#RCK           0.3219   0.2014  1.5985 1.103e-01
#CAN          -0.6220   0.1646 -3.7793 1.681e-04 ***

pop.glob <- fitme(ext_pop ~ CWD + GND + WDY + RCK + CAN + Matern(1|long+lat), data=datsp, family = "binomial")
summary(pop.glob, details = list(p_value = TRUE))

#------------ Fixed effects (beta) ------------
#            Estimate Cond. SE t-value  p-value
#(Intercept) -47.6307   17.746 -2.6840 0.007413
#CWD           2.4567    2.435  1.0088 0.313347
#GND          -1.0469    1.678 -0.6238 0.532906
#WDY           0.2703    1.441  0.1876 0.851265
#RCK           0.3766    3.162  0.1191 0.905237
#CAN          -0.8567    1.898 -0.4514 0.651838
#--------------- Random effects ---------------
#  Family: gaussian( link = identity ) 
#--- Correlation parameters:
#  1.nu      1.rho 
#6.39022503 0.09520131 
#--- Variance parameters ('lambda'):
#  lambda = var(u) for u ~ Gaussian; 
#long + lat  :  3016  
# of obs: 873; # of groups: long + lat, 871 
#------------- Likelihood values  -------------
#  logLik
#logL       (p_v(h)): -65.14082

##### Code for pop-level volume contours (for figures) #####

kde <- kernelUD(data_geom, h = 31, grid = r.sp, kern = "bivnorm") ##kernel UD using optimized bandwidth (equivalent to 'final' k of for loop)

pop.99 <- getverticeshr(kde, percent = 99) ## min "low" bin
pop.66 <- getverticeshr(kde, percent = 66.67) ## max "low" bin
pop.low <- gDifference(pop.99, pop.66) ## pop "low" area (99 to 66% volume contours)

pop.high <- getverticeshr(kde, percent = 33.33) ## min "high" bin

plot(pop.low, col="steelblue")
plot(pop.high, col = "red", add = TRUE)

pop.low.sf <- st_as_sf(pop.low)
pop.high.sf <- st_as_sf(pop.high)

name.lower <- paste("pop_UD_lower3rd.shp", sep = "")
name.upper <- paste("pop_UD_upper3rd.shp", sep = "")

st_write(pop.low.sf, name.lower, driver = "ESRI Shapefile", append = FALSE)
st_write(pop.high.sf, name.upper, driver = "ESRI Shapefile", append = FALSE)

ssns <- levels(data_UD$snake.season)

library(rgeos)

for (s in ssns) {
  #s <-  ssns[1]
  ssn <- data_UD %>% dplyr::filter(snake.season == s) # subset data_UD; only snake.season s observations
  index <- which(ssns == s) #which element of ssns is s? Use this as an index for bw.
  bw <- as.numeric(bw_df[index, 2]) # Index the bw from bw_df
  ssn.geom <- as(ssn$geometry, "Spatial") #convert ssn to SpatialPoints for kde
  ssn.kde <- kernelUD(ssn.geom, h = bw, grid = r.sp, kern = "bivnorm") ##kernel UD using ind-optimized bandwidth
  
  ind.99 <- getverticeshr(ssn.kde, percent = 99) ## min "low" bin
  ind.66 <- getverticeshr(ssn.kde, percent = 66.67) ## max "low" bin
  
  ind.low <- gDifference(ind.99, ind.66) ## ind "low" area (99 to 66% volume contours)
  ind.high <- getverticeshr(ssn.kde, percent = 33.33) ## min "high" bin
  
  plot(ind.low, col="steelblue")
  plot(ind.high, col = "red", add = TRUE)
  
  ind.low.sf <- st_as_sf(ind.low)
  ind.high.sf <- st_as_sf(ind.high)
  
  name.lower <- paste(s, "_UD_lower3rd.shp", sep = "")
  name.upper <- paste(s, "_UD_upper3rd.shp", sep = "")
  
  st_write(ind.low.sf, name.lower, driver = "ESRI Shapefile", append = FALSE)
  st_write(ind.high.sf, name.upper, driver = "ESRI Shapefile", append = FALSE)

}

