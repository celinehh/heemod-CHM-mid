# MID-STAGE PATIENT 
# 45 years old, starts in hindered. need to shift all eqns by 25

rm(list = ls())      # clear memory (removes all the variables from the workspace)

#### 01 Load packages ####
library(ggplot2)
library(dampack)
library(dplyr)
library(scales)
library(ellipse)

#### 02 Load Functions ####
source("/Users/celine/Dropbox/DPhil/Winter 2018:2019/CEA and decision modeling in R/Participant Material/Functions.R")

#### 03 Input Model Parameters ####
## Strategy names
v.names.str <- c("No Treatment", "Treatment")  
## Number of strategies
n.str <- length(v.names.str)

## Markov model parameters
age     <- 45                                #NOT AGE # age at baseline 
max.age <- 100                                 # maximum age of follow up
n.t  <- max.age - age                         # time horizon, number of cycles
v.n  <- c("Get_GT", "Healthy", "GT_Healthy", "Hindered1", "GT_Hindered1", "Hindered2","FBlind1","FBlind2","Blind1","Blind2", "Death")               # the 4 states of the model: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
# the 4 states of the model: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
# NOT STATE NUMBERS        # the 4 states of the model: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
n.s <- length(v.n)                            # number of health states 

## Transition probabilities (per cycle) and hazard ratios
# Update!
###
# Read age-specific mortality rates from csv file
lt.usa.2005 <- read.csv("/Users/celine/Dropbox/DPhil/Winter 2018:2019/CEA and decision modeling in R/Participant Material/Day 3/HMD_USA_Mx_2015.csv")
v.r.HD <- lt.usa.2005 %>% 
  filter(Age >= age & Age <= (max.age-1)) %>%
  select(Total) %>%
  as.matrix()
v.r.HD
lt.usa.2005
as.matrix( #turns a data table into a matrix, format is [datatable],
  select( #selects the variables we care about only, deletes the rest!
    filter( #finds rows/cases where a certain statement is true
      lt.usa.2005, Age >= age & Age <= (max.age-1)),Total))  #data table

#create a matrix
a=45:99
one <- 1^a
one
as.matrix(one)
#i will need to get all the values over time of my variables, then convert into a matrix, then use this in the trans probs matrix
###
p.GT.success<- 0.93
p.GT.fail <- 0.07
pRA1_CHM=0.5
pRA2_CHM=0.0354
pRA3_CHM=0.0354
age_range=45:99
eqn1=(-.00000010374)*age_range ^ 4 + (.000005428)*age_range ^ 3 + (.0002224)*age_range ^ 2 - (.0009174)*age_range + 0.04618
CHMVAloss1=ifelse(age_range<=59,eqn1,ifelse(age_range<=75,.7,.85)) #CHM associated VA loss NOT MATRIX
CHMVAloss<-as.matrix(CHMVAloss1)
nonCHMVA1=ifelse(age_range<=29,.00003318,ifelse(age_range<=39,.00012,ifelse(age_range<=49,.0001499,ifelse(age_range<=59,.00008081,ifelse(age_range<=69,.00116,ifelse(age_range<=79,.002453,.01153)))))) #non-CHM associated VA loss
nonCHMVA<-as.matrix(nonCHMVA1)
nonCHMVA
pDead <- lt.usa.2005 %>% 
  filter(Age >= age & Age <= (max.age-1)) %>%
  select(Total) %>%
  as.matrix()
nonCHMblind1=ifelse(age_range<=29,.000043,ifelse(age_range<=39,.0000926,ifelse(age_range<=49,.0001186,ifelse(age_range<=59,.00004818,ifelse(age_range<=69,.0001519,ifelse(age_range<=79,.0003448,.002651)))))) #non-CHM blindness
nonCHMblind<-as.matrix(nonCHMblind1)
#trans probs, Healthy
p.Healthy_Healthy <- (one-pDead)*(one-nonCHMblind)*(one-(CHMVAloss+nonCHMVA))*(one-pRA1_CHM)
p.Healthy_Hindered1 <- (one-pDead)*(one-nonCHMblind)*(CHMVAloss+nonCHMVA)
p.Healthy_Hindered2 <- (one-pDead)*(one-nonCHMblind)*(one-(CHMVAloss+nonCHMVA))*(pRA1_CHM)
p.Healthy_Blind2 <- (one-pDead)*(nonCHMblind)
#GT trans probs, GT-Healthy
p.GT.GTHealthy_GTHealthy <- (one-pDead)*(one-nonCHMblind)*(one-nonCHMVA)
p.GT.GTHealthy_GTHindered1 <- (one-pDead)*(one-nonCHMblind)*(nonCHMVA)
p.GT.GTHealthy_Blind2 <- (one-pDead)*(nonCHMblind)

#trans probs, Hindered1
p.Hindered1_Hindered1 <- (one-pDead)*(one-nonCHMblind)*(one-pRA1_CHM)
p.Hindered1_FBlind1 <- (one-pDead)*(one-nonCHMblind)*(pRA1_CHM)
p.Hindered1_Blind2 <- (one-pDead)*(nonCHMblind)
#GT trans probs, GT-Hindered1
p.GT.GTHindered1_GTHindered1 <- (one-pDead)*(one-nonCHMblind)
p.GT.GTHindered1_Blind2 <- (one-pDead)*(nonCHMblind)
#trans probs, Hindered2
p.Hindered2_Hindered2 <- (one-pDead)*(one-nonCHMblind)*(one-(CHMVAloss+nonCHMVA))*(one-pRA2_CHM)
p.Hindered2_FBlind1 <- (one-pDead)*(one-nonCHMblind)*(CHMVAloss+nonCHMVA)
p.Hindered2_FBlind2 <- (one-pDead)*(one-nonCHMblind)*(one-(CHMVAloss+nonCHMVA))*(pRA2_CHM)
p.Hindered2_Blind2 <- (one-pDead)*(nonCHMblind)
#trans probs, FBlind1
p.FBlind1_FBlind1 <- (one-pDead)*(one-nonCHMblind)*(one-pRA2_CHM)
p.FBlind1_Blind1 <- (one-pDead)*(one-nonCHMblind)*(pRA2_CHM)
p.FBlind1_Blind2 <- (one-pDead)*(nonCHMblind)
#trans probs, FBlind2
p.FBlind2_FBlind2 <- (one-pDead)*(one-nonCHMblind)*(one-(CHMVAloss+nonCHMVA))*(one-pRA3_CHM)
p.FBlind2_Blind1 <- (one-pDead)*(one-nonCHMblind)*(CHMVAloss+nonCHMVA)
p.FBlind2_Blind2 <- (one-pDead)*(nonCHMblind)+ (one-pDead)*(one-nonCHMblind)*(one-(CHMVAloss+nonCHMVA))* (pRA3_CHM)
#trans probs, Blind1
p.Blind1_Blind1 <- (one-pDead)*(one-pRA3_CHM)
p.Blind1_Blind2 <- (one-pDead)*(pRA3_CHM)
#trans probs, Blind2
p.Blind2_Blind2 <- one-pDead
#trans probs, dead
p.Dead_Dead <- 1

## Cost and utility inputs 
u_healthy=.87
u_hindered=.715
u_fblind=.6
u_blind=.39

c_healthy_hc=0
c_hindered_hc=2852
c_fblind_hc=5778
c_blind_hc=7883
c_surgery=4876
c_gt=100000

#### Change Discount Rate ####
d.r <- 0.015                                 # equal discount of costs and QALYs by 3%
v.dwc <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weights for costs for each cycle based on discount rate d.r
v.dwe <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weights for effectiveness for each cycle based on discount rate d.r

#### 04 Define and initialize matrices and vectors ####
#### 04.1 Cohort trace ####
# create the markov trace matrix M capturing the proportion of the cohort in each state at each cycle
m.M_no_trt <- m.M_trt <- matrix(NA, 
                                nrow = n.t + 1, ncol = n.s,
                                dimnames = list(paste("cycle", 0:n.t, sep = " "), v.n))
m.M_no_trt 
head(m.M_no_trt) # show first 6 rows of the matrix 

# The cohort starts as "getting GT"
m.M_no_trt[1, ] <- m.M_trt[1, ] <- c(1000, 0, 0, 0, 0, 0,0,0,0,0,0) # initiate first cycle of cohort trace 

#### 04.2 Transition probability ARRAY #### -----------------------------------
# create transition probability array for NO treatment
a.P_notrt <- array(0,                                    # Create 3-D array
                   dim = c(n.s, n.s, n.t),
                   dimnames = list(v.n, v.n, 0:(n.t-1))) # name dimensions of the transition probability array
a.P_notrt 
getOption("max.print")
options(max.print = 1000)
# fill in the transition probability array
### From Healthy
v.n

### NO TREARTMENT transition probability array ###
### Start
a.P_notrt["Get_GT", "Hindered1",] <- 1
a.P_notrt["Get_GT", "Healthy",] <- 0

### From Healthy (no one will be in this state
a.P_notrt["Healthy", "Healthy",]  <- 1

### From Hindered1

a.P_notrt["Hindered1", "Hindered1",] <- p.Hindered1_Hindered1
a.P_notrt["Hindered1", "FBlind1",] <- p.Hindered1_FBlind1
a.P_notrt["Hindered1", "Blind2",] <- p.Hindered1_Blind2
a.P_notrt["Hindered1", "Death",] <- pDead

### From Hindered2

a.P_notrt["Hindered2", "Hindered2",] <- p.Hindered2_Hindered2
a.P_notrt["Hindered2", "FBlind1",] <- p.Hindered2_FBlind1
a.P_notrt["Hindered2", "FBlind2",] <- p.Hindered2_FBlind2
a.P_notrt["Hindered2", "Blind2",] <- p.Hindered2_Blind2
a.P_notrt["Hindered2", "Death",] <- pDead

### From FBlind1
a.P_notrt["FBlind1", "FBlind1",] <- p.FBlind1_FBlind1
a.P_notrt["FBlind1", "Blind1",] <- p.FBlind1_Blind1
a.P_notrt["FBlind1", "Blind2",] <- p.FBlind1_Blind2
a.P_notrt["FBlind1", "Death",] <- pDead

### From FBlind2
a.P_notrt["FBlind2", "FBlind2",] <- p.FBlind2_FBlind2
a.P_notrt["FBlind2", "Blind1",] <- p.FBlind2_Blind1
a.P_notrt["FBlind2", "Blind2",] <- p.FBlind2_Blind2
a.P_notrt["FBlind2", "Death",] <- pDead

### From Blind1
a.P_notrt["Blind1", "Blind1",] <- p.Blind1_Blind1 
a.P_notrt["Blind1", "Blind2",] <- p.Blind1_Blind2
a.P_notrt["Blind1", "Death",] <- pDead


### From Blind2
a.P_notrt["Blind2", "Blind2",] <- p.Blind2_Blind2 
a.P_notrt["Blind2", "Death",] <- pDead

a.P_notrt["GT_Healthy", "Death",]  <- 0

### From Death
a.P_notrt["Death", "Death",] <- 1

## OTHER
a.P_notrt["GT_Healthy","GT_Healthy",] <- 1
a.P_notrt["GT_Hindered1","GT_Hindered1",] <- 1

a.P_notrt


# Check if transition matrix is valid (i.e., each row should add up to 1)
valid <- apply(a.P_notrt, 3, function(x) sum(rowSums(x))==n.s)
if (!isTRUE(all.equal(as.numeric(sum(valid)), as.numeric(n.t)))) {
  stop("This is not a valid transition Matrix")
} 


# create transition probability matrix for treatment same as NO treatment
### TREATMENT TRANSITION PROBABILITY ###
a.P_trt <- a.P_notrt

array.diff <- a.P_trt - a.P_notrt
ifelse(array.diff == 0, "NULL", "DIFF") 

### Start 
a.P_trt["Get_GT", "Hindered1",] <- p.GT.fail
a.P_trt["Get_GT", "GT_Hindered1",] <- p.GT.success
a.P_trt["Get_GT", "Healthy",] <- 0
a.P_trt["Get_GT", "Healthy",] <- 0

### From GT_Hindered 
a.P_trt["GT_Hindered1", "GT_Hindered1",] <- p.GT.GTHindered1_GTHindered1
a.P_trt["GT_Hindered1", "Blind2",] <- p.GT.GTHindered1_Blind2
a.P_trt["GT_Hindered1", "Death",] <- pDead

a.P_trt
# create the markov trace matrix M capturing the proportion of the cohort in each state at each cycle
m.M_no_trt <- m.M_trt <- matrix(NA, 
                                nrow = n.t + 1, ncol = n.s,
                                dimnames = list(paste("cycle", 0:n.t, sep = " "), v.n))

head(m.M_no_trt) # show first 6 rows of the matrix 

# The cohort starts as healthy
m.M_no_trt[1, ] <- m.M_trt[1, ] <- c(1000, 0, 0, 0, 0, 0,0,0,0,0,0) # initiate the Markov trace 

#### 05 Run Markov model ####
for (t in 1:n.t){                                               # loop through the number of cycles
  m.M_no_trt[t + 1, ] <- t(m.M_no_trt[t, ]) %*% a.P_notrt[,, t] # estimate the Markov trace for cycle the next cycle (t + 1)
  m.M_trt[t + 1, ]    <- t(m.M_trt[t, ])    %*% a.P_trt[,, t]   # estimate the Markov trace for cycle the next cycle (t + 1)
} # close the loop

head(m.M_no_trt)  # show the first 6 lines of the matrix
head(m.M_trt)
m.M_no_trt

#### 06 Compute and Plot Epidemiological Outcomes ####
#### 06.1 Cohort trace #####
matplot(m.M_no_trt, type = 'l', 
        ylab = "Number of Patients",
        xlab = "Cycle (Age - 20 years)",
        main = "Cohort Trace - No Treatment")              # create a plot of the data

matplot(m.M_trt, type = 'l', 
        ylab = "Number of Patients",
        xlab = "Cycle (Age - 20 years)",
        main = "Cohort Trace - Gene Therapy")              # create a plot of the data

legend("topright", v.n, col = 1:n.s,lty = 1:n.s, bty = "n")  # add a legend to the graph
dev.off()

#### 07 Compute Cost-Effectiveness Outcomes RESTART HERE ####  
### Vectors with costs and utilities by treatment
v.u_no_trt <- v.u_trt <- c(0, u_healthy, u_healthy, u_hindered, u_hindered, u_hindered, u_fblind, u_fblind, u_blind, u_blind, 0)

v.c_no_trt <- c(0, c_healthy_hc, c_healthy_hc, c_hindered_hc, c_hindered_hc, c_hindered_hc, c_fblind_hc, c_fblind_hc, c_blind_hc, c_blind_hc, 0)
v.c_trt    <- c(c_surgery+c_gt, c_healthy_hc, c_healthy_hc, c_hindered_hc, c_hindered_hc, c_hindered_hc, c_fblind_hc, c_fblind_hc, c_blind_hc, c_blind_hc, 0)

#### 07.1 Mean Costs and QALYs for Treatment and NO Treatment ####
# must be in line with health states
# estimate mean QALYs ...
v.tu_no_trt <- m.M_no_trt %*% v.u_no_trt
v.tu_trt    <- m.M_trt %*% v.u_trt
v.tu_incr   <- v.tu_trt - v.tu_no_trt
# Net Monetary Benefit
net.monetary.benefit <- (v.tu_incr/1000)*50000

# ... and costs
v.tc_no_trt <- m.M_no_trt %*% v.c_no_trt
v.tc_trt    <- m.M_trt %*% v.c_trt

#### 07.2 Discounted Mean Costs and QALYs ####
### discount costs and QALYs
tu.d_no_trt <- t(v.tu_no_trt) %*% v.dwe  # 1x31 %*% 31x1 -> 1x1
tu.d_trt    <- t(v.tu_trt) %*% v.dwe

tc.d_no_trt <- t(v.tc_no_trt) %*% v.dwc
tc.d_trt    <- t(v.tc_trt)    %*% v.dwc
#sum(t(v.dwc*t(v.tc_trt)))
#tc.d_trt
#SOC incremental 
GT_discounted_cycle_utility_015_mid <- t(v.dwe*t(v.tu_trt))
GT_discounted_cycle_cost_015_mid <- t(v.dwc*t(v.tc_trt))
SOC_discounted_cycle_utility_015_mid <- t(v.dwe*t(v.tu_no_trt))
SOC_discounted_cycle_cost_015_mid <- t(v.dwc*t(v.tc_no_trt))

write.csv(GT_discounted_cycle_utility_015_mid, file = "GT_discounted_cycle_utility_015_mid")
write.csv(GT_discounted_cycle_cost_015_mid, file = "GT_discounted_cycle_cost_015_mid")
write.csv(SOC_discounted_cycle_utility_015_mid, file = "SOC_discounted_cycle_utility_015_mid")
write.csv(SOC_discounted_cycle_cost_015_mid, file = "SOC_discounted_cycle_cost_015_mid")

### Vector
v.tc.d <- c(tc.d_no_trt, tc.d_trt)
v.tu.d <- c(tu.d_no_trt, tu.d_trt)

#### Plot Traces ####
# mean QALY per Capita GT
matplot(v.tu_trt/1000, type = 'l', 
        ylab = "QALY",
        xlab = "Cycle (Age - 20 years)",
        main = "Mean QALY Value per Cycle - Gene Therapy")              # create a plot of the data
# mean QALY per Capita SOC
matplot(v.tu_no_trt/1000, type = 'l', 
        ylab = "QALY",
        xlab = "Cycle (Age - 20 years)",
        main = "Mean QALY Value per Cycle - SOC")              # create a plot of the data
#incremental mean QALY per Capita
matplot(v.tu_incr/1000, type = 'l', 
        ylab = "Incremental QALY",
        xlab = "Cycle (Age - 20 years)",
        main = "Incremental QALY (GT-SOC)")              # create a plot of the data
#costs SOC
matplot(v.tc_no_trt/1000, type = 'l', 
        ylab = "Cost (USD) no discount",
        xlab = "Cycle (Age - 20 years)",
        main = "Costs per year, SOC")              # create a plot of the data
# costs GT
matplot(v.tc_trt/1000, type = 'l', 
        ylab = "Cost (USD) no discount",
        xlab = "Cycle (Age - 20 years)",
        main = "Costs per year, GT")              # create a plot of the data


# Matrix with discounted costs and effectiveness
m.ce <- data.frame(Strategy = v.names.str,
                   Cost     = v.tc.d,
                   Effect   = v.tu.d)
m.ce
#### 05 Compute ICERs of Decision Tree ####
m.cea <- calculate_icers(m.ce)
m.cea

#### 06 Plot frontier of Decision Tree ####
front.ce <- getFrontier(m.ce, plot = F)
plot.frontier(CEmat = m.ce, frontier = front.ce)

#### SHORTER TIME HORIZON CALCULATIONS ####




