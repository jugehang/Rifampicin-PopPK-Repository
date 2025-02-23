##Author: Gehang Ju
##Address: Central South University, Xiangya Hospital, Hunan Province, China
##Email: 218101055@csu.edu.cn
##Aim: Establish PopPK repository of RIF

#setwd and library R packages------------
rm(list=ls())
#set working directory to current folder
curr.dir<-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(curr.dir)

# load R packages
library(tidyverse) # for data visualisation and manipulation
library(mrgsolve) # for simulation
library(rxode2) # for simulation for complex model
library(cowplot) # for combine figures
library(ggplot2)
library(ggpubr)
library(dplyr)
library(PopED)

#create fold for figure output
output_dir <- "simPK1017"
if (!file.exists(output_dir)) {dir.create (output_dir)}

# Vitrual patients #
# Adult: BW=70kg, FFM=56kg, AMT=450mg QD, N.Dose=8
# Children: AGE=2Years, PMA=2.75Years, BW=14kg, FFM=11kg, AMT=150mg QD, N.Dose=8

# 1.##Abdelgawad (2022)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 8.82*(FFM/43)^0.75*(BRC/6.0)^(-0.333)
# Vd  (L)    = 56.8*(FFM/43)
# Ka  (1/h)  = 1.38
# F   (%)    = 1
# MTT (h)    = 0.342
# NN         = 12 FIX
# BSV (CV%): CL = 42.4%
# BOV (CV%): Ka = 119%
# BOV (CV%): F  = 21.3%
# BOV (CV%): MTT= 93.8%
# prop.err= 17.2%
# add.err = 0.0234 FIX mg/L

# Mrgsolve code----
set.seed(92001)

cod1 <- '
$PARAM CL=8.82, VC=56.8, KA=1.38, MTT=0.342, NN=12, BIO=1,
$CMT DEPOT CENT
$GLOBAL 
int NDOSE = 0;
double dosetime[0];
double dose[450];
$MAIN
if(NEWIND < 2) NDOSE = 0; 

if(self.amt > 0 && self.cmt==1) {
 NDOSE = NDOSE + 1; 
 dosetime[NDOSE] = self.time;
 dose[NDOSE] = self.amt;
 }

F_DEPOT = 0;
double KTR = (NN+1)/MTT;
double NFAC = exp(lgamma(NN+1));
double KINPT = BIO * pow(KTR,(NN+1)) / NFAC; 

$ODE
double INPT = 0;
int i = 0;
while(i <= NDOSE) {
  double IPT = 0;
  if(SOLVERTIME >= dosetime[i]) {
    double delta = SOLVERTIME - dosetime[i];
    IPT = dose[i] * pow(delta, NN) * exp(-KTR * delta);  
  }
  INPT = INPT + IPT;
  ++i;
 }
dxdt_DEPOT = KINPT * INPT - KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT;

$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC MTT
'

mod1 <- mcode("transit", cod1, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff1 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]], 
                   KA = p[["KA"]], 
                   BIO = p[["BIO"]], 
                   MTT = p[["MTT"]],
                   NN = p[["NN"]])
  out <- mrgsim_q(mod1, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg1 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * (a[3]/43)^(0.75) * (a[4]/6)^(-0.333) * exp(b[1]),
    VC = bpop[2] * (a[3]/43),
    KA = bpop[3],
    BIO = bpop[4],
    MTT = bpop[5],
    NN = bpop[6],
    DOSE = a[1],
    TAU = a[2],
    FFM = a[3],
    BRC = a[4],
    AID = a[5]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model1 <- create.poped.database(
  ff_fun =ff1,
  fg_fun =fg1,
  fError_fun=feps1,
  bpop=c(CL=8.82,
         VC=56.8,
         KA=1.38,
         BIO=1,
         MTT=0.342,
         NN=12), 
  notfixed_bpop = c(1,1,0,0,0,0),
  d=c(CL=0.424^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,FFM=56,BRC=6,AID=1))

pl1 <- plot_model_prediction(poped_db_model1,
                             PI=T,
                             sample.times = F,
                             PRED = T) +
  labs(title = "Model1: Abdelgawad (2022)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl1

jpeg(filename = paste0(output_dir,"/1.Abdelgawad2022_Adult.jpg"), 
     width=1000, height=1000, res=300)
print(pl1)
dev.off()

# 2.##Chang (2015)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 6.10+((BMI/20.3)*6.22)
# Vd  (L)    = 48.0+(DM*16.2)
# Ka  (1/h)  = 1.31+(DM*1.56)
# BSV (CV%): CL = 53.7%
# BSV (CV%): Vd = 32.8%
# BSV (CV%): Ka = 49.9%
# prop.err= 12.0%
# add.err = 1.42 FIX mg/L

# Mrgsolve code----
cod2 <- '
$PARAM CL=31.4, VC=21.1, KA=1.7
$CMT DEPOT CENT
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT;
$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA
'
mod2 <- mcode("optim", cod2, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff2 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]],
                   KA = p[["KA"]])
  out <- mrgsim_q(mod2, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg2 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = (bpop[1] + (a[3]/20.3)*6.22) * exp(b[1]),
    VC = (bpop[2] + (a[4]*16.2)) * exp(b[2]),
    KA = (bpop[3] + (a[4]*1.56)) * exp(b[3]),
    DOSE = a[1],
    TAU = a[2],
    BMI = a[3],
    DM = a[4]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model2 <- create.poped.database(
  ff_fun =ff2,
  fg_fun =fg2,
  fError_fun=feps1,
  bpop=c(CL=6.10,
         VC=48.0,
         KA=1.31), 
  notfixed_bpop = c(1,1,1),
  d=c(CL=0.537^2,VC=0.328^2,KA=0.499^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,BMI=20.3,DM=0,AID=2))

pl2 <- plot_model_prediction(poped_db_model2,
                             PI=T,
                             sample.times = F,
                             PRED = T) +
  labs(title = "Model2: Chang (2015)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl2

jpeg(filename = paste0(output_dir,"/2.Chang2015_Adult.jpg"), 
     width=1000, height=1000, res=300)
print(pl2)
dev.off()

# 3.##Chirehwa (2016)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL0,int  (L/h)  = 93.2*(FFM/42)^0.75
# CLss,int  (L/h)  = 176*(FFM/42)^0.75
# t1/2ind (day) = 4.5 = 108 h
# Vd  (L)    = 50.1*(FFM/42)
# Ka  (1/h)  = 1.96
# MTT (h) = 0.71
# NN = 19.3
# F = 1
# VH = 1 FIX
# QH = 50 FIX
# fu = 0.2 FIX
# Km (mg/L) = 3.35
# BSV (CV%): CL0,int = 22.5%
# BSV (CV%): Vd = 14.2%
# prop.err= 10.8%
# add.err = 0.064 mg/L
# Mrgsolve code----
set.seed(92001)

cod3 <- '
$PARAM CL0=93.2, CLss=176, VC=50.1, KA=1.96, MTT=0.71, NN=19.3, BIO=1, VH=1, QH=50, fu=0.2, Km=3.35, th=4.5, TM=168
$CMT DEPOT LIVER CENT
$GLOBAL 
int NDOSE = 0;
double dosetime[0];
double dose[450];

$MAIN
if(NEWIND < 2) NDOSE = 0; 

if(self.amt > 0 && self.cmt==1) {
 NDOSE = NDOSE + 1; 
 dosetime[NDOSE] = self.time;
 dose[NDOSE] = self.amt;
 }

F_DEPOT = 0;
double KTR = (NN+1)/MTT;
double NFAC = exp(lgamma(NN+1));
double KINPT = BIO * pow(KTR,(NN+1)) / NFAC;

$ODE
double INPT = 0;
int i = 0;
while(i <= NDOSE) {
  double IPT = 0;
  if(SOLVERTIME >= dosetime[i]) {
    double delta = SOLVERTIME - dosetime[i];
    IPT = dose[i] * pow(delta, NN) * exp(-KTR * delta);  
  }
  INPT = INPT + IPT;
  ++i;
  }

dxdt_DEPOT = KINPT * INPT - KA*DEPOT;
dxdt_LIVER = KA*DEPOT - QH*(1-((((CL0+(CLss-CL0)*(1-exp(-log(2)*TM/24/th)))*CH)/(CH+Km))*fu)/((((CL0+(CLss-CL0)*(1-exp(-log(2)*TM/24/th)))*CH)/(CH+Km))*fu+QH))/VH*LIVER - QH*((((CL0+(CLss-CL0)*(1-exp(-log(2)*TM/24/th)))*CH)/(CH+Km))*fu)/((((CL0+(CLss-CL0)*(1-exp(-log(2)*TM/24/th)))*CH)/(CH+Km))*fu+QH)/VH*LIVER + QH/VC*CENT;
dxdt_CENT  = QH*(1-((((CL0+(CLss-CL0)*(1-exp(-log(2)*TM/24/th)))*CH)/(CH+Km))*fu)/((((CL0+(CLss-CL0)*(1-exp(-log(2)*TM/24/th)))*CH)/(CH+Km))*fu+QH))/VH*LIVER - QH/VC*CENT;

$TABLE 
double CP  = CENT/VC;
double CH  = LIVER/VH;
$CAPTURE CP
'

mod3 <- mcode("transit", cod3, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff3 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   TM=time,
                   CL0 = p[["CL0"]], 
                   CLss = p[["CLss"]], 
                   VC = p[["VC"]], 
                   KA = p[["KA"]])
  out <- mrgsim_q(mod3, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg3 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL0 = bpop[1] * (a[3]/42)^(0.75) * exp(b[1]),
    CLss = bpop[2] * (a[3]/42)^(0.75),
    VC = bpop[3] * (a[3]/42) * exp(b[2]),
    KA = bpop[4],
    DOSE = a[1],
    TAU = a[2],
    FFM = a[3],
    AID = a[4]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model3 <- create.poped.database(
  ff_fun =ff3,
  fg_fun =fg3,
  fError_fun=feps1,
  bpop=c(CL0=93.2,
         CLss=176,
         VC=50.1,
         KA=1.96), 
  notfixed_bpop = c(1,1,1,1),
  d=c(CL0=0.225^2,VC=0.142^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,FFM=56,AID=3))

pl3 <- plot_model_prediction(poped_db_model3,
                             PI=T,
                             sample.times = F,
                             PRED = T) +
  labs(title = "Model3: Chirehwa (2016)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl3

jpeg(filename = paste0(output_dir,"/3.Chirehwa2016_Adult.jpg"), 
     width=1000, height=1000, res=300)
print(pl3)
dev.off()

# 4.##Denti (2016)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 16.2*(BW/66)^0.75*(1-Prgancy*0.14)
# Vd  (L)    = 43.3*(BW/66)
# Ka  (1/h)  = 1.67
# F   (%)    = 1
# MTT (h)    = 1.31
# NN         = 54.6
# BSV (CV%): CL = 30.4%
# prop.err= 13.1%
# add.err = 0.0585 FIX mg/L
# Mrgsolve code----
set.seed(92001)

cod4 <- '
$PARAM CL=16.2, VC=43.3, KA=1.67, MTT=1.31, NN=54.6, BIO=1,
$CMT DEPOT CENT
$GLOBAL 
int NDOSE = 0;
double dosetime[0];
double dose[450];
$MAIN
if(NEWIND < 2) NDOSE = 0; 

if(self.amt > 0 && self.cmt==1) {
 NDOSE = NDOSE + 1; 
 dosetime[NDOSE] = self.time;
 dose[NDOSE] = self.amt;
 }

F_DEPOT = 0;
double KTR = (NN+1)/MTT;
double NFAC = exp(lgamma(NN+1));
double KINPT = BIO * pow(KTR,(NN+1)) / NFAC; 

$ODE
double INPT = 0;
int i = 0;
while(i <= NDOSE) {
  double IPT = 0;
  if(SOLVERTIME >= dosetime[i]) {
    double delta = SOLVERTIME - dosetime[i];
    IPT = dose[i] * pow(delta, NN) * exp(-KTR * delta);  
  }
  INPT = INPT + IPT;
  ++i;
 }
dxdt_DEPOT = KINPT * INPT - KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT;

$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC MTT
'

mod4 <- mcode("transit", cod4, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff4 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]], 
                   KA = p[["KA"]], 
                   BIO = p[["BIO"]], 
                   MTT = p[["MTT"]],
                   NN = p[["NN"]])
  out <- mrgsim_q(mod4, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg4 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * (a[3]/66)^(0.75) * exp(b[1]),
    VC = bpop[2] * (a[3]/66),
    KA = bpop[3],
    BIO = bpop[4],
    MTT = bpop[5],
    NN = bpop[6],
    DOSE = a[1],
    TAU = a[2],
    FFM = a[3],
    AID = a[4]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model4 <- create.poped.database(
  ff_fun =ff4,
  fg_fun =fg4,
  fError_fun=feps1,
  bpop=c(CL=16.2,
         VC=43.3,
         KA=1.67,
         BIO=1,
         MTT=1.31,
         NN=54.6), 
  notfixed_bpop = c(1,1,1,0,0,0),
  d=c(CL=0.304^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,FFM=56,AID=4))

pl4 <- plot_model_prediction(poped_db_model4,
                             PI=T,
                             sample.times = F,
                             PRED = T) +
  labs(title = "Model4: Denti (2016)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl4

jpeg(filename = paste0(output_dir,"/4.Denti2016_Adult.jpg"), 
     width=1000, height=1000, res=300)
print(pl4)
dev.off()


# 5.##Gao (2021)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 9.4*(BW/50)^0.76*(1-0.15*SEX)
# Vd  (L)    = 37.0*(BW/50)^0.66
# Ka  (1/h)  = 0.82
# BSV (CV%): CL = 28.5%
# BSV (CV%): Vd = 14.6%
# BSV (CV%): Ka = 14.2%
# prop.err= 6.3%
# add.err = 0.07 mg/L

# Mrgsolve code----
cod5 <- '
$PARAM CL=31.4, VC=21.1, KA=1.7
$CMT DEPOT CENT
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT;
$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA
'
mod5 <- mcode("optim", cod5, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff5 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]],
                   KA = p[["KA"]])
  out <- mrgsim_q(mod5, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg5 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * (a[3]/50)^0.75 * exp(b[1]),
    VC = bpop[2] * (a[3]/50) * exp(b[2]),
    KA = bpop[3] * exp(b[3]),
    DOSE = a[1],
    TAU = a[2],
    BW = a[3],
    ADI = a[4]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model5 <- create.poped.database(
  ff_fun =ff5,
  fg_fun =fg5,
  fError_fun=feps1,
  bpop=c(CL=31.4,
         VC=21.1,
         KA=1.7), 
  notfixed_bpop = c(1,1,1),
  d=c(CL=0.285^2,VC=0.146^2,KA=0.142^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,BW=70,AID=5))

pl5 <- plot_model_prediction(poped_db_model5,
                             PI=T,
                             sample.times = F,
                             PRED = T) +
  labs(title = "Model5: Gao (2021)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl5

jpeg(filename = paste0(output_dir,"/5.Gao2021_Adult.jpg"), 
     width=1000, height=1000, res=300)
print(pl5)
dev.off()

# 6.##Jeremiah (2014)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL (L/h)   = CL7+(CLss-CL7)*(1-e(-log(2)*t/t1/2))
# CL7        = 13.9*(FFM/43)^0.75
# CLss       = 16.5*(FFM/43)^0.75
# t1/2 (day) = 6 FIX
# Vd  (L)    = 55.8*(FFM/43)
# Ka  (1/h)  = 1.77
# F   (%)    = 1
# MTT (h)    = 1.50
# NN         = 27.6
# BSV (CV%): CL = 24%
# prop.err= 13.7%
# add.err = 0.0417 mg/L
# Mrgsolve code----
set.seed(92001)

cod6 <- '
$PARAM 
CL7=13.9, CLss=16.5, VC=43.3, KA=1.67, MTT=1.31, NN=54.6, BIO=1, FFM=56, th=144, TM=168
$CMT DEPOT CENT
$GLOBAL 
int NDOSE = 0;
double dosetime[0];
double dose[450];
$MAIN
if(NEWIND < 2) NDOSE = 0; 

if(self.amt > 0 && self.cmt==1) {
 NDOSE = NDOSE + 1; 
 dosetime[NDOSE] = self.time;
 dose[NDOSE] = self.amt;
 }

F_DEPOT = 0;
double KTR = (NN+1)/MTT;
double NFAC = exp(lgamma(NN+1));
double KINPT = BIO * pow(KTR,(NN+1)) / NFAC; 

$ODE
double INPT = 0;
int i = 0;
while(i <= NDOSE) {
  double IPT = 0;
  if(SOLVERTIME >= dosetime[i]) {
    double delta = SOLVERTIME - dosetime[i];
    IPT = dose[i] * pow(delta, NN) * exp(-KTR * delta);  
  }
  INPT = INPT + IPT;
  ++i;
  }
dxdt_DEPOT = KINPT * INPT - KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL7+(CLss-CL7)*(1-exp(-log(2)*TM/th)))/VC*CENT;

$TABLE double CP  = CENT/VC;
$CAPTURE CP CL7 CLss VC MTT TM
'

mod6 <- mcode("transit", cod6, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff6 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   TM = time,
                   CL7 = p[["CL7"]], 
                   CLss = p[["CLss"]], 
                   VC = p[["VC"]], 
                   KA = p[["KA"]], 
                   BIO = p[["BIO"]], 
                   MTT = p[["MTT"]],
                   NN = p[["NN"]])
  out <- mrgsim_q(mod6, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg6 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL7 = bpop[1] * (a[3]/43)^(0.75) * exp(b[1]),
    CLss= bpop[2] * (a[3]/43)^(0.75) * exp(b[2]),
    VC = bpop[3] * (a[3]/43),
    KA = bpop[4],
    BIO = bpop[5],
    MTT = bpop[6],
    NN = bpop[7],
    DOSE = a[1],
    TAU = a[2],
    FFM = a[3],
    AID = a[4]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model6 <- create.poped.database(
  ff_fun =ff6,
  fg_fun =fg6,
  fError_fun=feps1,
  bpop=c(CL7=13.9,
         CLss=16.5,
         VC=55.8,
         KA=1.77,
         BIO=1,
         MTT=1.5,
         NN=27.6), 
  notfixed_bpop = c(1,1,1,1,0,0,0),
  d=c(CL7=0.24^2,CLss=0.24^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,FFM=56,AID=6))

pl6 <- plot_model_prediction(poped_db_model6,
                             PI=T,
                             sample.times = F,
                             PRED = T) +
  labs(title = "Model6: Jeremiah (2014)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl6

jpeg(filename = paste0(output_dir,"/6.Jeremiah2014_Adult.jpg"), 
     width=1000, height=1000, res=300)
print(pl6)
dev.off()


# 7.##Jing (2016)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 4.02
# Vd  (L)    = 57.8
# Ka  (1/h)  = 1.61 FIX
# BSV (CV%): CL = 64.5%
# BSV (CV%): Vd = 20.9%
# add.err = 6.55 mg/L

# Mrgsolve code----
cod7 <- '
$PARAM CL=4.02, VC=57.8, KA=1.61
$CMT DEPOT CENT
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT;
$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA
'
mod7 <- mcode("optim", cod7, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff7 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]],
                   KA = p[["KA"]])
  out <- mrgsim_q(mod7, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg7 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * exp(b[1]),
    VC = bpop[2] * exp(b[2]),
    KA = bpop[3],
    DOSE = a[1],
    TAU = a[2],
    ADI = a[3]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model7 <- create.poped.database(
  ff_fun =ff7,
  fg_fun =fg7,
  fError_fun=feps1,
  bpop=c(CL=4.02,
         VC=57.8,
         KA=1.61), 
  notfixed_bpop = c(1,1,0),
  d=c(CL=0.645^2,VC=0.209^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,AID=7))

pl7 <- plot_model_prediction(poped_db_model7,
                             PI=T,
                             sample.times = F,
                             PRED = T) +
  labs(title = "Model7: Jing (2016)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl7

jpeg(filename = paste0(output_dir,"/7.Jing2016_Adult.jpg"), 
     width=1000, height=1000, res=300)
print(pl7)
dev.off()

# 8.##Karballaei (2022)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 0.08
# Vd  (L)    = 0.68
# Ka1 (1/h)  = 1.1
# Ka2 (1/h)  = 0.46
# Fr  (%)    = 68%
# Tlag2 (h)  = 2.92
# BSV (CV%): CL = 52.0%
# BSV (CV%): Vd = 9.7%
# BSV (CV%): Ka1 = 85.0%
# BSV (CV%): Ka2 = 46.0%
# BSV (CV%): Fr  = 357%
# BSV (CV%): Tlag2 = 299%
# prop.err= 49%
# add.err = 0.15 mg/L

# rxode2 code----
mod8 <- rxode2 ({
  alagDepot = TLAG2
  d/dt(DEPOT1) = -KA1*DEPOT1;
  d/dt(DEPOT2) = -KA2*DEPOT2;
  alag(DEPOT2) = alagDepot
  d/dt(CENT) = Fr*KA1*DEPOT1 + (1-Fr)*KA2*DEPOT2  - (CL/VC)*CENT;
  CP=CENT/VC;
})

## -- Model structure definition function----
ff8 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)
  et(amount.units = "mg", time.units = "hours") %>%
    add.dosing(dosing.to = "DEPOT1",
               dose      = p[["DOSE"]],
               dosing.interval = p[["TAU"]],
               nbr.doses = 8,
               start.time = 0) %>%
    add.dosing(dosing.to = "DEPOT2",
               dose      = p[["DOSE"]],
               dosing.interval = p[["TAU"]],
               nbr.doses = 8,
               start.time = 0) %>%
    et(times_xt) -> data
  
  out <- rxSolve(mod8, p, data, atol=1e-8, rtol=1e-8,maxsteps=50000,
                 returnType="data.frame")
  y <-  out$CP[match(times_xt,out$time)]
  return(list(y=matrix(y,ncol=1),poped.db=poped.db))
}

## -- parameter definition function----
fg8 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * a[3] * exp(b[1]),
    VC = bpop[2] * a[3] * exp(b[2]),
    KA1 = bpop[3] * exp(b[3]),
    KA2 = bpop[4] * exp(b[4]),
    Fr = bpop[5] * exp(b[5]),
    TLAG2 = bpop[6] * exp(b[6]),
    DOSE = a[1],
    TAU = a[2],
    BW = a[3],
    ADI = a[4]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model8 <- create.poped.database(
  ff_fun =ff8,
  fg_fun =fg8,
  fError_fun=feps1,
  bpop=c(CL=0.08,
         VC=0.68,
         KA1=1.1,
         KA2=0.46,
         Fr=0.68,
         TLAG2=2.92), 
  notfixed_bpop = c(1,1,1,1,0,0),
  d=c(CL=0.52^2,VC=0.097^2,KA1=0.85^2,KA2=0.46^2,Fr=0^2,TLAG2=0^2), #Fr=3.57^2,TLAG2=2.99^2
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,BW=70,AID=8))

pl8 <- plot_model_prediction(poped_db_model8,
                             PI=T,
                             sample.times = F,
                             PRED = T) +
  labs(title = "Model8: Karballaei (2022)",
       x = "Time from last dose (h)",
       y = "Concentration (mg/L)") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl8

jpeg(filename = paste0(output_dir,"/8.Karballaei2022_Adult.jpg"), 
     width=2000, height=2000, res=300)
print(pl8)
dev.off()

# 9.##Kim (2021)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 11.4*(BW/60)^1.14
# Vd  (L)    = 17.8
# Q   (L/h)  = 2.78
# Vp  (L)    = 80.7
# Ka  (1/h)  = 0.436 FIX
# BSV (CV%): CL = 64.2%
# BSV (CV%): Vd = 70.2%
# BSV correlation: CL_Vd = 0.927
# prop.err= 24.3.0%

# Mrgsolve code----
set.seed(92001)

cod9 <- '
$PARAM CL=11.4, VC=17.8, Q=2.78, VP=80.7, KA=1.7
$CMT DEPOT CENT PERI
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT - (Q/VC)*CENT + (Q/VP)*PERI;
dxdt_PERI = (Q/VC)*CENT - (Q/VP)*PERI;
$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC Q VP KA
'
mod9 <- mcode("optim", cod9, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff9 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]],
                   Q  = p[["Q"]], 
                   VP = p[["VP"]],
                   KA = p[["KA"]])
  out <- mrgsim_q(mod9, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg9 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * (a[3]/60)^1.14 * exp(b[1]),
    VC = bpop[2] * exp(b[2]),
    Q  = bpop[3],
    VP = bpop[4],
    KA = bpop[5],
    DOSE = a[1],
    TAU = a[2],
    BMI = a[3],
    DM = a[4]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps2 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1])
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model9 <- create.poped.database(
  ff_fun =ff9,
  fg_fun =fg9,
  fError_fun=feps.add.prop,
  bpop=c(CL=11.4,
         VC=17.8,
         Q=2.78,
         VP=80.7,
         KA=1.31), 
  notfixed_bpop = c(1,1,1,1,1),
  d=c(CL=0.642^2,VC=0.702^2),
  covd = 0.418,
  sigma=c(prop=0,add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,BW=70,AID=9))

pl9 <- plot_model_prediction(poped_db_model9,
                             PI=T,
                             sample.times = F,
                             PRED = T) +
  labs(title = "Model9: Kim (2021)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl9

jpeg(filename = paste0(output_dir,"/9.Kim2021_Adult.jpg"), 
     width=1000, height=1000, res=300)
print(pl9)
dev.off()

# 10.##Kloprogge (2020)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 16.9*(BW/70)^0.75*(1+0.183*SEX)
# Vd  (L)    = 31.3*(BW/70)
# Ka  (1/h)  = 0.277 FIX
# F   (%)    = 1
# MTT (h)    = 0.326 FIX
# NN         = 1.5 FIX
# BSV (CV%): CL = 42.4%
# BOV (CV%): Ka = 119%
# BOV (CV%): F  = 21.3%
# BOV (CV%): MTT= 93.8%
# prop.err= 19.4%

# Mrgsolve code----
set.seed(92001)

cod10 <- '
$PARAM CL=16.9, VC=31.3, KA=0.277, MTT=0.326, NN=1.5, BIO=1,
$CMT DEPOT CENT
$GLOBAL 
int NDOSE = 0;
double dosetime[24];
double dose[450];

$MAIN
if(NEWIND < 2) NDOSE = 0; 

if(self.amt > 0 && self.cmt==1) {
 NDOSE = NDOSE + 1; 
 dosetime[NDOSE] = self.time;
 dose[NDOSE] = self.amt;
 }

F_DEPOT = 0;
double KTR = (NN+1)/MTT;
double NFAC = exp(lgamma(NN+1));
double KINPT = BIO * pow(KTR,(NN+1)) / NFAC; 

$ODE
double INPT = 0;
int i = 0;
while(i <= NDOSE) {
  double IPT = 0;
  if(SOLVERTIME >= dosetime[i]) {
    double delta = SOLVERTIME - dosetime[i];
    IPT = dose[i] * pow(delta, NN) * exp(-KTR * delta);  
  }
  INPT = INPT + IPT;
  ++i;
 }
dxdt_DEPOT = KINPT * INPT - KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT;

$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC MTT
'

mod10 <- mcode("transit", cod10, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff10 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]], 
                   KA = p[["KA"]], 
                   BIO = p[["BIO"]], 
                   MTT = p[["MTT"]],
                   NN = p[["NN"]])
  out <- mrgsim_q(mod10, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg10 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * (a[3]/70)^(0.75) * (1+0.183*a[4]) * exp(b[1]),
    VC = bpop[2] * (a[3]/70) * exp(b[2]),
    KA = bpop[3],
    BIO = bpop[4],
    MTT = bpop[5] * exp(b[3]),
    NN = bpop[6],
    DOSE = a[1],
    TAU = a[2],
    BW = a[3],
    SEX = a[4],
    AID = a[5]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model10 <- create.poped.database(
  ff_fun =ff10,
  fg_fun =fg10,
  fError_fun=feps1,
  bpop=c(CL=16.9,
         VC=31.3,
         KA=0.277,
         BIO=1,
         MTT=0.326,
         NN=1.5), 
  notfixed_bpop = c(1,1,0,0,0,0),
  d=c(CL=0.4008^2,VC=0.8452^2,MTT=0.2705^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,BW=70,SEX=0,AID=10))

pl10 <- plot_model_prediction(poped_db_model10,
                             PI=T,
                             sample.times = F,
                             PRED = T) +
  labs(title = "Model10: Kloprogge (2020) Female",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl10

jpeg(filename = paste0(output_dir,"/10.Kloprogge2020_Adult_F.jpg"), 
     width=1000, height=1000, res=300)
print(pl10)
dev.off()

poped_db_model101 <- create.poped.database(
  ff_fun =ff10,
  fg_fun =fg10,
  fError_fun=feps1,
  bpop=c(CL=16.9,
         VC=31.3,
         KA=0.277,
         BIO=1,
         MTT=0.326,
         NN=1.5), 
  notfixed_bpop = c(1,1,0,0,0,0),
  d=c(CL=0.4008^2,VC=0.8452^2,MTT=0.2705^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,BW=70,SEX=1,AID=10))

pl101 <- plot_model_prediction(poped_db_model101,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model10: Kloprogge (2020) Male",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl101

jpeg(filename = paste0(output_dir,"/10.Kloprogge2020_Adult_M.jpg"), 
     width=1000, height=1000, res=300)
print(pl101)
dev.off()

# 11.##Marsot (2017)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 13.7-8.6*Fusidic acid
# Vd  (L)    = 61.1-37.3*Fusidic acid
# Ka  (1/h)  = 1.15 FIX
# BSV (CV%): CL = 53.1%
# BSV (CV%): Vd = 34.9%
# add.err = 2.256 FIX mg/L

# Mrgsolve code----
cod11 <- '
$PARAM CL=13.7, VC=61.1, KA=1.15
$CMT DEPOT CENT
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT;
$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA
'
mod11 <- mcode("optim", cod11, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff11 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]],
                   KA = p[["KA"]])
  out <- mrgsim_q(mod11, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg11 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = (bpop[1] - a[3]*8.6) * exp(b[1]),
    VC = (bpop[2] + a[3]*37.3) * exp(b[2]),
    KA = bpop[3],
    DOSE = a[1],
    TAU = a[2],
    FA = a[3],
    AID = a[4]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model11 <- create.poped.database(
  ff_fun =ff11,
  fg_fun =fg11,
  fError_fun=feps1,
  bpop=c(CL=13.7,
         VC=61.1,
         KA=1.15), 
  notfixed_bpop = c(1,1,0),
  d=c(CL=0.531^2,VC=0.349^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,FA=0,AID=11))

pl11 <- plot_model_prediction(poped_db_model11,
                             PI=T,
                             sample.times = F,
                             PRED = T) +
  labs(title = "Model11: Marsot (2017)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl11

jpeg(filename = paste0(output_dir,"/11.Masort2017_Adult.jpg"), 
     width=1000, height=1000, res=300)
print(pl11)
dev.off()

# 12.##Medellin (2020)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 5.96
# Vd  (L)    = 0.7*BW
# Ka  (1/h)  = 1.24
# D0  (h)    = 0.62
# Tlag (h)   = 0.24
# BSV (CV%): CL =38.5%
# BSV (CV%): Vd = 26.8%
# BSV (CV%): Ka = 110.5%
# BSV (CV%): D0 = 131.1%
# add.err = 2.256 FIX mg/L

# Mrgsolve code----
cod12 <- '
$PARAM CL=5.96, VC=49, KA=1.24, DUR=0.62, TLAG=0.24
$CMT DEPOT CENT
$MAIN
D_CENT = DUR;
ALAG_DEPOT = TLAG;
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT;
$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA
'

mod12<- mcode("optim", cod12, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff12 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]],
                   KA = p[["KA"]])
  out <- mrgsim_q(mod12, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg12 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * exp(b[1]),
    VC = bpop[2] * a[3] * exp(b[2]),
    KA = bpop[3] * exp(b[3]),
    DOSE = a[1],
    TAU = a[2],
    BW = a[3],
    AID = a[4]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model12 <- create.poped.database(
  ff_fun =ff12,
  fg_fun =fg12,
  fError_fun=feps1,
  bpop=c(CL=5.96,
         VC=0.7,
         KA=1.24), 
  notfixed_bpop = c(1,1,1),
  d=c(CL=0.385^2,VC=0.268^2,KA=1.105^2),
  sigma=c(prop=0,add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,BW=70,AID=14))

pl12 <- plot_model_prediction(poped_db_model12,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model12: Medellin (2020)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl12

jpeg(filename = paste0(output_dir,"/12.Medellin2020_Adult.jpg"), 
     width=1000, height=1000, res=300)
print(pl12)
dev.off()

# 13.##Milán (2013)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 8.17*(1+SEX*0.4)
# Vd  (L)    = 50.1*(1+SEX*0.29)
# Ka  (1/h)  = 2.70
# BSV (CV%): CL = 31.9%
# BSV (CV%): Vd = 16.7%
# BSV (CV%): Ka = 92.9%
# prop.err= 7.74%

# Mrgsolve code----
set.seed(92001)

cod13 <- '
$PARAM CL=8.17, VC=50.1, KA=2.7
$CMT DEPOT CENT
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT;
$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA
'
mod13 <- mcode("optim", cod13, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff13 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]],
                   KA = p[["KA"]])
  out <- mrgsim_q(mod13, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg13 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1]*(1+a[3]*0.4) * exp(b[1]),
    VC = bpop[2]*(1+a[3]*0.29) * exp(b[2]),
    KA = bpop[3] * exp(b[3]),
    DOSE = a[1],
    TAU = a[2],
    SEX = a[3],
    AID = a[4]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model13 <- create.poped.database(
  ff_fun =ff13,
  fg_fun =fg13,
  fError_fun=feps1,
  bpop=c(CL=8.17,
         VC=50.1,
         KA=2.70), 
  notfixed_bpop = c(1,1,1),
  d=c(CL=0.319^2,VC=0.167^2,KA=0.929^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,SEX=0,AID=13))

pl13 <- plot_model_prediction(poped_db_model13,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model13: Milán (2013) Female",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl13

jpeg(filename = paste0(output_dir,"/13.Milán2013_Adult_F.jpg"), 
     width=1000, height=1000, res=300)
print(pl13)
dev.off()

poped_db_model131 <- create.poped.database(
  ff_fun =ff13,
  fg_fun =fg13,
  fError_fun=feps1,
  bpop=c(CL=8.17,
         VC=50.1,
         KA=2.70), 
  notfixed_bpop = c(1,1,1),
  d=c(CL=0.319^2,VC=0.167^2,KA=0.929^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,SEX=1,AID=13))

pl131 <- plot_model_prediction(poped_db_model131,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model13: Milán (2013) Male",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl131

jpeg(filename = paste0(output_dir,"/13.Milán2013_Adult_M.jpg"), 
     width=1000, height=1000, res=300)
print(pl131)
dev.off()

# 14.##Mukonzo (2020)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 19.8831
# Vd  (L)    = 0.5383
# Q   (L/h)  = 19.6854
# Vp  (L)    = 19.3284
# Ka  (1/h)  = 0.4682
# Tlag (h)   = 0.7748
# BSV (CV%): CL = 38.22%
# BSV (CV%): Vd = 255.65%
# BSV (CV%): Ka = 9.22%
# prop.err= 13.86%
# add.err= 0.0032

# Mrgsolve code----
set.seed(92001)

cod14 <- '
$PARAM CL=19.8831, VC=0.5383, KA=0.4682, Q=19.6854, VP=19.3284, TLAG=0.7748
$CMT DEPOT CENT PERI
$MAIN
ALAG_DEPOT = TLAG;
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT - (Q/VC)*CENT + (Q/VP)*PERI;
dxdt_PERI = (Q/VC)*CENT - (Q/VP)*PERI;
$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA Q VP
'
mod14<- mcode("optim", cod14, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff14 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]],
                   Q  = p[["Q"]], 
                   VP = p[["VP"]],
                   KA = p[["KA"]])
  out <- mrgsim_q(mod14, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg14 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * exp(b[1]),
    VC = bpop[2] * exp(b[2]),
    Q  = bpop[3],
    VP = bpop[4],
    KA = bpop[5] * exp(b[3]),
    DOSE = a[1],
    TAU = a[2],
    AID = a[3]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model14 <- create.poped.database(
  ff_fun =ff14,
  fg_fun =fg14,
  fError_fun=feps1,
  bpop=c(CL=19.8831,
         VC=0.5383,
         Q=19.6854,
         VP=10.3284,
         KA=0.4682), 
  notfixed_bpop = c(1,1,1,1,1),
  d=c(CL=0.3822^2,VC=2.5565^2,KA=0.0922^2),
  sigma=c(prop=0,add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,AID=14))

pl14 <- plot_model_prediction(poped_db_model14,
                             PI=T,
                             sample.times = F,
                             PRED = T) +
  labs(title = "Model14: Mukonzo (2020)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl14

jpeg(filename = paste0(output_dir,"/14.Mukonzo2020_Adult.jpg"), 
     width=1000, height=1000, res=300)
print(pl14)
dev.off()

# 15.##Nishimura (2020)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 15.2*(BW/51.5)^0.978
# Vd  (L)    = 2.09
# Ka  (1/h)  = 0.248*(1+FOOD*0.21)
# D0  (h)    = 1.58*FOOD
# F          = 1*(1+FOOD*0.32)
# Tlag (h)   = 1.43
# BSV (CV%): CL =38.92%
# BSV (CV%): Vd = 151.66%
# BSV (CV%): Ka = 24.98%
# BSV (CV%): D0 = 78.23%
# BSV (CV%): Tlag = 12.96%
# BSV (CV%): F = 0.00063%
# prop.err = 4.95%

# Mrgsolve code----
cod15 <- '
$PARAM CL=15.2, VC=2.09, KA=0.248, DUR=1.58, TLAG=1.43, BIO=1
$CMT DEPOT CENT
$MAIN
D_CENT = DUR;
ALAG_DEPOT = TLAG;
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = BIO*KA*DEPOT - (CL/VC)*CENT;
$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA DUR BIO
'

mod15<- mcode("optim", cod15, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff15 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]],
                   KA = p[["KA"]],
                   DUR = p[["DUR"]],
                   BIO = p[["BIO"]])
  out <- mrgsim_q(mod15, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg15 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * (a[3]/51.5)^0.978 * exp(b[1]),
    VC = bpop[2] * exp(b[2]),
    KA = bpop[3] * (a[4]*0.21+1) * exp(b[3]),
    DUR = bpop[4] * a[4] * exp(b[4]),
    BIO = bpop[5] * (a[4]*0.32+1) * exp(b[5]),
    DOSE = a[1],
    TAU = a[2],
    BW = a[3],
    FOOD = a[4],
    AID = a[5]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model15 <- create.poped.database(
  ff_fun =ff15,
  fg_fun =fg15,
  fError_fun=feps1,
  bpop=c(CL=15.2,
         VC=2.09,
         KA=0.248,
         DUR=1.58,
         BIO=1), 
  notfixed_bpop = c(1,1,1,0,0),
  d=c(CL=0.3892^2,VC=1.5166^2,KA=0.2498^2,DUR=0.7823^2,BIO=0.00063^2),
  sigma=c(prop=0,add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,BW=70,FOOD=0,AID=15))

pl15 <- plot_model_prediction(poped_db_model15,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model15: Nishimura (2020) Fasted",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl15

jpeg(filename = paste0(output_dir,"/15.Nishimura2020_Adult_NF.jpg"), 
     width=1000, height=1000, res=300)
print(pl15)
dev.off()

poped_db_model151 <- create.poped.database(
  ff_fun =ff15,
  fg_fun =fg15,
  fError_fun=feps1,
  bpop=c(CL=15.2,
         VC=2.09,
         KA=0.248,
         DUR=1.58,
         BIO=1), 
  notfixed_bpop = c(1,1,1,0,0),
  d=c(CL=0.3892^2,VC=1.5166^2,KA=0.2498^2,DUR=0.7823^2,BIO=0.00063^2),
  sigma=c(prop=0,add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,BW=70,FOOD=1,AID=15))

pl151 <- plot_model_prediction(poped_db_model151,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model15: Nishimura (2020) Fed",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl151

jpeg(filename = paste0(output_dir,"/15.Nishimura2020_Adult_F.jpg"), 
     width=1000, height=1000, res=300)
print(pl151)
dev.off()

# 16.##Perumal (2022)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 25.5*(FFM/49.6)^0.75
# Vd  (L)    = 90.1*(FFM/49.6)
# Ka  (1/h)  = 1.64
# Tlag (h)   = 0.38
# F.         = 1-0.253*PXR+0.193*MOX
# BSV (CV%): CL = 9.1%
# BSV (CV%): F  = 43.1%
# BSV (CV%): Ka = 89.2%
# BSV (CV%): Tlag = 88.9%
# prop.err= 36.3%
# add.err= 0.008

# Mrgsolve code----
set.seed(92001)

cod16 <- '
$PARAM CL=25.5, VC=90.1, KA=1.64, TLAG=0.38, BIO=1
$CMT DEPOT CENT
$MAIN
ALAG_DEPOT = TLAG;
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT;
$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA TLAG
'
mod16<- mcode("optim", cod16, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff16 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]],
                   KA = p[["KA"]],
                   TLAG = p[["TLAG"]])
  out <- mrgsim_q(mod16, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg16 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * (a[3]/49.6)^0.75 * exp(b[1]),
    VC = bpop[2] * (a[3]/49.6),
    KA = bpop[3] * exp(b[2]),
    TLAG = bpop[4] * exp(b[3]),
    DOSE = a[1],
    TAU = a[2],
    FFM = a[3],
    AID = a[4]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model16 <- create.poped.database(
  ff_fun =ff16,
  fg_fun =fg16,
  fError_fun=feps1,
  bpop=c(CL=25.5,
         VC=90.1,
         KA=1.642,
         TLAG=0.38), 
  notfixed_bpop = c(1,1,1,0),
  d=c(CL=0.091^2,KA=0.892^2,TLAG=0.889^2),
  sigma=c(prop=0,add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,FFM=56,AID=14))

pl16 <- plot_model_prediction(poped_db_model16,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model16: Perumal (2022)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl16

jpeg(filename = paste0(output_dir,"/16.Perumal2022_Adult.jpg"), 
     width=1000, height=1000, res=300)
print(pl16)
dev.off()

# 17.##Savic (2015)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 5.71
# Vd  (L)    = 24.9
# Q   (L/h)  = 9.46
# Vp  (L)    = 12.4
# Ka  (1/h)  = 0.644
# F          = 0.6
# BSV (CV%): CL = 34%
# BSV (CV%): Vd = 60%
# BSV correlation: CL_V = 0.53
# prop.err= 28%
# add.err= 2.18

# Mrgsolve code----
set.seed(92001)

cod17 <- '
$PARAM CL=5.71, VC=24.9, KA=0.644, Q=9.46, VP=12.4, BIO=0.6
$CMT DEPOT CENT PERI
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = BIO*KA*DEPOT - (CL/VC)*CENT - (Q/VC)*CENT + (Q/VP)*PERI;
dxdt_PERI = (Q/VC)*CENT - (Q/VP)*PERI;
$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA Q VP
'
mod17<- mcode("optim", cod17, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff17 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]],
                   Q  = p[["Q"]], 
                   VP = p[["VP"]],
                   KA = p[["KA"]])
  out <- mrgsim_q(mod17, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg17 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * exp(b[1]),
    VC = bpop[2] * exp(b[2]),
    Q  = bpop[3],
    VP = bpop[4],
    KA = bpop[5],
    DOSE = a[1],
    TAU = a[2],
    AID = a[3]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model17 <- create.poped.database(
  ff_fun =ff17,
  fg_fun =fg17,
  fError_fun=feps1,
  bpop=c(CL=5.71,
         VC=24.9,
         Q=9.46,
         VP=12.4,
         KA=0.644), 
  notfixed_bpop = c(1,1,1,1,1),
  d=c(CL=0.34^2,VC=0.6^2),
  covd = 0.10812,
  sigma=c(prop=0,add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,AID=17))

pl17 <- plot_model_prediction(poped_db_model17,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model17: Savic (2015)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl17

jpeg(filename = paste0(output_dir,"/17.Savic2015_Adult.jpg"), 
     width=1000, height=1000, res=300)
print(pl17)
dev.off()

# 18.##Sekaggya (2019)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 23.9*(BW/70)^0.75*(AGE/33)^0.517
# Vd  (L)    = 44.6*(BW/70)
# Ka  (1/h)  = 0.236
# F          = 1*0.517^Child
# BSV (CV%): CL = 46.6%
# BSV (CV%): Vd = 87.4%
# prop.err= 48%

# Mrgsolve code----
set.seed(92001)

cod18 <- '
$PARAM CL=23.9, VC=44.6, KA=0.236, BIO=1
$CMT DEPOT CENT
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = BIO*KA*DEPOT - (CL/VC)*CENT;
$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA BIO
'
mod18 <- mcode("optim", cod18, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff18 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]],
                   KA = p[["KA"]],
                   BIO = p[["BIO"]])
  out <- mrgsim_q(mod18, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg18 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * (a[3]/70)^0.75 * (a[4]/33)^0.517 * exp(b[1]),
    VC = bpop[2] * (a[3]/70) * exp(b[2]),
    KA = bpop[3],
    BIO = bpop[4] * 0.517^a[5],
    DOSE = a[1],
    TAU = a[2],
    BW = a[3],
    AGE = a[4],
    POP = a[5],
    AID = a[6]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model18 <- create.poped.database(
  ff_fun =ff18,
  fg_fun =fg18,
  fError_fun=feps1,
  bpop=c(CL=23.9,
         VC=44.6,
         KA=0.236,
         BIO=1), 
  notfixed_bpop = c(1,1,1,0),
  d=c(CL=0.466^2,VC=0.874^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,BW=70,AGE=33,POP=0,AID=18))

pl18 <- plot_model_prediction(poped_db_model18,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model18: Schipani (2016)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl18

jpeg(filename = paste0(output_dir,"/18.Schipani2016_Adult.jpg"), 
     width=1000, height=1000, res=300)
print(pl18)
dev.off()

# 19.##Sekaggya (2019)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = (12.2+20^week) * (FFM/43)^0.75
# Vd  (L)    = 58*(FFM/43)
# Ka  (1/h)  = 1.99
# Tlag (h)   = 0.83
# F.         = 1 FIX
# BSV (CV%): CL = 46.6%
# BSV (CV%): Vd = 87.4%
# prop.err= 48%

# Mrgsolve code----
set.seed(92001)

cod19 <- '
$PARAM CL=250, VC=58, KA=1.99, TLAG=0.83
$CMT DEPOT CENT
$MAIN
ALAG_DEPOT = TLAG;
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT;
$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA
'
mod19<- mcode("optim", cod19, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff19 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]],
                   KA = p[["KA"]])
  out <- mrgsim_q(mod19, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg19 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = (bpop[1]+20^a[4]) * (a[3]/43)^0.75 * exp(b[1]),
    VC = bpop[2] * (a[3]/43) * exp(b[2]),
    KA = bpop[3],
    DOSE = a[1],
    TAU = a[2],
    FFM = a[3],
    WEEK = a[4],
    AID = a[5]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model19 <- create.poped.database(
  ff_fun =ff19,
  fg_fun =fg19,
  fError_fun=feps1,
  bpop=c(CL=12.2,
         VC=58,
         KA=1.99), 
  notfixed_bpop = c(1,1,1),
  d=c(CL=0.466^2,VC=0.874^2),
  sigma=c(prop=0,add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,FFM=56,WEEK=0,AID=19))

pl19 <- plot_model_prediction(poped_db_model19,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model19: Sekaggya (2019) 2wk",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl19

jpeg(filename = paste0(output_dir,"/19.Sekaggya2019_Adult_week2.jpg"), 
     width=1000, height=1000, res=300)
print(pl19)
dev.off()

poped_db_model191 <- create.poped.database(
  ff_fun =ff19,
  fg_fun =fg19,
  fError_fun=feps1,
  bpop=c(CL=12.2,
         VC=58,
         KA=1.99), 
  notfixed_bpop = c(1,1,1),
  d=c(CL=0.466^2,VC=0.874^2),
  sigma=c(prop=0,add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,FFM=56,WEEK=1,AID=19))

pl191 <- plot_model_prediction(poped_db_model191,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model19: Sekaggya (2019) 8wk",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl191

jpeg(filename = paste0(output_dir,"/19.Sekaggya2019_Adult_week8.jpg"), 
     width=1000, height=1000, res=300)
print(pl191)
dev.off()

# 20.##Seng (2015)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 10.3*(BW/70)^0.75
# Vd  (L)    = 30.9*(BW/70)
# Ka  (1/h)  = 2.15
# F   (%)    = 1
# NN         = 2 FIX
# BSV (CV%): CL = 30.13%
# BOV (CV%): Ka = 57.71%
# BOV (CV%): F  = 51.09%
# prop.err= 26.4%

# Mrgsolve code----
set.seed(92001)

cod20 <- '
$PARAM CL=10.3, VC=30.9, KA=2.15, NN=2, BIO=1,
$CMT DEPOT CENT
$GLOBAL 
int NDOSE = 0;
double dosetime[24];
double dose[450];

$MAIN
if(NEWIND < 2) NDOSE = 0; 

if(self.amt > 0 && self.cmt==1) {
 NDOSE = NDOSE + 1; 
 dosetime[NDOSE] = self.time;
 dose[NDOSE] = self.amt;
 }

F_DEPOT = 0;
double KTR = KA;
double NFAC = exp(lgamma(NN+1));
double KINPT = BIO * pow(KTR,(NN+1)) / NFAC; 

$ODE
double INPT = 0;
int i = 0;
while(i <= NDOSE) {
  double IPT = 0;
  if(SOLVERTIME >= dosetime[i]) {
    double delta = SOLVERTIME - dosetime[i];
    IPT = dose[i] * pow(delta, NN) * exp(-KTR * delta);  
  }
  INPT = INPT + IPT;
  ++i;
 }
dxdt_DEPOT = KINPT * INPT - KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT;

$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA
'

mod20 <- mcode("transit", cod20, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff20 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]], 
                   KA = p[["KA"]])
  out <- mrgsim_q(mod20, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg20 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * (a[3]/70)^(0.75) * exp(b[1]),
    VC = bpop[2] * (a[3]/70),
    KA = bpop[3] * exp(b[2]),
    DOSE = a[1],
    TAU = a[2],
    BW = a[3],
    AID = a[4]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model20 <- create.poped.database(
  ff_fun =ff20,
  fg_fun =fg20,
  fError_fun=feps1,
  bpop=c(CL=10.3,
         VC=30.9,
         KA=2.15), 
  notfixed_bpop = c(1,1,1),
  d=c(CL=0.3013^2,KA=0.5^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,BW=70,AID=20))

pl20 <- plot_model_prediction(poped_db_model20,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model20: Seng (2015)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl20

jpeg(filename = paste0(output_dir,"/20.Seng2015_Adult.jpg"), 
     width=1000, height=1000, res=300)
print(pl20)
dev.off()

# 21.##Sloan (2017)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 19.6*1.2^SEX
# Vd  (L)    = 23.6
# Ka  (1/h)  = 0.277 FIX
# F   (%)    = 1
# MTT (h)    = 0.326 FIX
# NN         = 1.5 FIX
# BSV (CV%): CL = 27.57%
# BOV (CV%): Vd = 63.01%
# BOV (CV%): MTT= 26.57%
# prop.err= 22%

# Mrgsolve code----
set.seed(92001)

cod21 <- '
$PARAM CL=19.6, VC=23.6, KA=0.277, MTT=0.326, NN=1.5, BIO=1,
$CMT DEPOT CENT
$GLOBAL 
int NDOSE = 0;
double dosetime[24];
double dose[450];

$MAIN
if(NEWIND < 2) NDOSE = 0; 

if(self.amt > 0 && self.cmt==1) {
 NDOSE = NDOSE + 1; 
 dosetime[NDOSE] = self.time;
 dose[NDOSE] = self.amt;
 }

F_DEPOT = 0;
double KTR = (NN+1)/MTT;
double NFAC = exp(lgamma(NN+1));
double KINPT = BIO * pow(KTR,(NN+1)) / NFAC; 

$ODE
double INPT = 0;
int i = 0;
while(i <= NDOSE) {
  double IPT = 0;
  if(SOLVERTIME >= dosetime[i]) {
    double delta = SOLVERTIME - dosetime[i];
    IPT = dose[i] * pow(delta, NN) * exp(-KTR * delta);  
  }
  INPT = INPT + IPT;
  ++i;
 }
dxdt_DEPOT = KINPT * INPT - KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT;

$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA MTT
'

mod21 <- mcode("transit", cod21, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff21 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]], 
                   KA = p[["KA"]], 
                   MTT = p[["MTT"]])
  out <- mrgsim_q(mod21, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg21 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * 1.2^a[3] * exp(b[1]),
    VC = bpop[2] * exp(b[2]),
    KA = bpop[3],
    MTT = bpop[4] * exp(b[3]),
    DOSE = a[1],
    TAU = a[2],
    SEX = a[3],
    AID = a[4]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model21 <- create.poped.database(
  ff_fun =ff21,
  fg_fun =fg21,
  fError_fun=feps1,
  bpop=c(CL=19.6,
         VC=23.6,
         KA=0.277,
         MTT=0.326), 
  notfixed_bpop = c(1,1,0,0),
  d=c(CL=0.2757^2,VC=0.6301^2,MTT=0.2657^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,SEX=1,AID=21))

pl21 <- plot_model_prediction(poped_db_model21,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model21: Sloan (2017) Male",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl21

jpeg(filename = paste0(output_dir,"/21.Sloan2017_Adult_M.jpg"), 
     width=1000, height=1000, res=300)
print(pl21)
dev.off()

poped_db_model211 <- create.poped.database(
  ff_fun =ff21,
  fg_fun =fg21,
  fError_fun=feps1,
  bpop=c(CL=19.6,
         VC=23.6,
         KA=0.277,
         MTT=0.326), 
  notfixed_bpop = c(1,1,0,0),
  d=c(CL=0.2757^2,VC=0.6301^2,MTT=0.2657^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,SEX=0,AID=21))

pl211 <- plot_model_prediction(poped_db_model211,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model21: Sloan (2017) Female",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl211

jpeg(filename = paste0(output_dir,"/21.Sloan2017_Adult_F.jpg"), 
     width=1000, height=1000, res=300)
print(pl211)
dev.off()

# 22.##Soedarsono (2023)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 7.85*(AGE/40)^(-0.55)*(1+SLCO*0.24)*(FFM/44)^0.75
# Vd  (L)    = 30.9*(FFM/44)
# Ka  (1/h)  = 0.37
# BSV (CV%): CL = 78.1%
# add.err= 0.29

# Mrgsolve code----
set.seed(92001)

cod22 <- '
$PARAM CL=7.85, VC=30.9, KA=0.37
$CMT DEPOT CENT
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT;
$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA
'
mod22 <- mcode("optim", cod22, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff22 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]],
                   KA = p[["KA"]])
  out <- mrgsim_q(mod22, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg22 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1]*(a[3]/40)^(-0.55)*(1+a[4]*0.24)*(a[5]/44)^0.75 * exp(b[1]),
    VC = bpop[2]*(a[5]/44),
    KA = bpop[3],
    DOSE = a[1],
    TAU = a[2],
    AGE = a[3],
    SLCO= a[4],
    FFM = a[5],
    AID = a[6]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model22 <- create.poped.database(
  ff_fun =ff22,
  fg_fun =fg22,
  fError_fun=feps1,
  bpop=c(CL=7.85,
         VC=30.9,
         KA=0.37), 
  notfixed_bpop = c(1,1,1),
  d=c(CL=0.781^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,AGE=40,SLCO=0,FFM=56,AID=22))

pl22 <- plot_model_prediction(poped_db_model22,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model22: Soedarsono (2023) NS",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl22

jpeg(filename = paste0(output_dir,"/22.Soedarsono2023_Adult_NS.jpg"), 
     width=1000, height=1000, res=300)
print(pl22)
dev.off()

poped_db_model221 <- create.poped.database(
  ff_fun =ff22,
  fg_fun =fg22,
  fError_fun=feps1,
  bpop=c(CL=7.85,
         VC=30.9,
         KA=0.37), 
  notfixed_bpop = c(1,1,1),
  d=c(CL=0.781^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,AGE=40,SLCO=1,FFM=56,AID=22))

pl221 <- plot_model_prediction(poped_db_model221,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model22: Soedarsono (2023) S",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl221

jpeg(filename = paste0(output_dir,"/22.Soedarsono2023_Adult_S.jpg"), 
     width=1000, height=1000, res=300)
print(pl221)
dev.off()

# 24.##Wilkins (2008)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 19.2
# Vd  (L)    = 53.2
# Ka  (1/h)  = 1.15
# F   (%)    = 1
# MTT (h)    = 0.424 
# NN         = 7.13 
# BSV (CV%): CL = 22.5%
# BOV (CV%): MTT= 67.9%
# prop.err= 22.2%
# add.err= 0.0923

# Mrgsolve code----
set.seed(92001)

cod24 <- '
$PARAM CL=19.2, VC=53.2, KA=1.15, MTT=0.424, NN=7.13, BIO=1,
$CMT DEPOT CENT
$GLOBAL 
int NDOSE = 0;
double dosetime[24];
double dose[450];

$MAIN
if(NEWIND < 2) NDOSE = 0; 

if(self.amt > 0 && self.cmt==1) {
 NDOSE = NDOSE + 1; 
 dosetime[NDOSE] = self.time;
 dose[NDOSE] = self.amt;
 }

F_DEPOT = 0;
double KTR = (NN+1)/MTT;
double NFAC = exp(lgamma(NN+1));
double KINPT = BIO * pow(KTR,(NN+1)) / NFAC; 

$ODE
double INPT = 0;
int i = 0;
while(i <= NDOSE) {
  double IPT = 0;
  if(SOLVERTIME >= dosetime[i]) {
    double delta = SOLVERTIME - dosetime[i];
    IPT = dose[i] * pow(delta, NN) * exp(-KTR * delta);  
  }
  INPT = INPT + IPT;
  ++i;
 }
dxdt_DEPOT = KINPT * INPT - KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT;

$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA MTT NN
'

mod24 <- mcode("transit", cod24, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff24 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]], 
                   KA = p[["KA"]], 
                   NN = p[["NN"]], 
                   MTT = p[["MTT"]])
  out <- mrgsim_q(mod24, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg24 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * exp(b[1]),
    VC = bpop[2] * exp(b[2]),
    KA = bpop[3] * exp(b[3]),
    NN = bpop[4] * exp(b[4]),
    MTT = bpop[5] * exp(b[5]),
    DOSE = a[1],
    TAU = a[2],
    AID = a[3]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model24 <- create.poped.database(
  ff_fun =ff24,
  fg_fun =fg24,
  fError_fun=feps1,
  bpop=c(CL=19.2,
         VC=53.2,
         KA=1.15,
         NN=7.13,
         MTT=0.424), 
  notfixed_bpop = c(1,1,1,0,0),
  d=c(CL=0.528^2,VC=0.434^2,KA=0.663^2,NN=1.56^2,MTT=0.601^2),
  covd = 0.217,
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,AID=24))

pl24 <- plot_model_prediction(poped_db_model24,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model24: Wilkins (2008)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl24

jpeg(filename = paste0(output_dir,"/24.Wilkins2008_Adult.jpg"), 
     width=1000, height=1000, res=300)
print(pl24)
dev.off()

# 25.##Naidoo (2019)# -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 22.8*(FFM/47)^0.75
# Vd  (L)    = 77.4*(FFM/47)
# Ka  (1/h)  = 1.57
# F   (100%) = 1
# MTT (h)    = 0.53
# NN         = 34.6 
# BSV (CV%): CL = 23.45%
# prop.err= 37.4%
# add.err= 0.008 FIX

# Mrgsolve code----
set.seed(92001)

cod25 <- '
$PARAM CL=22.8, VC=77.4, KA=1.57, MTT=0.53, NN=34.6, BIO=1,
$CMT DEPOT CENT
$GLOBAL 
int NDOSE = 0;
double dosetime[24];
double dose[450];

$MAIN
if(NEWIND < 2) NDOSE = 0; 

if(self.amt > 0 && self.cmt==1) {
 NDOSE = NDOSE + 1; 
 dosetime[NDOSE] = self.time;
 dose[NDOSE] = self.amt;
 }

F_DEPOT = 0;
double KTR = (NN+1)/MTT;
double NFAC = exp(lgamma(NN+1));
double KINPT = BIO * pow(KTR,(NN+1)) / NFAC; 

$ODE
double INPT = 0;
int i = 0;
while(i <= NDOSE) {
  double IPT = 0;
  if(SOLVERTIME >= dosetime[i]) {
    double delta = SOLVERTIME - dosetime[i];
    IPT = dose[i] * pow(delta, NN) * exp(-KTR * delta);  
  }
  INPT = INPT + IPT;
  ++i;
 }
dxdt_DEPOT = KINPT * INPT - KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT;

$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA MTT NN
'

mod25 <- mcode("transit", cod25, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff25 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]], 
                   KA = p[["KA"]], 
                   NN = p[["NN"]], 
                   MTT = p[["MTT"]])
  out <- mrgsim_q(mod25, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg25 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * (a[3]/47)^0.75 * exp(b[1]),
    VC = bpop[2] * (a[3]/47),
    KA = bpop[3],
    NN = bpop[4],
    MTT = bpop[5],
    DOSE = a[1],
    TAU = a[2],
    FFM = a[3],
    AID = a[4]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model25 <- create.poped.database(
  ff_fun =ff25,
  fg_fun =fg25,
  fError_fun=feps1,
  bpop=c(CL=22.8,
         VC=77.4,
         KA=1.57,
         NN=34.6,
         MTT=0.53), 
  notfixed_bpop = c(1,1,1,0,0),
  d=c(CL=0.2345^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=600,TAU=24,FFM=56,AID=24))

pl25 <- plot_model_prediction(poped_db_model25,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model25: Naidoo (2019)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl25

jpeg(filename = paste0(output_dir,"/25.Naidoo2019_Adult.jpg"), 
     width=1000, height=1000, res=300)
print(pl25)
dev.off()

# 26.##Aruldhas (2019)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 8.11*(BW/19.5)^0.75
# Vd  (L)    = 44.7*(BW/19.5)
# Ka  (1/h)  = KTR
# F   (%)    = 1
# MTT (h)    = 0.932
# NN         = 9 FIX
# BSV (CV%): Vd = 42%
# BSV (CV%): F  = 68%
# BSV (CV%): MTT= 52.2%
# add.err = 0.0967 FIX mg/L

# Mrgsolve code----
set.seed(92001)

cod26 <- '
$PARAM CL=8.11, VC=44.7, MTT=0.932, NN=9, BIO=1,
$CMT DEPOT CENT
$GLOBAL 
int NDOSE = 0;
double dosetime[0];
double dose[150];
$MAIN
if(NEWIND < 2) NDOSE = 0; 

if(self.amt > 0 && self.cmt==1) {
 NDOSE = NDOSE + 1; 
 dosetime[NDOSE] = self.time;
 dose[NDOSE] = self.amt;
 }

F_DEPOT = 0;
double KTR = (NN+1)/MTT;
double NFAC = exp(lgamma(NN+1));
double KINPT = BIO * pow(KTR,(NN+1)) / NFAC; 

$ODE
double INPT = 0;
int i = 0;
while(i <= NDOSE) {
  double IPT = 0;
  if(SOLVERTIME >= dosetime[i]) {
    double delta = SOLVERTIME - dosetime[i];
    IPT = dose[i] * pow(delta, NN) * exp(-KTR * delta);  
  }
  INPT = INPT + IPT;
  ++i;
 }
dxdt_DEPOT = KINPT * INPT - KTR*DEPOT;
dxdt_CENT = KTR*DEPOT - (CL/VC)*CENT;

$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC BIO MTT
'

mod26 <- mcode("transit", cod26, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff26 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]], 
                   BIO = p[["BIO"]], 
                   MTT = p[["MTT"]])
  out <- mrgsim_q(mod26, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg26 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * (a[3]/19.5)^(0.75),
    VC = bpop[2] * (a[3]/19.5) * exp(b[1]),
    BIO = bpop[3] * exp(b[2]),
    MTT = bpop[4] * exp(b[3]),
    DOSE = a[1],
    TAU = a[2],
    BW = a[3],
    AID = a[4]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model26 <- create.poped.database(
  ff_fun =ff26,
  fg_fun =fg26,
  fError_fun=feps1,
  bpop=c(CL=8.11,
         VC=44.7,
         BIO=1,
         MTT=0.932), 
  notfixed_bpop = c(1,1,0,0),
  d=c(VC=0.42^2,BIO=0.68^2,MTT=0.522^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=150,TAU=24,BW=14,AID=26))

pl26 <- plot_model_prediction(poped_db_model26,
                             PI=T,
                             sample.times = F,
                             PRED = T) +
  labs(title = "Model26: Aruldhas (2019)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl26

jpeg(filename = paste0(output_dir,"/26.Aruldhas2019_Children.jpg"), 
     width=1000, height=1000, res=300)
print(pl26)
dev.off()

# 27.##Denti (2022)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 54.5*(FFM/9)^0.75*(PMAy^3.22/(PMAy^3.22+1.04^3.22))
# Vd  (L)    = 12.3*(FFM/9)
# Ka  (1/h)  = 1.82
# F   (%)    = 0.655+((1-0.655)/2.72)*Age
# MTT (h)    = 0.589
# NN         = 9.7
# BSV (CV%): CL = 41.8%
# BOV (CV%): Ka = 111%
# BOV (CV%): F  = 45.1%
# BOV (CV%): MTT= 58.8%
# prop.err = 13.7 %
# add.err = 0.0967 FIX mg/L

# RxODE code----
set.seed(92001)

mod27 <- rxode2({
  MTT = 0.589;
  NN  = 9.7;
  K   = CL/VC;
  KTR = (NN+1)/MTT;
 
  d/dt(DEPOT) = exp(log(BIO*podo(DEPOT))+log(KTR)+NN*log(KTR*tad(DEPOT))-KTR*tad(DEPOT)-lgammafn(NN+1))-KA*DEPOT;
  d/dt(CENT) = KA*DEPOT-K*CENT;
  
  CP=CENT/VC;
})

## -- Model structure definition function----
ff27 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)
  et(0,amt=p[["DOSE"]], ii=p[["TAU"]], until=max(times_xt)) %>%
    et(times_xt) -> data
  
  out <- rxSolve(mod27, p, data, atol=1e-8, rtol=1e-8, maxsteps=50000,
                 returnType="data.frame")
  
  y <-  out$CP[match(times_xt,out$time)]
  
  return(list(y=matrix(y,ncol=1),poped.db=poped.db))
  
}

## -- parameter definition function----
fg27 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * (a[3]/9)^(0.75)*(a[4]^3.22/(a[4]^3.22+1.04^3.22)) * exp(b[1]),
    VC = bpop[2] * (a[3]/9),
    KA = bpop[3],
    BIO = bpop[4]+((1-0.655)/2.72)*a[5],
    DOSE = a[1],
    TAU = a[2],
    FFM = a[3],
    PMA = a[4],
    AGE = a[5],
    AID = a[6]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model27 <- create.poped.database(
  ff_fun =ff27,
  fg_fun =fg27,
  fError_fun=feps1,
  bpop=c(CL=54.5,
         VC=12.3,
         KA=1.82,
         BIO=0.655), 
  notfixed_bpop = c(1,1,1,0),
  d=c(CL=0.418^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=150,TAU=24,FFM=11,PMA=2.75,AGE=2,AID=27))

pl27 <- plot_model_prediction(poped_db_model27,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model27: Denti (2022)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl27

jpeg(filename = paste0(output_dir,"/27.Denti2022_Children.jpg"), 
     width=1000, height=1000, res=300)
print(pl27)
dev.off()

# 28.##Horita (2018)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 5.96
# Vd  (L)    = 0.7*BW
# Ka  (1/h)  = 1.24
# D0  (h)    = 0.342
# BSV (CV%): CL =38.5%
# BSV (CV%): Vd = 26.8%
# BSV (CV%): Ka = 110.5%
# BSV (CV%): D0 = 131.1%
# add.err = 2.256 FIX mg/L

# Mrgsolve code----
cod28 <- '
$PARAM CL=5.96, VC=49, KA=1.24, DUR=0.62, TLAG=0.24
$CMT DEPOT CENT
$MAIN
D_CENT = DUR;
ALAG_DEPOT = TLAG;
$ODE
dxdt_DEPOT = -KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT;
$TABLE double CP  = CENT/VC;
$CAPTURE CP CL VC KA
'

mod28<- mcode("optim", cod28, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff28 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]],
                   KA = p[["KA"]])
  out <- mrgsim_q(mod28, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg28 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * exp(b[1]),
    VC = bpop[2] * a[3] * exp(b[2]),
    KA = bpop[3] * exp(b[3]),
    DOSE = a[1],
    TAU = a[2],
    BW = a[3],
    AID = a[4]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model28 <- create.poped.database(
  ff_fun =ff28,
  fg_fun =fg28,
  fError_fun=feps1,
  bpop=c(CL=5.96,
         VC=0.7,
         KA=1.24), 
  notfixed_bpop = c(1,1,1),
  d=c(CL=0.385^2,VC=0.268^2,KA=1.105^2),
  sigma=c(prop=0,add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=450,TAU=24,BW=70,AID=14))

pl28 <- plot_model_prediction(poped_db_model28,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model28: Horita (2018)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl28

jpeg(filename = paste0(output_dir,"/28.Horita2018_Children.jpg"), 
     width=1000, height=1000, res=300)
print(pl28)
dev.off()

# 29.##Panjasawatwong (2020)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 3.22*(BW/10.9)^0.75*(PMAm^1.38/(PMAm^1.28+6.81^1.38))
# Vd  (L)    = 12.3*(BW/10.9)
# Ka  (1/h)  = 1.24
# MTT (h)    = 1.25
# NN         = 2 FIX
# BSV (CV%): CL = 19.4%
# BSV (CV%): Vd = 23%
# BSV (CV%): MTT= 85%
# add.err = 0.513 mg/L

# Mrgsolve code----
set.seed(92001)

cod29 <- '
$PARAM CL=3.22, VC=12.3, KA=1.24, MTT=1.25, NN=2, BIO=1
$CMT DEPOT CENT
$GLOBAL 
int NDOSE = 0;
double dosetime[150];
double dose[150];
$MAIN
if(NEWIND < 2) NDOSE = 0; 

if(self.amt > 0 && self.cmt==1) {
 NDOSE = NDOSE + 1; 
 dosetime[NDOSE] = self.time;
 dose[NDOSE] = self.amt;
 }

F_DEPOT = 0;
double KTR = (NN+1)/MTT;
double NFAC = exp(lgamma(NN+1));
double KINPT = BIO * pow(KTR,(NN+1)) / NFAC; 

$ODE
double INPT = 0;
int i = 0;
while(i <= NDOSE) {
  double IPT = 0;
  if(SOLVERTIME >= dosetime[i]) {
    double delta = SOLVERTIME - dosetime[i];
    IPT = dose[i] * pow(delta, NN) * exp(-KTR * delta);  
  }
  INPT = INPT + IPT;
  ++i;
 }
dxdt_DEPOT = KINPT * INPT - KA*DEPOT;
dxdt_CENT = KA*DEPOT - (CL/VC)*CENT;

$TABLE double CP  = CENT/VC;
$CAPTURE CP
'

mod29 <- mcode("transit", cod29, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff29 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]], 
                   KA = p[["KA"]], 
                   MTT = p[["MTT"]])
  out <- mrgsim_q(mod29, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg29 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * (a[3]/10.9)^(0.75)*(a[4]^1.38/(a[4]^1.38+6.81^1.38)) * exp(b[1]),
    VC = bpop[2] * (a[3]/10.9) * exp(b[2]),
    KA = bpop[3],
    MTT = bpop[4] * exp(b[3]),
    DOSE = a[1],
    TAU = a[2],
    BW = a[3],
    PMA = a[4],
    AID = a[5]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model29 <- create.poped.database(
  ff_fun =ff29,
  fg_fun =fg29,
  fError_fun=feps1,
  bpop=c(CL=3.22,
         VC=12.3,
         KA=1.24,
         MTT=1.25), 
  notfixed_bpop = c(1,1,1,0),
  d=c(CL=0.194^2,VC=0.23^2,MTT=0.85^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=150,TAU=24,FFM=14,PMA=33,AID=29))

pl29 <- plot_model_prediction(poped_db_model29,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model29: Panjasawatwong (2020)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl29

jpeg(filename = paste0(output_dir,"/29.Panjasawatwong2020_Children.jpg"), 
     width=1000, height=1000, res=300)
print(pl29)
dev.off()

# 30.##Zvada (2014)## -----------------------------------------------------
# Define model ------------------------------------------------------------
# CL  (L/h)  = 8.15*(BW/12.5)^0.75*(1/(1+(PMAw/58.2)^-2.21))
# Vd  (L)    = 16.2*(BW/12.5)
# Ka  (1/h)  = KTR
# MTT (h)    = 1.04
# NN         = 8.04
# BSV (CV%): CL = 57.1%
# BSV (CV%): Vd = 65.88%
# prop.err = 23.4%
# add.err = 0.122 mg/L

# Mrgsolve code----
set.seed(92001)

cod30 <- '
$PARAM CL=8.15, VC=16.2, MTT=1.04, NN=8.04, BIO=1
$CMT DEPOT CENT
$GLOBAL 
int NDOSE = 0;
double dosetime[150];
double dose[150];
$MAIN
if(NEWIND < 2) NDOSE = 0; 

if(self.amt > 0 && self.cmt==1) {
 NDOSE = NDOSE + 1; 
 dosetime[NDOSE] = self.time;
 dose[NDOSE] = self.amt;
 }

F_DEPOT = 0;
double KTR = (NN+1)/MTT;
double NFAC = exp(lgamma(NN+1));
double KINPT = BIO * pow(KTR,(NN+1)) / NFAC; 

$ODE
double INPT = 0;
int i = 0;
while(i <= NDOSE) {
  double IPT = 0;
  if(SOLVERTIME >= dosetime[i]) {
    double delta = SOLVERTIME - dosetime[i];
    IPT = dose[i] * pow(delta, NN) * exp(-KTR * delta);  
  }
  INPT = INPT + IPT;
  ++i;
 }
dxdt_DEPOT = KINPT * INPT - KTR*DEPOT;
dxdt_CENT = KTR*DEPOT - (CL/VC)*CENT;

$TABLE double CP  = CENT/VC;
$CAPTURE CP
'

mod30 <- mcode("transit", cod30, atol=1e-8, rtol=1e-8, maxsteps=50000)

## -- Model structure definition function----
ff30 <- function(model_switch, xt, p, poped.db){
  times_xt <- drop(xt)  
  dose_times <- seq(from=0,to=max(times_xt),by=p[["TAU"]])
  time <- sort(unique(c(0,times_xt,dose_times)))
  is.dose <- time %in% dose_times
  data <- 
    tibble::tibble(ID = 1,
                   time = time,
                   amt = ifelse(is.dose,p[["DOSE"]], 0), 
                   cmt = ifelse(is.dose, 1, 0), 
                   evid = cmt,
                   CL = p[["CL"]], 
                   VC = p[["VC"]], 
                   MTT = p[["MTT"]])
  out <- mrgsim_q(mod30, data=data)
  y <-  out$CP
  y <- y[match(times_xt,out$time)]
  
  return(list(y=y,poped.db=poped.db))
}

## -- parameter definition function----
fg30 <- function(x, a, bpop, b, bocc){
  parameters = c(
    CL = bpop[1] * (a[3]/12.5)^(0.75)*(1/(1+(a[4]/58.2)^-2.21)) * exp(b[1]),
    VC = bpop[2] * (a[3]/12.5) * exp(b[2]),
    MTT = bpop[3],
    DOSE = a[1],
    TAU = a[2],
    BW = a[3],
    PMA = a[4],
    AID = a[5]
  )
  return(parameters)
}

## -- Residual unexplained variablity (RUV) function----
## prop.err & add.err
feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
  y <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
  y = y * (1+epsi[,1]) + epsi[,2]
  return(list(y=y,poped.db=poped.db)) 
}

## Create PopED databases----
poped_db_model30 <- create.poped.database(
  ff_fun =ff30,
  fg_fun =fg30,
  fError_fun=feps1,
  bpop=c(CL=8.15,
         VC=16.2,
         MTT=1.04), 
  notfixed_bpop = c(1,1,0),
  d=c(CL=0.571^2,VC=0.6588^2),
  sigma=c(prop=0, add=0),
  m=1,
  groupsize=1000,
  xt=seq(168,192,0.1),
  a=cbind(DOSE=150,TAU=24,BW=14,PMA=140,AID=30))

pl30 <- plot_model_prediction(poped_db_model30,
                              PI=T,
                              sample.times = F,
                              PRED = T) +
  labs(title = "Model30: Zvada (2014)",
       x = " ",
       y = " ") +
  scale_y_continuous(limits = c(-5,35))+
  annotate("rect", xmin = 168, xmax = 192,  ymin = 8, ymax = 24, alpha = 0.1, fill = "blue")+
  theme_bw(base_size = 10)+
  theme (legend.position = "none",
         plot.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())
pl30

jpeg(filename = paste0(output_dir,"/30.Zvada2014_Children.jpg"), 
     width=1000, height=1000, res=300)
print(pl30)
dev.off()

