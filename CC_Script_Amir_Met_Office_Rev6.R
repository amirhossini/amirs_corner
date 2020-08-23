## Install Packages 
# install.packages('TSstudio')
# install.packages('forecast')
# install.packages('roll')

## Import Libraries -------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(dplyr)
library(mice)
library(TSstudio)
library(forecast)
library(seasonal)
library(tseries)
library(roll)

## Read in files ----------------------------------------------------------
infl_metoffice <- '../../01-Data/CRU.csv'
df_mtof <- read.csv(infl_metoffice, header=TRUE)

## Clean the Data & Create Absolute Temperatures -------------------------
df_mtof <- na.omit(df_mtof)# Remove the coverage info from Met Office Data  
df_mtof <- df_mtof[,-c(14)]
df_mtof <- df_mtof[-nrow(df_mtof),]
df_mtof <- gather(df_mtof, Month, DTemp, Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec)
df_mtof <- df_mtof[order(df_mtof$Year),]
rownames(df_mtof) <- 1:nrow(df_mtof)
df_mtof$ATemp <- df_mtof$DTemp + 14 + 273.15
qplot(df_mtof$ATemp)

## Time Series: Set-up & Visualization -----------------------------------
mtof_total = ts(df_mtof$ATemp,start=c(1850,1), frequency=12)
autoplot(mtof_total)+
  labs(title='Met Office Monthly Temperature', 
       y='Temperature (K)', x=NULL)

ts_heatmap(mtof_total)

ggseasonplot(mtof_total, year.labels=TRUE, continuous=TRUE, polar = T)+
  labs(title = "Seasonal plots Met Office Temperature Data ",
       y="Temperature (K)",
       x=NULL)

## Auto-correlation Function ---------------------------------------------
ggAcf(mtof_total, lag=120)
ggPacf(mtof_total, lag=120)

## Decompostion ----------------------------------------------------------
### STL
frequency(mtof_total)
SWIN <- as.integer(frequency(mtof_total)/2)*2+1 
#TWIN <- 30
mtof_fitStl <-stl(mtof_total ,s.window=SWIN, robust=TRUE)
autoplot(mtof_fitStl)
mtof_fitStl_ts <- mtof_fitStl$time.series
write.csv(mtof_fitStl$time.series, 'STL_decomp.csv') 
roll_sd(mtof_fitStl_ts,120)


## Total Data Split ------------------------------------------------------
mtof.split.size <- as.integer(length(mtof_total)/10)
mtof_split <- ts_split(mtof_total, sample.out = mtof.split.size)
mtof.train <- mtof_split$train
mtof.test <- mtof_split$test

## Stationarity & Differencing -------------------------------------------
### Stationarity *********************************************************
kpss.test(mtof_total) # --> not-stationary
adf.test(mtof_total)  # --> stationary         ----> We need differencing 

### Differencing *********************************************************
ndiffs(mtof_total)  # --> lag-1 difference is needed
mtof_dff<-diff(mtof_total)
autoplot(mtof_dff)
kpss.test(mtof_dff) # --> stationary
adf.test(mtof_dff)  # --> stationary

## Total Data Forecasting & Hold-out Cross-validation --------------------
### Random Walk with Drift ***********************************************
MetOffice.RWD = rwf(mtof.train, h=mtof.split.size, drift = TRUE) 
autoplot(mtof.train) + 
  ggtitle("Met Office - Random Walk with Drift Forecast") +
  xlab("Year") +
  ylab("Temperature (K)") +
  autolayer(MetOffice.RWD, series="Random Walk \nwith Drift") +
  autolayer(mtof_total, series="Actual")+ 
  scale_colour_manual(values=c( "red", "black"),breaks=c("Random Walk \nwith Drift","Actual")) 

autoplot(mtof.train) + 
  ggtitle("Met Office - Random Walk with Drift Forecast") +
  xlab("Year") +
  ylab("Temperature (K)") +
  autolayer(MetOffice.RWD$mean, series="Random Walk \nwith Drift") +
  autolayer(mtof_total, series="Actual")+ 
  scale_colour_manual(values=c( "red", "black"),breaks=c("Random Walk \nwith Drift","Actual")) 

checkresiduals(MetOffice.RWD)  
summary(residuals(MetOffice.RWD))
accuracy(MetOffice.RWD, mtof.test)

### ETS: ANN ************************************************************
f_mtof <- mtof.split.size + 973 # Forecast until 2100
options(max.print=1000000)
MTOF_ETS_ANN <- ets(diff(mtof.train), model="ANN") 
MTOF_ETS_ANN_pred <- forecast(MTOF_ETS_ANN, h=f_mtof, level=c(0.8, 0.9)) 
MTOF_ETS_ANN
autoplot(MTOF_ETS_ANN_pred, xlab="Year", ylab="Temperature (K)") + 
  autolayer(diff(mtof.test))

### ETS: AAN ************************************************************
MTOF_ETS_AAN <- ets(mtof.train, model="AAN") 
MTOF_ETS_AAN_pred <- forecast(MTOF_ETS_AAN, h=f_mtof, level=c(0.8, 0.9)) 
MTOF_ETS_AAN
autoplot(mtof.train, series='Actual') + 
  autolayer(mtof.test, series='Actual') +
  autolayer(MTOF_ETS_AAN_pred$mean, series = 'AAN') + 
  scale_colour_manual(values=c('black','red'), breaks=c('Actual','AAN'))+
  ggtitle("Met Office, AAN")+
  xlab("Year") + ylab("Temperature (K)") + 
  theme(legend.position = "bottom", legend.title = element_blank())+
  coord_cartesian(xlim=c(1960,2030))
accuracy(MTOF_ETS_AAN_pred, mtof.test)

### ETS: AAA ************************************************************
MTOF_ETS_AAA <- ets(mtof.train, model="AAA") 
MTOF_ETS_AAA_pred <- forecast(MTOF_ETS_AAA, h=f_mtof, level=c(0.8, 0.9)) 
MTOF_ETS_AAA
autoplot(mtof.train, series='Actual') + 
  autolayer(MTOF_ETS_AAA_pred$mean, series = 'AAA') + 
  autolayer(mtof.test, series='Actual') +
  scale_colour_manual(values=c('black','red'), breaks=c('Actual','AAA'))+
  ggtitle("Met Office, ETS AAA")+
  xlab("Year") + ylab("Temperature (K)") + 
  theme(legend.position = "bottom", legend.title = element_blank())+
  coord_cartesian(xlim=c(1960,2030))  

autoplot(mtof.train, series='Actual') + 
  autolayer(MTOF_ETS_AAA_pred, series = 'AAA') + 
  autolayer(mtof.test, series='Actual') +
  scale_colour_manual(values=c('black','red'), breaks=c('Actual','AAA'))+
  ggtitle("Met Office, ETS AAA")+
  xlab("Year") + ylab("Temperature (K)") + 
  theme(legend.position = "bottom", legend.title = element_blank())+
  coord_cartesian(xlim=c(1850,2100),ylim=c(285,291))  

accuracy(MTOF_ETS_AAA_pred, mtof.test)
checkresiduals(MTOF_ETS_AAA)  
summary(residuals(MTOF_ETS_AAA))
kpss.test(residuals(MTOF_ETS_AAA)) # --> stationary
adf.test(residuals(MTOF_ETS_AAA))  # --> stationary    

### ETS: MAA ************************************************************
MTOF_ETS_MAA <- ets(mtof.train, model="MAA") 
MTOF_ETS_MAA_pred <- forecast(MTOF_ETS_MAA, h=f_mtof, level=c(0.8, 0.9)) 
MTOF_ETS_MAA
autoplot(mtof.train, series='Actual') + 
  autolayer(mtof.test, series='Actual') +
  autolayer(MTOF_ETS_MAA_pred$mean, series = 'MAA') + 
  scale_colour_manual(values=c('black','red'), breaks=c('Actual','MAA'))+
  ggtitle("Met Office, MAA")+
  xlab("Year") + ylab("Temperature (K)") + 
  theme(legend.position = "bottom", legend.title = element_blank())+
  coord_cartesian(xlim=c(1960,2030))  
accuracy(MTOF_ETS_MAA_pred, mtof.test)

### ETS: MAN ************************************************************
MTOF_ETS_MAN <- ets(mtof.train, model="MAN") 
MTOF_ETS_MAN_pred <- forecast(MTOF_ETS_MAN, h=f_mtof, level=c(0.8, 0.9)) 
MTOF_ETS_MAN
autoplot(mtof.train, series='Actual') + 
  autolayer(mtof.test, series='Actual') +
  autolayer(MTOF_ETS_MAN_pred$mean, series = 'MAN') + 
  scale_colour_manual(values=c('black','red'), breaks=c('Actual','MAN'))+
  ggtitle("Met Office, MAN")+
  xlab("Year") + ylab("Temperature (K)") + 
  theme(legend.position = "bottom", legend.title = element_blank())+
  coord_cartesian(xlim=c(1960,2030))  
accuracy(MTOF_ETS_MAN_pred, mtof.test)

### ETS: MAM ************************************************************
MTOF_ETS_MAM <- ets(mtof.train, model="MAM") 
MTOF_ETS_MAM_pred <- forecast(MTOF_ETS_MAM, h=f_mtof, level=c(0.8, 0.9)) 
MTOF_ETS_MAM
autoplot(mtof.train, series='Actual') + 
  autolayer(mtof.test, series='Actual') +
  autolayer(MTOF_ETS_MAM_pred$mean, series = 'MAM') + 
  scale_colour_manual(values=c('black','red'), breaks=c('Actual','MAM'))+
  ggtitle("Met Office, MAM")+
  xlab("Year") + ylab("Temperature (K)") + 
  theme(legend.position = "bottom", legend.title = element_blank())+
  coord_cartesian(xlim=c(1960,2030))  
accuracy(MTOF_ETS_MAM_pred, mtof.test)

### ETS: MMM ************************************************************
MTOF_ETS_MMM <- ets(mtof.train, model="MMM") 
MTOF_ETS_MMM_pred <- forecast(MTOF_ETS_MMM, h=f_mtof, level=c(0.8, 0.9)) 
MTOF_ETS_MMM
autoplot(mtof.train, series='Actual') + 
  autolayer(mtof.test, series='Actual') +
  autolayer(MTOF_ETS_MMM_pred$mean, series = 'MMM') + 
  scale_colour_manual(values=c('black','red'), breaks=c('Actual','MMM'))+
  ggtitle("Met Office, MMM")+
  xlab("Year") + ylab("Temperature (K)") + 
  theme(legend.position = "bottom", legend.title = element_blank())+
  coord_cartesian(xlim=c(1960,2030))  
accuracy(MTOF_ETS_MMM_pred, mtof.test)

### ETS: MNA ************************************************************
MTOF_ETS_MNA <- ets(mtof.train, model="MNA") 
MTOF_ETS_MNA_pred <- forecast(MTOF_ETS_MNA, h=f_mtof, level=c(0.8, 0.9)) 
MTOF_ETS_MNA
autoplot(mtof.train, series='Actual') + 
  autolayer(mtof.test, series='Actual') +
  autolayer(MTOF_ETS_MNA_pred$mean, series = 'MNA') + 
  scale_colour_manual(values=c('black','red'), breaks=c('Actual','MNA'))+
  ggtitle("Met Office, MNA")+
  xlab("Year") + ylab("Temperature (K)") + 
  theme(legend.position = "bottom", legend.title = element_blank())+
  coord_cartesian(xlim=c(1960,2030))  
accuracy(MTOF_ETS_MNA_pred, mtof.test)

### TBATS ***************************************************************
MTOF_TBATS <- tbats(mtof.train)
MTOF_TBATS_pred <- forecast(MTOF_TBATS, h=f_mtof, level=c(0.8, 0.9)) 
MTOF_TBATS
autoplot(mtof.train, series='Actual') + 
  autolayer(mtof.test, series='Actual') +
  autolayer(MTOF_TBATS_pred$mean, series = 'TBATS') + 
  scale_colour_manual(values=c('black','red'), breaks=c('Actual','TBATS'))+
  ggtitle("Met Office, TBATS")+
  xlab("Year") + ylab("Temperature (K)") + 
  theme(legend.position = "bottom", legend.title = element_blank())+
  coord_cartesian(xlim=c(1960,2030))  
accuracy(MTOF_TBATS_pred, mtof.test)

### (S)ARIMA Models ++++++++++++++++++++++++++++++++++++++++++++++++++++++
#### Exploratory Data Analysis: Need for Transformation
autoplot(mtof_total, series = 'Met Office') + 
  scale_colour_manual(values = c('black'),breaks=('Met Office'))
autoplot(log(mtof_total), series = 'Log Met Office') + 
  scale_colour_manual(values = c('black'),breaks=('Log Met Office'))

#### Decomposition
BoxCox.lambda(mtof_total) #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
SWIN <- as.integer(frequency(mtof_total)/2)*2+1 
TWIN <- 960
mtof_fitStl <-stl(log(mtof_total) ,s.window=SWIN, t.window=TWIN, robust=TRUE)
autoplot(mtof_fitStl)

#### Differencing & Stationarity
ndiffs(log(mtof_total)) # --> lag differencing: 1
nsdiffs(log(mtof_total)) # --> seasonal differencing: 0 

kpss.test(diff(log(mtof_total))) # --> Stationary
adf.test(diff(log(mtof_total))) # --> Stationary

#### Selecting p, q, P & Q
ggtsdisplay(diff(log(mtof_total)), lag=60)
    ## Non-seasonal behaviour
    ## --> ACF: MA(q=1),MA(q=2)
    ## --> PACF: AR(p=1),AR(p=2),AR(p=3)
    ## Seasonal behaviour
    ## --> ACF: MA(Q=0).MA(Q=1)
    ## --> PACF: AR(P=2),AR(P=3)
    ##
    ## Models to try
    ## 1  ARIMA(1,1,1)(2,0,0)12
    ## 2  ARIMA(2,1,1)(2,0,0)12
    ## 3  ARIMA(3,1,1)(2,0,0)12
    ## 4  ARIMA(1,1,1)(2,0,1)12
    ## 5  ARIMA(2,1,1)(2,0,1)12
    ## 6  ARIMA(3,1,1)(2,0,1)12
    ## 7  ARIMA(1,1,1)(3,0,0)12
    ## 8  ARIMA(2,1,1)(3,0,0)12
    ## 9  ARIMA(3,1,1)(3,0,0)12
    ## 10 ARIMA(1,1,1)(3,0,1)12
    ## 11 ARIMA(2,1,1)(3,0,1)12
    ## 12 ARIMA(3,1,1)(3,0,1)12

### Auto ARIMA ***********************************************************
MTOF_ARIMAA <- auto.arima(mtof.train, lambda = NULL , approximation = F, 
                          stepwise = F, trace = T)
MTOF_ARIMAA # Auto ARIMA <- ARIMA(1,1,2)(2,0,0)[12] with drift 
MTOF_ARIMAA_pred = forecast(MTOF_ARIMAA, h = f_mtof)
autoplot(mtof.train, series='Actual') + 
  autolayer(MTOF_ARIMAA_pred, series = 'ARIMAA') + 
  scale_colour_manual(values=c('black','red'), breaks=c('Actual','ARIMAA'))+
  autolayer(mtof.test, series='Actual') +
  ggtitle("Met Office, Auto-ARIMA (1,1,2)(2,0,0)[12] with drift ")+
  xlab("Year") + ylab("Temperature (K)") + 
  theme(legend.position = "bottom", legend.title = element_blank())+
  coord_cartesian(xlim=c(1960,2100),ylim=c(285,291))  
checkresiduals(MTOF_ARIMAA_pred)
accuracy(MTOF_ARIMAA_pred, mtof.test)

### ARIMA Alternatives ****************************************************
## 1  ARIMA(1,1,1)(2,0,0)12
## 2  ARIMA(2,1,1)(2,0,0)12
## 3  ARIMA(3,1,1)(2,0,0)12
## 4  ARIMA(1,1,1)(2,0,1)12
## 5  ARIMA(2,1,1)(2,0,1)12
## 6  ARIMA(3,1,1)(2,0,1)12
## 7  ARIMA(1,1,1)(3,0,0)12
## 8  ARIMA(2,1,1)(3,0,0)12
## 9  ARIMA(3,1,1)(3,0,0)12
## 10 ARIMA(1,1,1)(3,0,1)12
## 11 ARIMA(2,1,1)(3,0,1)12
## 12 ARIMA(3,1,1)(3,0,1)12
#### ARIMA1 ---------------------------------------------------------------
MTOF_ARIMA1<-Arima(mtof.train,order=c(1,1,1),seasonal=c(2,0,0),lambda=0,include.drift=TRUE)
MTOF_ARIMA1
MTOF_ARIMA1_pred = forecast(MTOF_ARIMA1, h = f_mtof)
checkresiduals(MTOF_ARIMA1_pred)
accuracy(MTOF_ARIMA1_pred, mtof.test)
#### ARIMA2 ---------------------------------------------------------------
MTOF_ARIMA2<-Arima(mtof.train,order=c(2,1,1),seasonal=c(2,0,0),lambda=0,include.drift=TRUE)
MTOF_ARIMA2
MTOF_ARIMA2_pred = forecast(MTOF_ARIMA2, h = f_mtof)
checkresiduals(MTOF_ARIMA2_pred)
accuracy(MTOF_ARIMA2_pred, mtof.test)
#### ARIMA3 ---------------------------------------------------------------
MTOF_ARIMA3<-Arima(mtof.train,order=c(3,1,1),seasonal=c(2,0,0),lambda=0,include.drift=TRUE)
MTOF_ARIMA3
MTOF_ARIMA3_pred = forecast(MTOF_ARIMA3, h = f_mtof)
checkresiduals(MTOF_ARIMA3_pred)
accuracy(MTOF_ARIMA3_pred, mtof.test)
#### ARIMA4 ---------------------------------------------------------------
MTOF_ARIMA4<-Arima(mtof.train,order=c(1,1,1),seasonal=c(2,0,1),lambda=0,include.drift=TRUE)
MTOF_ARIMA4
MTOF_ARIMA4_pred = forecast(MTOF_ARIMA4, h = f_mtof)
checkresiduals(MTOF_ARIMA4_pred)
accuracy(MTOF_ARIMA4_pred, mtof.test)
#### ARIMA5 ---------------------------------------------------------------
MTOF_ARIMA5<-Arima(mtof.train,order=c(2,1,1),seasonal=c(2,0,1),lambda=0,include.drift=TRUE)
MTOF_ARIMA5
MTOF_ARIMA5_pred = forecast(MTOF_ARIMA5, h = f_mtof)
checkresiduals(MTOF_ARIMA5_pred)
accuracy(MTOF_ARIMA5_pred, mtof.test)
#### ARIMA6 ---------------------------------------------------------------
MTOF_ARIMA6<-Arima(mtof.train,order=c(3,1,1),seasonal=c(2,0,1),lambda=0,include.drift=TRUE)
MTOF_ARIMA6
MTOF_ARIMA6_pred = forecast(MTOF_ARIMA6, h = f_mtof)
checkresiduals(MTOF_ARIMA6_pred)
accuracy(MTOF_ARIMA6_pred, mtof.test)
#### ARIMA7 ---------------------------------------------------------------
MTOF_ARIMA7<-Arima(mtof.train,order=c(1,1,1),seasonal=c(3,0,0),lambda=0,include.drift=TRUE)
MTOF_ARIMA7
MTOF_ARIMA7_pred = forecast(MTOF_ARIMA7, h = f_mtof)
checkresiduals(MTOF_ARIMA7_pred)
accuracy(MTOF_ARIMA7_pred, mtof.test)
#### ARIMA8 ---------------------------------------------------------------
MTOF_ARIMA8<-Arima(mtof.train,order=c(2,1,1),seasonal=c(3,0,0),lambda=0,include.drift=TRUE)
MTOF_ARIMA8
MTOF_ARIMA8_pred = forecast(MTOF_ARIMA8, h = f_mtof)
checkresiduals(MTOF_ARIMA8_pred)
accuracy(MTOF_ARIMA8_pred, mtof.test)
#### ARIMA9 ---------------------------------------------------------------
MTOF_ARIMA9<-Arima(mtof.train,order=c(3,1,1),seasonal=c(3,0,0),lambda=0,include.drift=TRUE)
MTOF_ARIMA9
MTOF_ARIMA9_pred = forecast(MTOF_ARIMA9, h = f_mtof)
checkresiduals(MTOF_ARIMA9_pred)
accuracy(MTOF_ARIMA9_pred, mtof.test)
#### ARIMA10 ---------------------------------------------------------------
MTOF_ARIMA10<-Arima(mtof.train,order=c(1,1,1),seasonal=c(3,0,0),lambda=0,include.drift=TRUE)
MTOF_ARIMA10
MTOF_ARIMA10_pred = forecast(MTOF_ARIMA10, h = f_mtof)
checkresiduals(MTOF_ARIMA10_pred)
accuracy(MTOF_ARIMA10_pred, mtof.test)
#### ARIMA11 ---------------------------------------------------------------
MTOF_ARIMA11<-Arima(mtof.train,order=c(2,1,1),seasonal=c(3,0,0),lambda=0,include.drift=TRUE)
MTOF_ARIMA11
MTOF_ARIMA11_pred = forecast(MTOF_ARIMA11, h = f_mtof)
checkresiduals(MTOF_ARIMA11_pred)
accuracy(MTOF_ARIMA11_pred, mtof.test)
#### ARIMA12 ---------------------------------------------------------------
MTOF_ARIMA12<-Arima(mtof.train,order=c(3,1,1),seasonal=c(3,0,0),lambda=0,include.drift=TRUE)
MTOF_ARIMA12
MTOF_ARIMA12_pred = forecast(MTOF_ARIMA12, h = f_mtof)
checkresiduals(MTOF_ARIMA12_pred)
accuracy(MTOF_ARIMA12_pred, mtof.test)
#### Visualize the best ARIMA models
#### Auto-SARIMA, SARIMA4, SARIMA5, SARIMA6
autoplot(mtof.train, series='Actual') +
  autolayer(MTOF_ARIMAA_pred$mean, series = 'AutoARIMA') +
  autolayer(MTOF_ARIMA4_pred$mean, series = 'ARIMA4') +
  autolayer(MTOF_ARIMA5_pred$mean, series = 'ARIMA5') +
  autolayer(MTOF_ARIMA6_pred$mean, series = 'ARIMA6') +
  autolayer(mtof.test, series='Actual')+
  ggtitle('Best ARIMA Models')+
  xlab('Year')+ylab('Temperature (K)')+
  scale_colour_manual(values=c('black','red','blue','green','orange'),
                      breaks = c('Actual','AutoARIMA','ARIMA4','ARIMA5','ARIMA6'))+
  coord_cartesian(xlim = c(1960,2030))
#### Auto-ARIMA, ARIMA4  
autoplot(mtof.train, series='Actual') +
  autolayer(MTOF_ARIMAA_pred$mean, series = 'AutoARIMA') +
  autolayer(MTOF_ARIMA4_pred$mean, series = 'ARIMA4') +
  autolayer(mtof.test, series='Actual')+
  ggtitle('Best ARIMA Models')+
  xlab('Year')+ylab('Temperature (K)')+
  scale_colour_manual(values=c('black','red','blue'),
                      breaks = c('Actual','AutoARIMA','ARIMA4'))+
  coord_cartesian(xlim = c(1900,2100))  

### STLF ***************************************************************
MTOF_STLF <- stlf(mtof.train, h=f_mtof , s.window=SWIN,
                  t.window = TWIN, lambda=0,level=c(0.8, 0.9))
autoplot(mtof.train, series='Actual') + 
  autolayer(mtof.test, series='Actual') +
  autolayer(MTOF_STLF$mean, series = 'STLF') + 
  scale_colour_manual(values=c('black','red'), breaks=c('Actual','STLF'))+
  ggtitle("Met Office, STLF")+
  xlab("Year") + ylab("Temperature (K)") + 
  theme(legend.position = "bottom", legend.title = element_blank())+
  coord_cartesian(xlim=c(1960,2100))  
accuracy(MTOF_STLF, mtof.test)

## Ensemble with the best models ----------------------------------------------------------------
## The best models: Random Walk with Drift, ETS(AAA), Auto_SARIMA, SARIMA4
### Visulization - multi method
MetOffice.RWD = rwf(mtof.train, h=f_mtof, drift = TRUE) 

autoplot(mtof.train, series='Actual') +
  autolayer(MTOF_ARIMAA_pred$mean, series = 'Auto-SARIMA') +
  autolayer(MTOF_ARIMA4_pred$mean, series = 'SARIMA4') +
  autolayer(MetOffice.RWD$mean, series = 'Random-Walk-wDrift') +
  autolayer(MTOF_ETS_AAA_pred$mean, series = 'ETS(AAA)') +
  autolayer(MTOF_STLF$mean, series = 'STLF') +
  autolayer(mtof.test, series='Actual')+
  ggtitle('Met Office: Best Forecasting Models')+
  xlab('Year')+ylab('Temperature (K)')+
  scale_colour_manual(values=c('black','red','blue','green','grey','cyan'),
                      breaks = c('Actual','Auto-SARIMA','SARIMA4',
                                 'Random-Walk-wDrift','ETS(AAA)','STLF' ))+
  coord_cartesian(xlim = c(1900,2100))+
  theme(legend.title = element_blank(), legend.position = 'bottom')

### Tabulated Results: Best Models
ARIMA4_acc <- accuracy(MTOF_ARIMA4_pred, mtof.test)
ARIMA5_acc <- accuracy(MTOF_ARIMA5_pred, mtof.test)
ARIMA6_acc <- accuracy(MTOF_ARIMA6_pred, mtof.test)
ARIMA10_acc <- accuracy(MTOF_ARIMA10_pred, mtof.test)

ARIMA4_acc <- ARIMA4_acc[-1,]
ARIMA5_acc <- ARIMA5_acc[-1,]
ARIMA6_acc <- ARIMA6_acc[-1,]
ARIMA10_acc <- ARIMA10_acc[-1,]

multi_acc <- rbind(ARIMA4_acc,ARIMA5_acc,ARIMA6_acc,ARIMA10_acc)
multi_acc

### Tabulated Results: All Models
ARIMAA_acc <- accuracy(MTOF_ARIMAA_pred, mtof.test)
ARIMA1_acc <- accuracy(MTOF_ARIMA1_pred, mtof.test)
ARIMA2_acc <- accuracy(MTOF_ARIMA2_pred, mtof.test)
ARIMA3_acc <- accuracy(MTOF_ARIMA3_pred, mtof.test)
ARIMA4_acc <- accuracy(MTOF_ARIMA4_pred, mtof.test)
ARIMA5_acc <- accuracy(MTOF_ARIMA5_pred, mtof.test)
ARIMA6_acc <- accuracy(MTOF_ARIMA6_pred, mtof.test)
ARIMA7_acc <- accuracy(MTOF_ARIMA7_pred, mtof.test)
ARIMA8_acc <- accuracy(MTOF_ARIMA8_pred, mtof.test)
ARIMA9_acc <- accuracy(MTOF_ARIMA9_pred, mtof.test)
ARIMA10_acc <- accuracy(MTOF_ARIMA10_pred, mtof.test)
ARIMA11_acc <- accuracy(MTOF_ARIMA11_pred, mtof.test)
ARIMA12_acc <- accuracy(MTOF_ARIMA12_pred, mtof.test)
RWF_acc <- accuracy(MetOffice.RWD, mtof.test)
ETS_ANN_acc <- accuracy(MTOF_ETS_ANN_pred, mtof.test)
ETS_AAN_acc <- accuracy(MTOF_ETS_AAN_pred, mtof.test)
ETS_AAA_acc <- accuracy(MTOF_ETS_AAA_pred, mtof.test)
ETS_MAN_acc <- accuracy(MTOF_ETS_MAN_pred, mtof.test)
ETS_MAA_acc <- accuracy(MTOF_ETS_MAA_pred, mtof.test)
ETS_MAM_acc <- accuracy(MTOF_ETS_MAM_pred, mtof.test)
ETS_MMM_acc <- accuracy(MTOF_ETS_MMM_pred, mtof.test)
ETS_MNA_acc <- accuracy(MTOF_ETS_MNA_pred, mtof.test)
STLF_acc <- accuracy(MTOF_STLF, mtof.test)
TBATS_acc <- accuracy(MTOF_TBATS_pred, mtof.test)

ARIMAA_acc <- ARIMAA_acc[-1,]
ARIMA1_acc <- ARIMA1_acc[-1,]
ARIMA2_acc <- ARIMA2_acc[-1,]
ARIMA3_acc <- ARIMA3_acc[-1,]
ARIMA4_acc <- ARIMA4_acc[-1,]
ARIMA5_acc <- ARIMA5_acc[-1,]
ARIMA6_acc <- ARIMA6_acc[-1,]
ARIMA7_acc <- ARIMA7_acc[-1,]
ARIMA8_acc <- ARIMA8_acc[-1,]
ARIMA9_acc <- ARIMA9_acc[-1,]
ARIMA10_acc <- ARIMA10_acc[-1,]
ARIMA11_acc <- ARIMA11_acc[-1,]
ARIMA12_acc <- ARIMA12_acc[-1,]
RWF_acc <- RWF_acc[-1,]
ETS_ANN_acc <- ETS_ANN_acc[-1,]
ETS_AAN_acc <- ETS_AAN_acc[-1,]
ETS_AAA_acc <- ETS_AAA_acc[-1,]
ETS_MAN_acc <- ETS_MAN_acc[-1,]
ETS_MAA_acc <- ETS_MAA_acc[-1,]
ETS_MAM_acc <- ETS_MAM_acc[-1,]
ETS_MMM_acc <- ETS_MMM_acc[-1,]
ETS_MNA_acc <- ETS_MNA_acc[-1,]
STLF_acc <- STLF_acc[-1,]

all_acc <- rbind(RWF_acc,ETS_ANN_acc,ETS_AAN_acc,ETS_AAA_acc,
                   ETS_MAN_acc,ETS_MAA_acc,ETS_MAM_acc,ETS_MMM_acc,
                   ETS_MNA_acc,STLF_acc,ARIMAA_acc,ARIMA1_acc,
                   ARIMA2_acc,ARIMA3_acc,ARIMA4_acc,ARIMA5_acc,
                   ARIMA6_acc,ARIMA7_acc,ARIMA8_acc,ARIMA9_acc,
                   ARIMA10_acc,ARIMA11_acc,ARIMA12_acc)
all_acc

### Forecast
Ensmbl_frcst <- (MTOF_ARIMA4_pred$mean +
                   MTOF_ARIMA5_pred$mean +
                   MTOF_ARIMA6_pred$mean +
                   MTOF_ARIMA10_pred$mean)/4

Ensmbl_acc <- accuracy(Ensmbl_frcst, mtof.test)
Ensmbl_acc

### Visulization - Ensemble
autoplot(mtof.train, series='Actual') +
  autolayer(MTOF_ARIMA4_pred$mean, series = 'SARIMA4') +
  autolayer(MTOF_ARIMA5_pred$mean, series = 'SARIMA5') +
  autolayer(MTOF_ARIMA6_pred$mean, series = 'SARIMA6') +
  autolayer(MTOF_ARIMA10_pred$mean, series = 'SARIMA10') +
  autolayer(Ensmbl_frcst, series = 'Ensemble') +
  autolayer(mtof.test, series='Actual')+
  ggtitle('Met Office: Best Forecasting Models')+
  xlab('Year')+ylab('Temperature (K)')+
  scale_colour_manual(values=c('black','grey','grey','grey','grey','red'),
                      breaks = c('Actual','SARIMA4','SARIMA5',
                                 'SARIMA6','SARIMA10','Ensemble' ))+
  coord_cartesian(xlim = c(1840,2100),ylim = c(286,289))+
  theme(legend.title = element_blank(), legend.position = 'right')

### Best Model Uncertainty
MTOF_ARIMA5_pred = forecast(MTOF_ARIMA5, h = f_mtof, level = c(0.5,0.7,0.9))
autoplot(mtof.train, series='Actual') +
  autolayer(MTOF_ARIMA5_pred, series = 'SARIMA5') +
  autolayer(mtof.test, series='Actual')+
  ggtitle('Met Office: Point Forecast and Confidence Intervals')+
  xlab('Year')+ylab('Temperature (K)')+
  scale_colour_manual(values=c('black','grey','blue','grey','grey','red'),
                      breaks = c('Actual','Auto-SARIMA','SARIMA5',
                                 'Random-Walk-wDrift','ETS(AAA)','Ensemble' ))+
  coord_cartesian(xlim = c(1840,2100))+
  theme(legend.title = element_blank(), legend.position = 'right')

### Write out the results
write.csv(Ensmbl_frcst, 'MetOffice_ensemble.csv')
write.csv(MTOF_ARIMA5_pred, 'MetOffice_bestmodel(ARIMA).csv')

##################################################
## Questions 6 & 7
##################################################

## Create Train & Test data

# Per the climate bet, we'll train our model on data up to and
# including 2007 and make predictions on the 10 year period of
# 2008 to 2017 (inclusive).

mtof_total67 <- ts(df_mtof$ATemp,start = c(1850,1), 
                   end = c(2017,12), frequency=12)

mtof.train67 <- window(mtof_total67, start = c(1850,1),end = c(2007,10))

mtof.test67 <- window(mtof_total67, start = c(2008,1),end = c(2017,12))

#############################
## Armstrong's Model - Naive
#############################

# NAIVE model using training data to fit model, h=prediction horizon
ARM_Naive = naive(mtof.train67, h=120) 

autoplot(mtof.train67) + ggtitle("Naive forecast for Met Office Temperatures 10 Yrs") +
  autolayer(ARM_Naive$mean, series="Naive") +  #adding prediction line with 80% and 95% prediction interval
  xlab("Year") +
  autolayer(mtof_total67, series="Actual")+ 
  ylab("Temperature (K)") +
  scale_colour_manual(values=c("red","black"),
                      breaks=c("Naive", "Actual"))+
  coord_cartesian(xlim = c(2006,2017))+
  theme(legend.title = element_blank(), 
        legend.position = 'bottom')

checkresiduals(ARM_Naive)
accuracy(ARM_Naive, mtof.test67)

#write.csv(ARM_Naive$mean, 'MetOff_Naive_2008-2017.csv') 
#write.csv(mtof_total67, 'MetOff_Actual_1850-2017.csv') 

##################################
## Our Best Model - RWD with Drift
##################################

### Random Walk with Drift ***********************************************
MetOffice.RWD67 = rwf(mtof.train67, h=120, drift = TRUE) 

autoplot(mtof.train67) + 
  ggtitle("Met Office - Random Walk with Drift Forecast 10 Yrs") +
  xlab("Year") +
  ylab("Temperature (K)") +
  autolayer(MetOffice.RWD67, series="Random Walk \nwith Drift") +
  autolayer(mtof_total67, series="Actual")+ 
  scale_colour_manual(values=c( "red", "black"),
                      breaks=c("Random Walk \nwith Drift","Actual")) 

autoplot(mtof.train67) + 
  ggtitle("Met Office - Random Walk with Drift Forecast 10 Yrs") +
  xlab("Year") +
  ylab("Temperature (K)") +
  autolayer(MetOffice.RWD67$mean, series="Random Walk \nwith Drift") +
  autolayer(mtof_total67, series="Actual")+ 
  scale_colour_manual(values=c( "red", "black"),
                      breaks=c("Random Walk \nwith Drift","Actual")) 

checkresiduals(MetOffice.RWD67)   
summary(residuals(MetOffice.RWD67))
accuracy(MetOffice.RWD67, mtof.test67) 

#### Visualization of the model with the best RMSE and MAE on the testing set --> Random Walk 
#### together with Armstrong's Naive model
autoplot(mtof_total67, series='Actual') +
  autolayer(ARM_Naive$mean, series="Armstrong's Naive") +
  autolayer(MetOffice.RWD67$mean, series="Random Walk \nwith Drift") +
  autolayer(mtof.test67, series='Actual')+
  ggtitle("Climate Bet 10 Yrs: \nArmstrong's Naive Model vs. Random Walk with Drift")+
  xlab('Year')+ylab('Temperature (K)')+
  scale_colour_manual(values=c('black','blue','red'),
                      breaks = c('Actual',"Armstrong's Naive",'Random Walk \nwith Drift'))+
  coord_cartesian(xlim = c(2005,2020),ylim=c(287,288.5))+
  theme(legend.title = element_blank(), 
        legend.position = 'bottom')



###########################################################
## Q7 - Repeat Predictions for a 30 year period (1988-2017)
###########################################################


## Recreate Train & Test data

mtof_total67_1 <- ts(df_mtof$ATemp,start = c(1850,1), 
                     end = c(2017,12), frequency=12)

mtof.train67_1 <- window(mtof_total67_1, start = c(1850,1),end = c(1987,11))

mtof.test67_1 <- window(mtof_total67_1, start = c(1988,1),end = c(2017,12))


#############################
## Armstrong's Model - Naive
#############################

ARM_Naive_1 = naive(mtof.train67_1, h=360) 

autoplot(mtof.train67_1) + ggtitle("Naive forecast for Met Office Temperatures 30 Yrs") +
  autolayer(ARM_Naive_1$mean, series="Naive") + 
  xlab("Year") +
  autolayer(mtof_total67_1, series="Actual")+ 
  ylab("Temperature (K)") +
  scale_colour_manual(values=c("red","black"),
                      breaks=c("Naive", "Actual"))+
  coord_cartesian(xlim = c(1986,2017))+
  theme(legend.title = element_blank(), 
        legend.position = 'bottom')

checkresiduals(ARM_Naive_1) 
accuracy(ARM_Naive_1$mean, mtof.test67_1)

##################################
## Our Best Model - RWD with Drift
##################################


MetOffice.RWD67_1 = rwf(mtof.train67_1, h=360, drift = TRUE) 

autoplot(mtof.train67_1) + 
  ggtitle("Met Office - Random Walk with Drift Forecast 30 Yrs") +
  xlab("Year") +
  ylab("Temperature (K)") +
  autolayer(MetOffice.RWD67_1, series="Random Walk \nwith Drift") +
  autolayer(mtof_total67_1, series="Actual")+ 
  scale_colour_manual(values=c( "red", "black"),
                      breaks=c("Random Walk \nwith Drift","Actual")) 

autoplot(mtof.train67_1) + 
  ggtitle("Met Office - Random Walk with Drift Forecast 30 Yrs") +
  xlab("Year") +
  ylab("Temperature (K)") +
  autolayer(MetOffice.RWD67_1$mean, series="Random Walk \nwith Drift") +
  autolayer(mtof_total67_1, series="Actual")+ 
  scale_colour_manual(values=c( "red", "black"),
                      breaks=c("Random Walk \nwith Drift","Actual")) 

checkresiduals(MetOffice.RWD67_1)  
summary(residuals(MetOffice.RWD67_1))
accuracy(MetOffice.RWD67_1$mean, mtof.test67_1) 

##################################
## Our Best Model - SARIMA 
##################################
#### ARIMA4 ---------------------------------------------------------------
MTOF_BET_ARIMA4<-Arima(mtof.train67_1,order=c(1,1,1),seasonal=c(2,0,1),lambda=0,include.drift=TRUE)
MTOF_BET_ARIMA4
MTOF_BET_ARIMA4_pred = forecast(MTOF_BET_ARIMA4, h = 360)
checkresiduals(MTOF_BET_ARIMA4_pred)
accuracy(MTOF_BET_ARIMA4_pred, mtof.test67_1)



#### Visualization of the model with the best RMSE and MAE on the testing set --> Random Walk 
#### together with Armstrong's Naive model
autoplot(mtof_total67_1, series='Actual') +
  autolayer(ARM_Naive_1$mean, series="Armstrong's Naive") +
  autolayer(MetOffice.RWD67_1$mean, series="Random Walk \nwith Drift") +
  autolayer(mtof_total67_1, series='Actual')+
  ggtitle("Climate Bet: \nArmstrong's Naive Model vs. Random Walk with Drift 30Yrs")+
  xlab('Year')+ylab('Temperature (K)')+
  scale_colour_manual(values=c('black','blue','red'),
                      breaks = c('Actual',"Armstrong's Naive",'Random Walk \nwith Drift'))+
  coord_cartesian(xlim = c(1985,2020))+
  theme(legend.title = element_blank(), legend.position = 'bottom')

#### Visualization of the model with the best RMSE and MAE on the testing set --> Random Walk 
#### together with Armstrong's Naive model
autoplot(mtof_total67_1, series='Actual') +
  autolayer(ARM_Naive_1$mean, series="Armstrong's Naive") +
  autolayer(MetOffice.RWD67_1$mean, series="Random Walk \nwith Drift") +
  autolayer(MTOF_BET_ARIMA4_pred$mean, series="SARIMA4") +
  autolayer(mtof_total67_1, series='Actual')+
  ggtitle("Climate Bet: \nArmstrong's Naive Model vs. Random Walk with Drift vs. SARIMA - 30Yrs")+
  xlab('Year')+ylab('Temperature (K)')+
  scale_colour_manual(values=c('black','blue','red','darkgreen'),
                      breaks = c('Actual',"Armstrong's Naive",
                                 'Random Walk \nwith Drift', 'SARIMA4'))+
  coord_cartesian(xlim = c(1985,2020))+
  theme(legend.title = element_blank(), legend.position = 'bottom')
