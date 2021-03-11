library(ggplot2)
library(tseries)
options(stringsAsFactors = FALSE)
#setWD
setwd(getwd())

#Read data
AirportData <- read.csv("E:/HKU Study/Time series forcasting/Project/covid_impact_on_airport_traffic.csv")

#Pre-processing
airports<- unique(AirportData$AirportName)
sum(lengths(airports))
plotAllAirports <- function(data){
  airports<- unique(data$AirportName)
  for (airport in airports){
    Red <- subset(data, select=c('Date', 'PercentOfBaseline','AirportName'))
    anaData <- Red[Red[, "AirportName"] == airport,]
    ordData <- anaData[order(anaData$Date),]
    ordData$Date <- as.Date(ordData$Date)
    plot(PercentOfBaseline ~ Date, ordData, xaxt = "n", type = "l",main=airport)
  }
}
plotAllAirports(AirportData)
#Dallas/Fort Worth International is selected since the trend seems follow good
#property
Reduced <- subset(AirportData, select=c('Date', 'PercentOfBaseline','AirportName'))
anaData <- Reduced[Reduced[, "AirportName"] == "Dallas/Fort Worth International ",]
ordData <- subset(anaData, select=c('Date', 'PercentOfBaseline'))
ordData$Date <- as.Date(ordData$Date)
ordData <- ordData[order(ordData$Date),]
plot(PercentOfBaseline ~ Date, ordData, xaxt = "n", type = "l",main='Dallas/Fort Worth International ')
ordSeries <- ts(ordData$PercentOfBaseline)

#Split train-test Data
train <- head(ordSeries,-5)
test <- tail(ordSeries,5)

#According to plot, the time series may not be stationary since it has a weak upward
#trend. First plot ACF
acf(ordSeries,ci.type='ma')
pacf(ordSeries)
adf.test(ordSeries, alternative = "stationary")
#The plot also shows strong seasonality exists in the data set.

#Differencing
firstOrder <- diff(ordSeries , 1)
plot(firstOrder, xaxt = "n", type = "l",main="1-st order Differencing")
seventhOrder <- diff(ordSeries , 7)
plot(seventhOrder, xaxt = "n", type = "l",main="7-th order Differencing")


acf(firstOrder,ci.type='ma')
pacf(firstOrder)
adf.test(firstOrder, alternative = "stationary")

acf(seventhOrder,ci.type='ma')
pacf(seventhOrder)
adf.test(seventhOrder, alternative = "stationary")

diff71 <- diff(seventhOrder , 1)
diff77 <- diff(seventhOrder, 7)
acf(diff71)
pacf(diff71)
acf(diff77)
pacf(diff77)

diff11 <- diff(firstOrder,1)
diff17 <- diff(firstOrder,7)
acf(diff11)
pacf(diff11)
acf(diff17,ci.type='ma')
pacf(diff17)


plot(diff17, xaxt = "n", type = "l",main="After Differencing")

#Model 1
fitSample <- arima(train,order=c(0,1,1),
                   seasonal = (list ( order =c(1,1,0), period = 7)),method="ML")
fitSample
residualF <- fitSample$residuals
acf(fitSample$residuals, ci.type = 'ma')
pacf(fitSample$residuals)

BoxStat <- numeric(14)
for (i in 1:14){
  BoxStat[i] = Box.test(residualF, lag = i, type = "Ljung", fitdf = 16)
}

#Model 2
fitII <- arima(train,order=c(0,1,1),
                   seasonal = (list ( order =c(1,1,1), period = 7)),method="ML")
fitII
residualF2 <- fitII$residuals
acf(residualF2, ci.type = 'ma')
pacf(residualF2)

BoxStat2 <- numeric(14)
for (i in 1:14){
  BoxStat2[i] = Box.test(residualF2, lag = i, type = "Ljung", fitdf = 16)
}

#Model 3
fitIII <- arima(train,order=c(0,1,1),
               seasonal = (list ( order =c(0,1,1), period = 7)),method="ML")
fitIII
residualF3 <- fitIII$residuals
acf(residualF3, ci.type = 'ma')
pacf(residualF3)

fitIII$aic

BoxStat3 <- numeric(14)
for (i in 1:14){
  BoxStat3[i] = Box.test(residualF3, lag = i, type = "Ljung", fitdf = 16)
}

#BFS
#Generating all combinations of (p,d,q)
createChunk <- function(k) {
  chunks <- list()
  for (i in 0:k){
    ca = c(i,0,0)
    for (j in 0:k){
      cb = ca
      cb[2] = j
      for (m in 0:k){
        cc = cb
        cc[3] = m
        n<-length(chunks)
        chunks[[n+1]] = cc
      }
    }
  }
  return(chunks)
}

#BFS, printing failed combinations
fitModels <- function(k){
  chunk <- createChunk(k)
  chunk2 <- createChunk(k)
  totalN <- length(chunk)*length(chunk2)
  results <- list()
  AIC0 <- 999999
  bestARIMA <- c(0,0,0)
  bestSARIMA <- c(0,0,0)
  for (i in 1:length(chunk)){
    for (j in 1:length(chunk2)){
      tryCatch({
        fit <- arima(train,order=chunk[[i]],
                     seasonal = (list ( order =chunk2[[j]], period = 7)),method="ML")
        AIC <- fit$aic
        if (AIC <= AIC0) {
          AIC0=AIC
          bestARIMA = chunk[[i]]
          bestSARIMA = chunk2[[j]]
        }
      }, error=function(e){
        print(chunk[[i]],chunk2[[j]])
      })
    }
  }
  return(c(AIC0,bestARIMA,bestSARIMA))
}

getBest <- fitModels(2)

#Forecasting

predictI <- predict ( fitSample , n.ahead = 5)
predictII <- predict ( fitII , n.ahead = 5)
predictIII <- predict ( fitIII , n.ahead = 5)

errorI <- predictI$pred - test
errorII <- predictII$pred -test
errorIII <- predictIII$pred - test

sumOfSquare <- function(x) {
  length(x)*mean(x^2)
}

SSI <- sumOfSquare(errorI)
SSII <- sumOfSquare(errorII)
SSIII <- sumOfSquare(errorIII)
