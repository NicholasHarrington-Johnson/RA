# Read in data sheet
load("truncated_data.Rda")
# tr[[i]] now refers to the data sheet for the i'th restaurant

# Pickup Method
pickup <- function(data,h){
  # Log+1 transformation (remove info beyond h in this step)
  log_data <- apply(data+1,2,log)[,1:(h+1)]
  # Split into training data and forecasting values
  top <- log_data[1:(nrow(log_data)-h),]
  # Calculate grossing up factors
  # Take the mean of each column and then take the difference of those means
  means <- apply(top,2,mean)
  devs <- rev(diff(rev(means)))
  # Values that are grossed up from when forecasting
  end <- tail(log_data,n=h)
  diag.end <- end*upper.tri(end,diag=FALSE)
  # Forecast
  # Extract diagonal values to forecast from
  forecast.v <- diag(end[,2:(h+1)])
  # Lower triangle matrix of gross up factors
  dev.mat<-t(matrix(rep(devs,(h)),nrow=h))
  dev.mat <- dev.mat*lower.tri(dev.mat,diag=TRUE)
  # Add column of zeros for grossing up purposes
  dev.mat.2 <- cbind(dev.mat,rep(0,h))
  # Add diagonal values to development factor matrix
  # Extract diagonal elements as matrix
  # Include 0 column
  beg.elem <- cbind(rep(0,h),diag(diag(end)))
  # Add
  tot <- dev.mat.2 + beg.elem
  # Cumulate gross up factors and forecasting values
  fin <- t(apply(tot,1,function(x){rev(cumsum(rev(x)))}))
  # Extract forecasted column
  fc <- fin[,1]
  # Exponentiate
  fc <- exp(fc)-1
  # Return forecast
  return(fc)
}
# Regression models with ARIMA Errors
booking.arima <- function(data,h,phols=0,full.fit=0){
  # Take training set
  new.data <- data[(1:(nrow(data)-h)),]
  # Apply time series properties
  new.data$pubd <- ts(new.data$pubd,end=tsp(data$pubd)[1]+((nrow(data)-h)/365))
  new.data$pubd <- ts(new.data$pubd,start=tsp(data$pubd)[1],frequency=tsp(data$pubd)[3])
  # Perform log transformation
  log.attend <- apply(new.data+1,2,log)[,1]
  # Fix frequency to 7
  log.attend <- ts(log.attend,start=1,frequency=7)
  # Remove zeros from data
  log.attend[log.attend<1] <- NA
  # Take Public Holiday Variables if required
  if(phols==1){
    xdums <- new.data[,(ncol(new.data)-2):ncol(new.data)]
    colnames(xdums) <- c("pubd","pubi","pubny")
    # Figure out if there are public holidays in sample
    total.xdum <- colSums(xdums)
    # Remove dummy variables if there are no public holidays in sample
    if (total.xdum[1]==0) {xdums$pubd <- NULL}
    if (total.xdum[2]==0) {xdums$pubi <- NULL}
    if (total.xdum[3]==0) {xdums$pubny <- NULL} 
    # Fit Regression Model with Public Holidays and ARIMA Errors
    fit <- auto.arima(log.attend,xreg=xdums)
  } else {
    # Fit ARIMA Model
    fit <- auto.arima(log.attend)
  }
  # Create x regressors as required
  if(phols==1){
    xfor <- data[((nrow(data)-h):nrow(data)),(ncol(data)-2):ncol(data)]
    colnames(xfor) <- c("pubd","pubi","pubny")
    # Remove x regressors as required
    if (total.xdum[1]==0) {xfor$pubd <- NULL} 
    if (total.xdum[2]==0) {xfor$pubi <- NULL}
    if (total.xdum[3]==0) {xfor$pubny <- NULL}
    # Forecast with Public Holiday Dummies
    fc <- forecast(fit,h,xreg=xfor)
  } else{
    # Forecast
    fc <- forecast(fit,h)
  }
  # Exponentiate
  fc$mean <- exp(fc$mean)-1
  fc$lower <- exp(fc$lower)-1
  fc$upper <- exp(fc$upper)-1
  fc$x <- ts(new.data[,1],frequency=365)
  tsp(fc$x)<-tsp(new.data$pubd)
  fc$mean <- ts(fc$mean, end =tsp(data$pubd)[2], frequency=365)
  fc$mean <- ts(fc$mean, start = (tsp(fc$x)[2]+(1/365)), frequency=365)
  tsp(fc$upper) <- tsp(fc$lower) <- tsp(fc$mean) 
  if(full.fit==1){
    obj <- list(fit,fc)
    return(obj)
  }
  return(fc)
}
# Next one with splines